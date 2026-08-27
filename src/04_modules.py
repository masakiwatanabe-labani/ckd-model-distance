"""モジュールレベル解析。

要点:
  - 手作りモジュールと MSigDB の両方で回し、結論が両者で一致するかを確認する。
  - **モジュール性の検定を必ず通す**。探索段階で「神経系モジュール」に見えた集合は、
    遺伝子間共発現がランダム集合と区別できなかった（平均r -0.009 vs +0.007±0.057）。
    共発現構造を持たない集合は「モジュール」として解釈してはいけない。
  - 競合的 z と、検体単位の固有値（eigengene）の両方を出す。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, competitive_z, append_summary  # noqa: E402

log = get_logger("04_modules")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
REF = ROOT / cfg["paths"]["ref"]
SEED = cfg["seed"]

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))


def load_modules() -> dict[str, list[str]]:
    with open(ROOT / "config" / "modules.yaml", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return {k: str(v).split() for k, v in raw.items()}


def load_gmt(path: Path, min_n=10, max_n=500) -> dict[str, list[str]]:
    sets = {}
    if not path.exists():
        return sets
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        genes = [g.strip().upper() for g in parts[2:] if g.strip()]
        if min_n <= len(genes) <= max_n:
            sets[parts[0]] = genes
    return sets


def lfc(dataset: str, to_mouse=False) -> pd.Series:
    s = DE.loc[dataset, "lfc"]
    s = s[np.isfinite(s)]
    if to_mouse:
        idx = pd.Series(s.index).map(CAT2MOUSE)
        s = s.set_axis(idx)
        s = s[s.index.notna()]
    s.index = [str(i).upper() for i in s.index]
    return s[~s.index.duplicated()]


def module_coherence(mat: pd.DataFrame, genes: list[str], n_null=200, seed=0):
    """集合内の平均相関を、発現量を合わせたランダム集合と比較する。"""
    m = mat.copy()
    m.index = [str(i).upper() for i in m.index]
    m = m[~m.index.duplicated()]
    g = [x for x in genes if x in m.index]
    sub = m.loc[g]
    sub = sub[sub.std(axis=1) > 0]
    if len(sub) < 4:
        return np.nan, np.nan, np.nan, len(sub)
    c = np.corrcoef(sub.to_numpy())
    iu = np.triu_indices_from(c, 1)
    obs = float(np.nanmean(c[iu]))
    # 発現量パーセンタイルを合わせたランダム集合
    mean_expr = m.mean(axis=1)
    target = mean_expr.loc[sub.index].rank(pct=True)
    rng = np.random.default_rng(seed)
    pool_rank = mean_expr.rank(pct=True)
    null = []
    for _ in range(n_null):
        pick = []
        for t in target:
            cand = pool_rank[(pool_rank - t).abs() < 0.05].index
            pick.append(rng.choice(cand) if len(cand) else rng.choice(pool_rank.index))
        s2 = m.loc[list(dict.fromkeys(pick))]
        s2 = s2[s2.std(axis=1) > 0]
        if len(s2) < 4:
            continue
        cc = np.corrcoef(s2.to_numpy())
        null.append(float(np.nanmean(cc[np.triu_indices_from(cc, 1)])))
    null = np.array(null)
    z = (obs - null.mean()) / null.std() if null.size > 5 else np.nan
    p = (np.sum(null >= obs) + 1) / (null.size + 1) if null.size > 5 else np.nan
    return obs, z, p, len(sub)


def eigengene(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    m = mat.copy()
    m.index = [str(i).upper() for i in m.index]
    m = m[~m.index.duplicated()]
    g = [x for x in genes if x in m.index]
    sub = m.loc[g]
    sub = sub[sub.std(axis=1) > 0]
    if len(sub) < 3:
        return pd.Series(dtype=float)
    return ((sub.T - sub.mean(axis=1)) / sub.std(axis=1)).mean(axis=1)


def main():
    mods = load_modules()
    # ネコのプロテオームを入れ忘れると、下の「モジュールレベルの種間相関」で
    # cat_prot_ctx_late の z が空になり、蛋白ペアが黙って skip される。
    datasets = ["cat_rna_ctx_early", "cat_rna_ctx_late", "cat_rna_ctx_prog",
                "cat_rna_med_early", "cat_rna_med_late", "cat_rna_med_prog",
                "cat_prot_ctx_early", "cat_prot_ctx_late", "cat_prot_ctx_prog",
                "cat_prot_med_early", "cat_prot_med_late", "cat_prot_med_prog",
                "m_rna_5d", "m_rna_2w", "m_rna_3w", "m_prot_d14", "m_prot_d21"]

    rows = []
    for ds in datasets:
        cat_side = ds.startswith("cat")
        v = lfc(ds, to_mouse=False if not cat_side else False)  # z は種内順位なので写像不要
        for name, genes in mods.items():
            z, mean_lfc, n = competitive_z(v, genes, cfg["thresholds"]["min_module_genes"])
            rows.append({"dataset": ds, "module": name, "z": z, "mean_lfc": mean_lfc, "n_genes": n})
    mz = pd.DataFrame(rows)
    mz.to_csv(RES / "module_competitive_z.csv", index=False)

    # モジュールレベルの種間相関
    xrows = []
    for cat_ds, mouse_ds in [("cat_rna_ctx_late", "m_rna_2w"),
                             ("cat_rna_med_late", "m_rna_2w"),
                             ("cat_prot_ctx_late", "m_prot_d21")]:
        a = mz[mz.dataset == cat_ds].set_index("module")["z"]
        b = mz[mz.dataset == mouse_ds].set_index("module")["z"]
        ok = a.notna() & b.notna()
        if ok.sum() < 5:
            log.warning("モジュール種間相関 %s x %s: 共通モジュールが %d 個しかないため skip",
                        cat_ds, mouse_ds, int(ok.sum()))
            continue
        rho, p = stats.spearmanr(a[ok], b[ok])
        xrows.append({"pair": f"{cat_ds} x {mouse_ds}", "rho": float(rho),
                      "p": float(p), "n_modules": int(ok.sum())})
    pd.DataFrame(xrows).to_csv(RES / "module_crossspecies.csv", index=False)

    # ---------- モジュール性の検定 ----------
    coh = []
    for tis, mat_name in [("ctx", "cat_rna_ctx"), ("med", "cat_rna_med")]:
        mat = pd.read_parquet(INT / f"{mat_name}.parquet")
        for name, genes in mods.items():
            obs, z, p, n = module_coherence(mat, genes, seed=SEED)
            coh.append({"tissue": tis, "module": name, "mean_r": obs,
                        "z_vs_null": z, "p": p, "n_genes": n,
                        "coherent": bool(np.isfinite(z) and z > 2)})
    cdf = pd.DataFrame(coh)
    cdf.to_csv(RES / "module_coherence.csv", index=False)
    weak = cdf[(~cdf.coherent) & cdf.z_vs_null.notna()]
    if len(weak):
        log.warning("共発現構造が弱いモジュール（解釈注意）:\n%s",
                    weak[["tissue", "module", "mean_r", "z_vs_null"]].to_string(index=False))

    # ---------- 検体単位の固有値（群間推移） ----------
    eig_rows = []
    for tis, mat_name, grp_name in [("ctx", "cat_rna_ctx", "cat_rna_ctx_grp"),
                                    ("med", "cat_rna_med", "cat_rna_med_grp")]:
        mat = pd.read_parquet(INT / f"{mat_name}.parquet")
        grp = pd.read_parquet(INT / f"{grp_name}.parquet").squeeze()
        for name, genes in mods.items():
            ev = eigengene(mat, genes)
            if ev.empty:
                continue
            rec = {"tissue": tis, "module": name}
            for g in ["Control", "CKD1/2", "CKD3/4"]:
                cols = grp[grp == g].index
                rec[g] = float(ev[cols].mean())
            # 進行（3/4 vs 1/2）の並べ替え検定
            a = ev[grp[grp == "CKD3/4"].index].to_numpy()
            b = ev[grp[grp == "CKD1/2"].index].to_numpy()
            obs = a.mean() - b.mean()
            pool = np.concatenate([a, b])
            rng = np.random.default_rng(SEED)
            null = np.array([(lambda p: p[:len(a)].mean() - p[len(a):].mean())(rng.permutation(pool))
                             for _ in range(cfg["permutation"]["n_perm_module"])])
            rec["prog_delta"] = float(obs)
            rec["prog_p_perm"] = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1))
            eig_rows.append(rec)
    edf = pd.DataFrame(eig_rows)
    edf.to_csv(RES / "module_eigengene_by_stage.csv", index=False)

    # ---------- MSigDB による独立検証 ----------
    for name in ["hallmark", "gobp", "kegg"]:
        gmt = load_gmt(REF / f"geneset_{name}.gmt")
        if not gmt:
            log.warning("%s の GMT がありません。00_fetch_refs.py を実行してください", name)
            continue
        rows = []
        for ds in ["cat_rna_ctx_late", "cat_rna_med_prog", "m_rna_2w", "m_prot_d21"]:
            v = lfc(ds)
            for term, genes in gmt.items():
                z, ml, n = competitive_z(v, genes, 10)
                if np.isfinite(z):
                    rows.append({"dataset": ds, "term": term, "z": z, "mean_lfc": ml, "n_genes": n})
        pd.DataFrame(rows).to_csv(RES / f"msigdb_{name}_competitive_z.csv", index=False)
        log.info("MSigDB %s: %d 行", name, len(rows))

    append_summary("04_modules", {
        "モジュール種間相関": xrows,
        "共発現が弱いモジュール": weak.module.unique().tolist() if len(weak) else "なし",
        "髄質で進行有意なモジュール(p<0.05)":
            edf.query("tissue=='med' and prog_p_perm<0.05")[["module", "prog_delta", "prog_p_perm"]]
               .round(4).to_dict("records"),
        "皮質で進行有意なモジュール(p<0.05)":
            edf.query("tissue=='ctx' and prog_p_perm<0.05")[["module", "prog_delta", "prog_p_perm"]]
               .round(4).to_dict("records"),
    }, cfg)


if __name__ == "__main__":
    main()
