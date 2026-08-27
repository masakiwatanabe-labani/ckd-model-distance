"""IRI 9時点の遺伝子ごとの単調性で層別し、距離行列の時間相関がどう変わるかを見る。

定義:
  最早時点(IRI 2h)を基準とした log2FC の軌跡が、log10(経過日数) に対して
  単調かを Spearman で評価。|rho| > 0.7 を単調群、それ以外を非単調群。
  （n=9 時点で |rho|=0.683 が p=0.05 相当なので 0.7 はほぼその水準）

**循環性の警告**:
  単調群は「IRI の時系列で時間と相関する遺伝子」として選ばれている。
  その遺伝子集合で IRI 状態を含む距離行列を組めば、IRI-IRI ペアの
  時間相関が上がるのは定義上ほぼ自明。
  したがって、遺伝子選択に使っていない PodTRECK と UUO だけで組んだ
  5状態版（非循環）を必ず併記する。こちらが本当の検証。
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("18_monotonic_genes")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
REF = ROOT / cfg["paths"]["ref"]
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]
RHO_CUT = 0.7

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m16", Path(__file__).parent / "16_entry_model_checks.py")
m16 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m16)
m14, m12, m09 = m16.m14, m16.m12, m16.m09
S, DAYS, MODEL, ENTRY = m16.S, m16.DAYS, m16.MODEL, m16.ENTRY
STATES = m16.STATES

IRI_ORDER = ["IRI2h", "IRI4h", "IRI24h", "IRI48h", "IRI72h", "IRI7d", "IRI14d", "IRI28d", "IRI12m"]


def trajectories() -> pd.DataFrame:
    """最早時点(2h)を基準とした各時点の log2FC（遺伝子 x 時点）。"""
    mat = pd.read_parquet(INT / "mouse_iri_matrix.parquet")
    grp = pd.read_parquet(INT / "mouse_iri_grp.parquet").squeeze()
    means = {}
    for lab in IRI_ORDER:
        cols = grp[grp == lab].index
        means[lab] = np.log2(mat[cols] + 1).mean(axis=1)
    m = pd.DataFrame(means)
    keep = m.max(axis=1) > 1          # 低発現を除く
    m = m[keep]
    return m.sub(m["IRI2h"], axis=0)  # 2h 基準


def classify(traj: pd.DataFrame) -> pd.DataFrame:
    x = np.log10([DAYS[s] for s in IRI_ORDER])
    v = traj[IRI_ORDER].values
    rho = np.array([stats.spearmanr(x, row)[0] for row in v])
    amp = np.nanmax(v, axis=1) - np.nanmin(v, axis=1)
    out = pd.DataFrame({"gene": traj.index, "rho_time": rho, "amplitude": amp})
    out["group"] = np.where(np.abs(out.rho_time) > RHO_CUT,
                            np.where(out.rho_time > 0, "monotonic_up", "monotonic_down"),
                            "non_monotonic")
    out["is_monotonic"] = out.group != "non_monotonic"
    return out


def to_human(genes) -> set:
    s = pd.Series(1.0, index=list(genes))
    return set(m09.to_human_space(s, "mouse").index)


def mantel_for(states, gene_set):
    """指定遺伝子集合に限定して距離行列を組み、Mantel を返す。"""
    sub = {}
    for k in states:
        v = S[k]
        sub[k] = v[v.index.isin(gene_set)]
    n = len(states)
    D = np.zeros((n, n)); T = np.zeros((n, n)); E = np.zeros((n, n))
    npairs = []
    for i, j in combinations(range(n), 2):
        a, b = states[i], states[j]
        idx = sub[a].index.intersection(sub[b].index)
        npairs.append(len(idx))
        D[i, j] = D[j, i] = 1 - stats.spearmanr(sub[a].loc[idx], sub[b].loc[idx])[0]
        T[i, j] = T[j, i] = abs(np.log10(DAYS[a]) - np.log10(DAYS[b]))
        E[i, j] = E[j, i] = float(ENTRY[a] != ENTRY[b])
    rt, pt = m14.mantel(D, T, seed=SEED)
    re, pe = m14.mantel(D, E, seed=SEED)
    rte, pte = m14.mantel(D, T, seed=SEED, partial=E)
    return {"n_states": n, "median_overlap": int(np.median(npairs)),
            "time_rho": round(rt, 3), "time_p": round(pt, 4),
            "entry_rho": round(re, 3), "entry_p": round(pe, 4),
            "time|entry_rho": round(rte, 3), "time|entry_p": round(pte, 4)}


def enrich(target: set, background: set, top=12) -> pd.DataFrame:
    rows = []
    for lib in ["gobp", "kegg", "hallmark"]:
        p = REF / f"geneset_{lib}.gmt"
        if not p.exists():
            continue
        for line in p.read_text(encoding="utf-8").splitlines():
            f = line.split("\t")
            term, genes = f[0], {g.strip().upper() for g in f[2:] if g.strip()}
            gs = genes & background
            if not (5 <= len(gs) <= 500):
                continue
            k = len(gs & target)
            if k < 3:
                continue
            pv = stats.hypergeom.sf(k - 1, len(background), len(gs), len(target))
            rows.append({"library": lib, "term": term, "n_set": len(gs), "n_hit": k,
                         "expected": round(len(gs) * len(target) / len(background), 1),
                         "fold": round(k / (len(gs) * len(target) / len(background)), 2), "p": pv})
    if not rows:
        return pd.DataFrame()
    t = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    # Benjamini-Hochberg。ライブラリ（gobp/kegg/hallmark）ごとに独立して補正する。
    # 3ライブラリは互いに大きく重複するので、まとめて補正すると保守的すぎるため。
    t["q"] = np.nan
    for lib, idx in t.groupby("library").groups.items():
        pv = t.loc[idx, "p"].values
        order = np.argsort(pv)
        m = len(pv)
        q = np.empty(m)
        q[order] = np.minimum.accumulate(
            (pv[order] * m / np.arange(1, m + 1))[::-1])[::-1]
        t.loc[idx, "q"] = np.minimum(q, 1.0)
        t.loc[idx, "n_tested_in_library"] = m
    return t.sort_values("p").reset_index(drop=True)


def main():
    traj = trajectories()
    cls = classify(traj)
    cls.to_csv(RES / "iri_gene_monotonicity.csv", index=False)
    log.info("軌跡を評価した遺伝子 %d", len(cls))
    log.info("群構成: %s", cls.group.value_counts().to_dict())
    log.info("単調 %d (%.1f%%) / 非単調 %d",
             int(cls.is_monotonic.sum()), 100 * cls.is_monotonic.mean(),
             int((~cls.is_monotonic).sum()))

    mono_h = to_human(cls.loc[cls.is_monotonic, "gene"])
    non_h = to_human(cls.loc[~cls.is_monotonic, "gene"])
    all_h = to_human(cls.gene)
    log.info("ヒト空間: 単調 %d / 非単調 %d / 全体 %d", len(mono_h), len(non_h), len(all_h))

    res = {}
    for lab, gs in [("全遺伝子", all_h), ("単調群のみ", mono_h), ("非単調群のみ", non_h)]:
        res[f"14状態: {lab}"] = mantel_for(STATES, gs)
        log.info("[14状態 %s] %s", lab, res[f"14状態: {lab}"])
    noiri = [s for s in STATES if MODEL[s] != "IRI"]
    for lab, gs in [("全遺伝子", all_h), ("単調群のみ", mono_h), ("非単調群のみ", non_h)]:
        res[f"5状態(非循環, PodTRECK+UUO): {lab}"] = mantel_for(noiri, gs)
        log.info("[5状態(非循環) %s] %s", lab, res[f"5状態(非循環, PodTRECK+UUO): {lab}"])

    pd.DataFrame(res).T.to_csv(RES / "monotonic_mantel_comparison.csv")

    en = enrich(mono_h, all_h)
    if len(en):
        en.to_csv(RES / "monotonic_enrichment.csv", index=False)
        log.info("\n単調群の濃縮 上位 (BH補正済み):\n%s",
                 en.head(15)[["library", "term", "n_hit", "n_set", "fold", "p", "q"]].to_string(index=False))
        log.info("q<0.05 の項目数: %s",
                 en[en.q < 0.05].library.value_counts().to_dict())

    up = to_human(cls.loc[cls.group == "monotonic_up", "gene"])
    dn = to_human(cls.loc[cls.group == "monotonic_down", "gene"])
    for lab, gs in [("up", up), ("dn", dn)]:
        e = enrich(gs, all_h)
        if len(e):
            e.to_csv(RES / f"monotonic_{lab}_enrichment.csv", index=False)
            log.info("\n単調%s 上位 (BH補正済み):\n%s", lab,
                     e.head(8)[["library", "term", "n_hit", "fold", "p", "q"]].to_string(index=False))

    append_summary("18_monotonic_genes / 単調性による層別", {
        "定義": f"IRI 2h基準の log2FC 軌跡 vs log10(日数) の |Spearman rho| > {RHO_CUT}",
        "群構成": cls.group.value_counts().to_dict(),
        "ヒト空間の遺伝子数": {"単調": len(mono_h), "非単調": len(non_h), "全体": len(all_h)},
        "Mantel比較": res,
        "循環性の警告": ("単調群は IRI の時系列で時間相関する遺伝子として選ばれているため、"
                         "IRI 状態を含む14状態版で時間相関が上がるのは定義上ほぼ自明。"
                         "遺伝子選択に使っていない PodTRECK+UUO の5状態版（非循環）が本当の検証。"),
        "濃縮": ("results/monotonic_enrichment.csv, monotonic_up_enrichment.csv, "
                 "monotonic_dn_enrichment.csv"),
        "多重比較補正": ("Benjamini-Hochberg。gobp/kegg/hallmark の各ライブラリ内で独立に補正"
                         "（3ライブラリは項目が大きく重複するため、まとめて補正すると保守的すぎる）。"
                         "列 q が BH-adjusted p、n_tested_in_library が各ライブラリの検定数。"),
    }, cfg)


if __name__ == "__main__":
    main()
