"""モジュールレベルでの距離行列と Mantel 検定。

Fig2C で「モジュールに集約すると種間一致が種内ベンチマークを上回る」ことが
見えたので、その一致がモデル間でも成立するか、そしてモジュールレベルでも
「入口ではなく時間」という構図が保たれるかを確認する。

手順:
  1. 全16状態（マウス14 + ネコ2）について、各モジュールの競合的 z を算出。
     z は各データセット内での順位パーセンタイルの標準化なので、
     状態間で直接比較できる（log2FC の絶対値には依存しない）。
  2. 状態間の Spearman rho（z ベクトル間）→ 距離 1 - rho。
  3. ペア分類（モデル間 / モデル内 / 種間）ごとの分布を遺伝子レベルと並置。
  4. 入口・モデル・時間の Mantel / 偏 Mantel。
  5. モジュール数が少ないので、検出力とモジュール定義依存性を確認する。
     MSigDB Hallmark（50セット）で組み直した版と比較する。
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, competitive_z, append_summary  # noqa: E402

log = get_logger("24_module_level")
cfg = load_config()
ROOT = Path(cfg["_root"])
REF = ROOT / cfg["paths"]["ref"]
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]
MIN_G = cfg["thresholds"]["min_module_genes"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m19", Path(__file__).parent / "19_distance_distributions.py")
_spec_only = _s  # 実行はしない（main が走ってしまうため）
_s2 = importlib.util.spec_from_file_location("m16", Path(__file__).parent / "16_entry_model_checks.py")
m16 = importlib.util.module_from_spec(_s2); _s2.loader.exec_module(m16)
m14 = m16.m14
S, DAYS, MODEL, ENTRY, STATES = m16.S, m16.DAYS, m16.MODEL, m16.ENTRY, m16.STATES
CAT = m16.CAT

POOL = {**S, **CAT}
ALLST = STATES + list(CAT)
MODEL_C = {**MODEL, "cat_natural_ctx": "cat", "cat_natural_med": "cat"}
ENTRY_C = {**ENTRY, "cat_natural_ctx": "tubular", "cat_natural_med": "tubular"}


def curated_modules() -> dict:
    with open(ROOT / "config" / "modules.yaml", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return {k: str(v).split() for k, v in raw.items()}


def hallmark_modules(min_n=10, max_n=500) -> dict:
    p = REF / "geneset_hallmark.gmt"
    out = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        genes = [g.strip().upper() for g in f[2:] if g.strip()]
        if min_n <= len(genes) <= max_n:
            out[f[0]] = genes
    return out


def zmatrix(mods: dict) -> pd.DataFrame:
    """状態 x モジュール の z 行列。"""
    rows = {}
    for st in ALLST:
        v = POOL[st]
        rows[st] = {name: competitive_z(v, genes, MIN_G)[0] for name, genes in mods.items()}
    Z = pd.DataFrame(rows).T
    return Z


def dist_matrix(Z: pd.DataFrame):
    st = list(Z.index)
    n = len(st)
    D = np.zeros((n, n))
    npair = []
    for i, j in combinations(range(n), 2):
        a, b = Z.iloc[i], Z.iloc[j]
        ok = a.notna() & b.notna()
        npair.append(int(ok.sum()))
        r = stats.spearmanr(a[ok], b[ok])[0] if ok.sum() >= 5 else np.nan
        D[i, j] = D[j, i] = 1 - r
    return D, st, int(np.median(npair))


def classify(st):
    n = len(st)
    kinds = {}
    for i, j in combinations(range(n), 2):
        a, b = st[i], st[j]
        ca, cb = MODEL_C[a] == "cat", MODEL_C[b] == "cat"
        if ca and cb:
            k = "within_cat"
        elif ca or cb:
            k = "cross_species"
        elif MODEL_C[a] == MODEL_C[b]:
            k = "within_model"
        else:
            k = "between_model"
        kinds[(i, j)] = k
    return kinds


def desc(x, name):
    x = np.asarray(x); x = x[np.isfinite(x)]
    return {"level": name, "n": len(x), "min": round(float(x.min()), 4),
            "Q1": round(float(np.percentile(x, 25)), 4),
            "median": round(float(np.median(x)), 4),
            "Q3": round(float(np.percentile(x, 75)), 4),
            "max": round(float(x.max()), 4)}


def analyse(mods: dict, tag: str):
    Z = zmatrix(mods)
    Z.to_csv(RES / f"module_z_matrix_{tag}.csv")
    usable = int(Z.notna().all(axis=0).sum())
    log.info("[%s] モジュール %d 個、全状態で z が出たもの %d 個", tag, Z.shape[1], usable)
    D, st, med_overlap = dist_matrix(Z)
    kinds = classify(st)
    iu = np.triu_indices(len(st), 1)

    rows = []
    for k, (i, j) in enumerate(zip(*iu)):
        rows.append({"a": st[i], "b": st[j], "kind": kinds[(i, j)], "distance": D[i, j]})
    t = pd.DataFrame(rows)
    t.to_csv(RES / f"module_distance_pairs_{tag}.csv", index=False)

    summ = [desc(t[t.kind == k].distance, k)
            for k in ["between_model", "within_model", "cross_species"]]
    summ.append({"level": "within_cat (reference)", "n": 1,
                 "median": round(float(t[t.kind == "within_cat"].distance.iloc[0]), 4)})

    # Mantel
    n = len(st)
    T = np.zeros((n, n)); E = np.zeros((n, n)); M = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        a, b = st[i], st[j]
        E[i, j] = E[j, i] = float(ENTRY_C[a] != ENTRY_C[b])
        M[i, j] = M[j, i] = float(MODEL_C[a] != MODEL_C[b])
    # 時間はマウスのみ（ネコに日数が無い）
    mouse_idx = [i for i, s in enumerate(st) if MODEL_C[s] != "cat"]
    Dm = D[np.ix_(mouse_idx, mouse_idx)]
    stm = [st[i] for i in mouse_idx]
    nm = len(stm)
    Tm = np.zeros((nm, nm)); Em = np.zeros((nm, nm)); Mm = np.zeros((nm, nm))
    for i, j in combinations(range(nm), 2):
        a, b = stm[i], stm[j]
        Tm[i, j] = Tm[j, i] = abs(np.log10(DAYS[a]) - np.log10(DAYS[b]))
        Em[i, j] = Em[j, i] = float(ENTRY_C[a] != ENTRY_C[b])
        Mm[i, j] = Mm[j, i] = float(MODEL_C[a] != MODEL_C[b])

    man = {}
    for nm_, X in [("onset compartment (16 states)", E), ("model identity (16 states)", M)]:
        r, p = m14.mantel(D, X, seed=SEED)
        man[nm_] = {"rho": round(r, 3), "p": round(p, 4)}
    for nm_, X, part in [("elapsed time (14 mouse states)", Tm, None),
                         ("onset compartment (14 mouse)", Em, None),
                         ("model identity (14 mouse)", Mm, None),
                         ("time | onset compartment", Tm, Em),
                         ("onset compartment | time", Em, Tm)]:
        r, p = m14.mantel(Dm, X, seed=SEED, partial=part)
        man[nm_] = {"rho": round(r, 3), "p": round(p, 4)}
    for k, v in man.items():
        log.info("  [%s] Mantel %-32s rho=%+.3f p=%.4f", tag, k, v["rho"], v["p"])
    return {"n_modules": Z.shape[1], "usable_all_states": usable,
            "median_modules_per_pair": med_overlap,
            "distribution": summ, "mantel": man}, t


def main():
    cur, t_cur = analyse(curated_modules(), "curated")
    hal, t_hal = analyse(hallmark_modules(), "hallmark")

    # 遺伝子レベル（19番の出力）と並置
    gl = pd.read_csv(RES / "distance_distributions_16state.csv")
    gene_summ = [desc(gl[gl.kind == k].distance, k)
                 for k in ["between_model", "within_model", "cross_species"]]

    comp = []
    for lab, rows in [("gene level", gene_summ),
                      ("module level (curated, 27 defined)", cur["distribution"]),
                      ("module level (MSigDB Hallmark, 50)", hal["distribution"])]:
        for r in rows:
            if r["level"].startswith("within_cat"):
                continue
            comp.append({"level": lab, **r})
    cdf = pd.DataFrame(comp)
    cdf.to_csv(RES / "module_vs_gene_distribution.csv", index=False)
    log.info("\n%s", cdf.to_string(index=False))

    # 2つのモジュール定義間で距離が一致するか（定義依存性）
    m1 = t_cur.set_index(["a", "b"]).distance
    m2 = t_hal.set_index(["a", "b"]).distance
    idx = m1.index.intersection(m2.index)
    rho_def = stats.spearmanr(m1.loc[idx], m2.loc[idx])[0]
    gl2 = gl.set_index(["a", "b"]).distance
    idx2 = m1.index.intersection(gl2.index)
    rho_gene = stats.spearmanr(m1.loc[idx2], gl2.loc[idx2])[0]
    log.info("定義依存性: curated vs Hallmark の距離 rho=%.3f (n=%d)", rho_def, len(idx))
    log.info("モジュール(curated) vs 遺伝子レベル の距離 rho=%.3f (n=%d)", rho_gene, len(idx2))

    append_summary("24_module_level / モジュールレベルの距離と Mantel", {
        "curated (config/modules.yaml)": cur,
        "MSigDB Hallmark": hal,
        "分布の並置": cdf.to_dict("records"),
        "定義依存性": {
            "curated vs Hallmark の距離相関": round(float(rho_def), 3),
            "curated モジュール vs 遺伝子レベル距離の相関": round(float(rho_gene), 3)},
        "検出力の注記": ("モジュールは 27 定義中、全状態で z が算出できたものだけが使われる。"
                         "1ペアあたりの実効次元は中央値で curated %d / Hallmark %d。"
                         "Spearman rho の推定は次元数が小さいほど不安定で、"
                         "Mantel の p も遺伝子レベルより弱くなる。"
                         % (cur["median_modules_per_pair"], hal["median_modules_per_pair"])),
    }, cfg)


if __name__ == "__main__":
    main()
