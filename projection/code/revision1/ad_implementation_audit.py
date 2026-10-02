# -*- coding: utf-8 -*-
"""AD の実装を確定させる。Methods に書けるレベルまで、実測で裏を取る。

答えるべきこと:
  1. 各手法の実装（パッケージ・バージョン・関数・引数）
  2. 入力に何を渡したか（Δ 行列そのものか）
  3. ssGSEA と singscore のスコア行列が完全一致するのは代数的帰結か。
     正規化オプションを変えても一致するか、どの設定で壊れるか

出力: results/roundR/ad_implementation_audit.tsv
      results/roundR/ad_ssgsea_options.tsv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402
import pathway_scoring_methods as PSM  # noqa: E402

OUT = HERE / "results" / "roundR"
N_PROBE = 60          # オプションを振る検査は先頭 60 集合で足りる（一致は集合数に依らない）


def versions() -> dict:
    import gseapy, scipy, sklearn  # noqa: E401
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy.__version__,
            "gseapy": gseapy.__version__, "scikit-learn": sklearn.__version__}


def gseapy_matrix(fn, D, sets, **kw):
    import gseapy
    f = getattr(gseapy, fn)
    r = f(data=D, gene_sets=sets, outdir=None, min_size=PSM.MIN_GENES, max_size=100000,
          threads=4, seed=20260826, **kw)
    res = r.res2d if hasattr(r, "res2d") else r
    col = "ES" if "ES" in res else "NES"
    piv = res.pivot(index="Term", columns="Name", values=col).astype(float)
    nes = (res.pivot(index="Term", columns="Name", values="NES").astype(float)
           if "NES" in res else None)
    return piv, nes


def main() -> int:
    ver = versions()
    print("環境: " + " / ".join(f"{k} {v}" for k, v in ver.items()))

    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    groups = [gi for _c, _n, gi in paths]
    genes = list(D.index)
    print(f"入力: Δ 行列 {X.shape[0]} 遺伝子 x {X.shape[1]} 状態、経路 {len(paths)}")

    rows = []

    def rec(method, package, func, args, note, U, Cn):
        md = float(np.abs(np.asarray(U) - np.asarray(Cn)).max())
        rows.append({"method": method, "package": package, "version": ver.get(package, ""),
                     "function": func, "arguments": args,
                     "input_matrix": f"delta matrix {X.shape[0]}x{X.shape[1]}",
                     "score_column": note,
                     "max_abs_score_difference": md,
                     "invariant_to_per_state_shift": bool(md == 0.0)})
        print(f"  {method:22s} 中心化前後の最大差 {md:.3e}  → 不変 {md == 0.0}")

    print("\n1) 自前実装")
    rec("pathway mean (main)", "numpy", "score_mean (this repository)",
        "arithmetic mean over set genes, axis=0", "mean of delta",
        PSM.score_mean(X, groups), PSM.score_mean(C, groups))
    rec("singscore-style", "numpy", "score_singscore (this repository)",
        "within-state ranks of all genes (mergesort, ties broken by order), divided by n_genes, "
        "mean over set genes, minus 0.5; single undirected set, no up/down split, "
        "no absolute-deviation option", "normalized mean rank",
        PSM.score_singscore(X, groups), PSM.score_singscore(C, groups))
    rec("PLAGE-style", "numpy", "score_plage (this repository)",
        "genes standardized within the set across the 16 states (ddof=1); first right singular "
        "vector of the set matrix via numpy.linalg.svd; sign fixed so the coefficients sum "
        "positive", "first right singular vector",
        PSM.score_plage(X, groups), PSM.score_plage(C, groups))

    print("\n2) gseapy")
    sets = {f"{c}|{n}": [genes[i] for i in gi] for c, n, gi in paths}
    Xd = pd.DataFrame(X, index=D.index, columns=D.columns)
    Cd = pd.DataFrame(C, index=D.index, columns=D.columns)
    for label, fn, args in (
            ("ssGSEA", "ssgsea",
             "sample_norm_method='rank', correl_norm_type='rank', weight=0.25, "
             "ascending=False, permutation_num=None, min_size=30, max_size=100000, "
             "threads=4, seed=20260826 (defaults except size, threads, seed)"),
            ("GSVA", "gsva",
             "kcdf='Gaussian', weight=1.0, mx_diff=True, abs_rnk=False, min_size=30, "
             "max_size=100000, threads=4, seed=20260826 (defaults except size, threads, seed)")):
        u, _ = gseapy_matrix(fn, Xd, sets)
        c, _ = gseapy_matrix(fn, Cd, sets)
        common = [s for s in D.columns if s in u.columns and s in c.columns]
        keep = [t for t in u.index if t in c.index]
        rec(label, "gseapy", f"gseapy.{fn}", args, "ES (unnormalized)",
            u.loc[keep, common].to_numpy(float), c.loc[keep, common].to_numpy(float))

    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "ad_implementation_audit.tsv", sep="\t", index=False)

    # 3) ssGSEA の正規化オプションを振る
    print("\n3) ssGSEA のオプションを振ったときの不変性")
    sub = {f"{c}|{n}": [genes[i] for i in gi] for c, n, gi in paths[:N_PROBE]}
    opts = [("sample_norm_method='rank', correl_norm_type='rank' (used)", {}),
            ("correl_norm_type='zscore'", {"correl_norm_type": "zscore"}),
            ("correl_norm_type='symrank'", {"correl_norm_type": "symrank"}),
            ("sample_norm_method='log_rank'", {"sample_norm_method": "log_rank"}),
            ("sample_norm_method='log' (values, not ranks)", {"sample_norm_method": "log"}),
            ("sample_norm_method=None (input used as the metric)", {"sample_norm_method": None})]
    orows = []
    for label, kw in opts:
        try:
            u, un = gseapy_matrix("ssgsea", Xd, sub, **kw)
            c, cn = gseapy_matrix("ssgsea", Cd, sub, **kw)
            common = [s for s in D.columns if s in u.columns and s in c.columns]
            keep = [t for t in u.index if t in c.index]
            de = float(np.abs(u.loc[keep, common].to_numpy(float)
                              - c.loc[keep, common].to_numpy(float)).max())
            dn = (float(np.abs(un.loc[keep, common].to_numpy(float)
                               - cn.loc[keep, common].to_numpy(float)).max())
                  if un is not None else np.nan)
            orows.append({"option": label, "n_sets": len(keep),
                          "max_abs_ES_difference": de, "max_abs_NES_difference": dn,
                          "ES_invariant": bool(de == 0.0),
                          "NES_invariant": bool(dn == 0.0)})
            print(f"  {label:52s} ES 最大差 {de:.3e} / NES 最大差 {dn:.3e}")
        except Exception as e:                                  # noqa: BLE001
            orows.append({"option": label, "n_sets": 0, "max_abs_ES_difference": np.nan,
                          "max_abs_NES_difference": np.nan, "ES_invariant": False,
                          "NES_invariant": False})
            print(f"  {label:52s} 失敗 {type(e).__name__}: {e}")
    pd.DataFrame(orows).to_csv(OUT / "ad_ssgsea_options.tsv", sep="\t", index=False)

    # 一致が代数的であることの直接確認: 列ごとの順位が中心化で変わらないこと
    ru = np.argsort(np.argsort(X, axis=0, kind="mergesort"), axis=0, kind="mergesort")
    rc = np.argsort(np.argsort(C, axis=0, kind="mergesort"), axis=0, kind="mergesort")
    same = bool((ru == rc).all())
    print(f"\n  中心化前後で各状態内の遺伝子順位が完全に同一: {same}")
    pd.DataFrame([{"quantity": "within-state gene ranks identical after centering",
                   "value": float(same)},
                  {"quantity": "genes", "value": float(X.shape[0])},
                  {"quantity": "states", "value": float(X.shape[1])}]).to_csv(
        OUT / "ad_rank_invariance.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'ad_implementation_audit.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
