# -*- coding: utf-8 -*-
"""AD: 経路スコアリング手法を変えて中心化の対比を見る（査読者 2 comment 3）。

単純平均に加えて ssGSEA と GSVA を回し、非中心化／中心化の AUC を並べる。
どちらも遺伝子をサンプル内で順位化するので、状態ごとの全体水準はスコアに残らない。
事前の予測は ~/Desktop/CKD_xspecies_review/AD_prediction.md に実行前に記録した。

singscore と PLAGE も、実装が軽いので同じ枠で計算する。
  singscore: 集合内遺伝子の順位平均を [-0.5, 0.5] に正規化
  PLAGE    : 集合の遺伝子×状態の行列を遺伝子ごとに標準化し、第1主成分の係数
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

OUT = HERE / "results" / "roundR"
MIN_GENES = 30
# gseapy 側の手法が Δ 行列をどう読むか。不変性そのものは実測から導く
MECHANISM = {
    "ssGSEA": "within-state gene ranks only (gseapy, sample_norm_method='rank')",
    "GSVA": "kernel CDF of the values across states, then within-state ranks (gseapy)",
}
# この表は Python 側の相互確認用。報告値は Bioconductor 参照実装（reference_impl_all.R）から取る
LABEL = {"ssGSEA": "ssGSEA (gseapy)", "GSVA": "GSVA (gseapy)"}


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def score_mean(X, groups):
    return np.vstack([X[g].mean(axis=0) for g in groups])


def score_singscore(X, groups):
    """状態ごとに全遺伝子を順位化し、集合内の平均順位を中心化して返す。"""
    n = X.shape[0]
    R = np.empty_like(X)
    for j in range(X.shape[1]):
        R[:, j] = (np.argsort(np.argsort(X[:, j], kind="mergesort"),
                              kind="mergesort") + 1) / n
    return np.vstack([R[g].mean(axis=0) - 0.5 for g in groups])


def score_plage(X, groups):
    """集合ごとに遺伝子を標準化し、第1主成分の状態側の係数を返す。"""
    out = []
    for g in groups:
        A = X[g]
        A = A - A.mean(axis=1, keepdims=True)
        sd = A.std(axis=1, ddof=1, keepdims=True)
        A = A / np.where(sd > 0, sd, 1.0)
        u, s, vt = np.linalg.svd(A, full_matrices=False)
        v = vt[0]
        if np.sum(v) < 0:                 # 符号の任意性を固定する
            v = -v
        out.append(v)
    return np.vstack(out)


def run_gseapy(fn, D, sets, name):
    import gseapy
    df = fn(data=D, gene_sets=sets, outdir=None, min_size=MIN_GENES, max_size=100000,
            threads=4, seed=20260826)
    res = df.res2d if hasattr(df, "res2d") else df
    piv = res.pivot(index="Term", columns="Name", values="ES" if "ES" in res else "NES")
    piv = piv.astype(float)
    print(f"  {name}: {piv.shape[0]} セット x {piv.shape[1]} 状態")
    return piv


def main() -> int:
    missing = [fn for fn in PCS.COLL_FULL.values() if not (PCS.REF / fn).exists()]
    if missing:
        print(f"GMT が無い: {missing}", file=sys.stderr)
        return 2

    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()
    groups = [gi for _c, _n, gi in paths]

    def cos_auc(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        c = (Mn.T @ Mn)[ia, ib]
        return auc(c[m & wi], c[m & ~wi]), float(np.median(c[m & wi])), \
            float(np.median(c[m & ~wi])), float(c[m].max() - c[m].min())

    # 遺伝子セットは gseapy 用に名前→遺伝子名の辞書にする
    genes = list(D.index)
    sets = {f"{c}|{n}": [genes[i] for i in gi] for c, n, gi in paths}

    rows = []

    def add(method, uncen, cen, mechanism=""):
        """note は手で書かない。実測した最大差から不変性を導く（AH1）。"""
        a_u, wu, cu, su = cos_auc(uncen)
        a_c, wc, cc, sc = cos_auc(cen)
        # 中心化前後でスコア行列そのものが動くか（0 なら手法が暗黙に中心化している）
        md = float(np.abs(np.asarray(uncen) - np.asarray(cen)).max())
        inv = md == 0.0
        note = (f"{mechanism}; " if mechanism else "") + (
            "invariant to per-state centering: score matrix identical (max abs difference 0)"
            if inv else
            f"not invariant to per-state centering (max abs difference {md:.4f})")
        rows.append({"method": method, "n_sets": uncen.shape[0],
                     "auc_uncentered": a_u, "auc_centered": a_c,
                     "difference": a_u - a_c,
                     "max_abs_score_difference": md,
                     "invariant_to_per_state_shift": inv,
                     "within_median_uncentered": wu, "cross_median_uncentered": cu,
                     "spread_uncentered": su,
                     "within_median_centered": wc, "cross_median_centered": cc,
                     "note": note})
        print(f"  {method:26s} 非中心化 {a_u:.4f} / 中心化 {a_c:.4f} / 差 {a_u - a_c:+.4f}"
              f"  / 不変 {inv}")

    print("経路スコアリング手法ごとの中心化の対比:")
    add("pathway mean (main)", score_mean(X, groups), score_mean(C, groups),
        "arithmetic mean of the values")
    add("singscore (own implementation)", score_singscore(X, groups), score_singscore(C, groups),
        "within-state gene ranks only; own implementation, not the Bioconductor package")
    add("PLAGE (own implementation)", score_plage(X, groups), score_plage(C, groups),
        "first right singular vector of the values standardized across states; own implementation, not the Bioconductor package")

    for label, fn in (("ssGSEA", "ssgsea"), ("GSVA", "gsva")):
        import gseapy
        f = getattr(gseapy, fn)
        try:
            u = run_gseapy(f, pd.DataFrame(X, index=D.index, columns=D.columns), sets, label)
            c = run_gseapy(f, pd.DataFrame(C, index=D.index, columns=D.columns), sets,
                           label + " (centered)")
            common = [s for s in D.columns if s in u.columns and s in c.columns]
            keep = [t for t in u.index if t in c.index]
            add(LABEL[label], u.loc[keep, common].to_numpy(float),
                c.loc[keep, common].to_numpy(float), MECHANISM[label])
        except Exception as e:                                   # noqa: BLE001
            print(f"  {label}: 失敗 {type(e).__name__}: {e}")
            rows.append({"method": LABEL.get(label, label), "n_sets": 0, "auc_uncentered": np.nan,
                         "auc_centered": np.nan, "difference": np.nan,
                         "invariant_to_per_state_shift": None,
                         "note": f"failed: {type(e).__name__}"})

    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "pathway_scoring_methods.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'pathway_scoring_methods.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
