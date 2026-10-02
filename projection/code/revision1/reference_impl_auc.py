# -*- coding: utf-8 -*-
"""AI2b: Bioconductor 参照実装のスコア行列から、87 ペアの集約 AUC を計算する。

AUC は主解析と同じ関数で計算する。R には計算させない。
PLAGE の第1特異ベクトルは符号が任意なので、符号の扱いを 2 通り出す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402

R = HERE / "results" / "roundR"
METHODS = [("ssGSEA (ES, normalize = FALSE)", "ssgsea"),
           ("ssGSEA (normalize = TRUE)", "ssgseanorm"),
           ("GSVA (kcdf = Gaussian)", "gsva"),
           ("GSVA (kcdf = auto)", "gsvaauto"),
           ("PLAGE", "plage"),
           ("singscore", "singscore")]


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def main() -> int:
    D, P = PCS.load()
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = ((P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))).to_numpy()

    def agg(A):
        An = A / np.linalg.norm(A, axis=0, keepdims=True)
        c = (An.T @ An)[ia, ib]
        return auc(c[m & wi], c[m & ~wi]), float(np.median(c[m & wi])), \
            float(np.median(c[m & ~wi]))

    rows = []
    for label, key in METHODS:
        mats = {}
        for tag in ("uncentered", "centered"):
            f = R / f"r_{key}_{tag}.tsv"
            if not f.exists():
                print(f"無い: {f.name}", file=sys.stderr)
                return 2
            mats[tag] = pd.read_csv(f, sep="\t", index_col=0)
        sets = [s for s in mats["uncentered"].index if s in mats["centered"].index]
        cols = [c for c in D.columns if c in mats["uncentered"].columns]
        U = mats["uncentered"].loc[sets, cols].to_numpy(float)
        C = mats["centered"].loc[sets, cols].to_numpy(float)
        md = float(np.abs(U - C).max())
        variants = [("as returned", U, C)]
        if key == "plage":
            # 符号を集合ごとに「係数の和が正」へ揃えた版も出す
            def fix(A):
                A = A.copy()
                s = A.sum(axis=1)
                A[s < 0] *= -1
                return A
            variants.append(("sign fixed so coefficients sum positive", fix(U), fix(C)))
        for vlab, Uv, Cv in variants:
            au, wu, cu = agg(Uv)
            ac, wc, cc = agg(Cv)
            rows.append({"method": label, "sign_convention": vlab, "n_sets": len(sets),
                         "auc_uncentered": au, "auc_centered": ac, "difference": au - ac,
                         "max_abs_score_difference": float(np.abs(Uv - Cv).max()),
                         "invariant_to_per_state_shift": float(np.abs(Uv - Cv).max()) == 0.0,
                         "within_median_uncentered": wu, "cross_median_uncentered": cu})
            print(f"  {label:34s} {vlab[:22]:22s} 非中心化 {au:.4f} / 中心化 {ac:.4f} / "
                  f"差 {au - ac:+.4f} / 最大差 {md:.3e}")

    T = pd.DataFrame(rows)
    T.to_csv(R / "reference_implementation_auc.tsv", sep="\t", index=False)
    print(f"\n書き出し: {R / 'reference_implementation_auc.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
