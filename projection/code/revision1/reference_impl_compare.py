# -*- coding: utf-8 -*-
"""AI2: R の参照実装と自前実装を突き合わせる。

R が書いたスコア行列を読み、Python 側の行列とセル単位で比べ、
同じ AUC 関数で 87 ペアの集約 AUC を両方について計算する。
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
IN = R / "r_input"


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def main() -> int:
    need = [R / f"r_{m}_{t}.tsv" for m in ("plage", "singscore")
            for t in ("uncentered", "centered")]
    missing = [p.name for p in need if not p.exists()]
    if missing:
        print(f"R の出力が無い: {missing}\n"
              f"  先に Rscript code/revision1/reference_impl_compare.R {HERE} を実行する",
              file=sys.stderr)
        return 2

    D, P = PCS.load()
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    P = P.copy()
    P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
    wi = P.wi.to_numpy()

    def agg_auc(A):
        An = A / np.linalg.norm(A, axis=0, keepdims=True)
        c = (An.T @ An)[ia, ib]
        return auc(c[m & wi], c[m & ~wi])

    rows = []
    for method in ("singscore", "plage"):
        mats = {}
        for tag in ("uncentered", "centered"):
            py = pd.read_csv(IN / f"py_{method}_{tag}.tsv", sep="\t", index_col=0)
            rr = pd.read_csv(R / f"r_{method}_{tag}.tsv", sep="\t", index_col=0)
            common_sets = [s for s in py.index if s in rr.index]
            cols = [c for c in py.columns if c in rr.columns]
            a = py.loc[common_sets, cols].to_numpy(float)
            b = rr.loc[common_sets, cols].to_numpy(float)
            # PLAGE の第1特異ベクトルは符号が任意。集合ごとに符号を揃えてから比べる
            flipped = 0
            if method == "plage":
                for i in range(a.shape[0]):
                    if np.dot(a[i], b[i]) < 0:
                        b[i] = -b[i]
                        flipped += 1
            mats[tag] = (a, b, common_sets, cols, flipped)
            d = np.abs(a - b)
            corr = float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
            rows.append({"method": method, "matrix": tag, "n_sets": len(common_sets),
                         "n_states": len(cols), "sets_missing_in_r": len(py) - len(common_sets),
                         "max_abs_difference": float(d.max()),
                         "median_abs_difference": float(np.median(d)),
                         "pearson_r_of_cells": corr,
                         "plage_sign_flipped_sets": flipped,
                         "auc_python": agg_auc(a), "auc_r": agg_auc(b)})
        au, ac = rows[-2], rows[-1]
        rows.append({"method": method, "matrix": "centering difference", "n_sets": au["n_sets"],
                     "n_states": au["n_states"], "sets_missing_in_r": au["sets_missing_in_r"],
                     "max_abs_difference": np.nan, "median_abs_difference": np.nan,
                     "pearson_r_of_cells": np.nan, "plage_sign_flipped_sets": np.nan,
                     "auc_python": au["auc_python"] - ac["auc_python"],
                     "auc_r": au["auc_r"] - ac["auc_r"]})

    T = pd.DataFrame(rows)
    T.to_csv(R / "reference_implementation_comparison.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(T.to_string(index=False))
    print(f"\n書き出し: {R / 'reference_implementation_comparison.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
