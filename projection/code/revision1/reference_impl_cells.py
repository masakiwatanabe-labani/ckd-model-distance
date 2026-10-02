# -*- coding: utf-8 -*-
"""AI2c: 参照実装（Bioconductor）と、これまで使ってきた実装をセル単位で突き合わせる。

比較の相手:
  ssGSEA    gseapy.ssgsea (ES)        ←→ ssgseaParam(normalize = FALSE)
  GSVA      gseapy.gsva               ←→ gsvaParam(kcdf = "Gaussian")
  PLAGE     自前実装                   ←→ plageParam
  singscore 自前実装                   ←→ singscore::simpleScore

行列は 422 x 16。PLAGE は集合ごとに符号が任意なので、符号を揃えてから比べる。
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

R = HERE / "results" / "roundR"
IN = R / "r_input"


def main() -> int:
    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    names = [f"{c}__{n}" for c, n, _g in paths]
    groups = [gi for _c, _n, gi in paths]
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    sets = {nm: [D.index[i] for i in gi] for nm, gi in zip(names, groups)}

    def gseapy_mat(fn, M):
        import gseapy
        f = getattr(gseapy, fn)
        r = f(data=pd.DataFrame(M, index=D.index, columns=D.columns), gene_sets=sets,
              outdir=None, min_size=PSM.MIN_GENES, max_size=100000, threads=4, seed=20260826)
        res = r.res2d if hasattr(r, "res2d") else r
        return res.pivot(index="Term", columns="Name", values="ES").astype(float)

    own = {}
    for tag, M in (("uncentered", X), ("centered", C)):
        own[("singscore", tag)] = pd.DataFrame(PSM.score_singscore(M, groups),
                                               index=names, columns=D.columns)
        own[("plage", tag)] = pd.DataFrame(PSM.score_plage(M, groups),
                                           index=names, columns=D.columns)
        own[("ssgsea", tag)] = gseapy_mat("ssgsea", M)
        own[("gsva", tag)] = gseapy_mat("gsva", M)

    LAB = {"ssgsea": ("ssGSEA", "gseapy.ssgsea (ES)", "ssgseaParam(normalize = FALSE)"),
           "gsva": ("GSVA", "gseapy.gsva", 'gsvaParam(kcdf = "Gaussian")'),
           "plage": ("PLAGE", "own implementation", "plageParam"),
           "singscore": ("singscore", "own implementation", "singscore::simpleScore")}
    rows = []
    for key, (label, a_lab, b_lab) in LAB.items():
        for tag in ("uncentered", "centered"):
            a = own[(key, tag)]
            b = pd.read_csv(R / f"r_{key}_{tag}.tsv", sep="\t", index_col=0)
            s = [x for x in a.index if x in b.index]
            c = [x for x in D.columns if x in a.columns and x in b.columns]
            A, B = a.loc[s, c].to_numpy(float), b.loc[s, c].to_numpy(float)
            flipped = 0
            if key == "plage":
                for i in range(A.shape[0]):
                    if np.dot(A[i], B[i]) < 0:
                        B[i] = -B[i]
                        flipped += 1
            d = np.abs(A - B)
            scale = max(np.abs(A).max(), np.abs(B).max())
            rows.append({"method": label, "compared": f"{a_lab} vs {b_lab}", "matrix": tag,
                         "n_sets": len(s), "n_states": len(c),
                         "max_abs_difference": float(d.max()),
                         "median_abs_difference": float(np.median(d)),
                         "max_relative_difference": float(d.max() / scale) if scale else 0.0,
                         "pearson_r_of_cells": float(np.corrcoef(A.ravel(), B.ravel())[0, 1]),
                         "spearman_r_of_cells": float(pd.Series(A.ravel()).corr(
                             pd.Series(B.ravel()), method="spearman")),
                         "plage_sign_flipped_sets": flipped})
            print(f"  {label:10s} {tag:11s} 最大差 {d.max():.3e} / 中央 {np.median(d):.3e} / "
                  f"相対 {rows[-1]['max_relative_difference']:.3e} / "
                  f"r {rows[-1]['pearson_r_of_cells']:.6f}")

    T = pd.DataFrame(rows)
    T.to_csv(R / "reference_implementation_cells.tsv", sep="\t", index=False)
    print(f"\n書き出し: {R / 'reference_implementation_cells.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
