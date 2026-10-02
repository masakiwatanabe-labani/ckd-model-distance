# -*- coding: utf-8 -*-
"""AI: Table S4 の行を一箇所で組み立てる。

報告値は Bioconductor の参照実装（reference_impl_all.R → reference_impl_auc.py）から取る。
経路平均だけは参照実装に対応物がないので Python 側の表から取る。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
R = HERE / "results" / "roundR"

REF_ROWS = [("ssGSEA", "ssGSEA (ES, normalize = FALSE)",
             "GSVA 2.6.6, ssgseaParam + gsva"),
            ("GSVA", "GSVA (kcdf = Gaussian)", "GSVA 2.6.6, gsvaParam + gsva"),
            ("singscore", "singscore", "singscore 1.32.0, simpleScore"),
            ("PLAGE", "PLAGE", "GSVA 2.6.6, plageParam + gsva")]
COLS = ["method", "n_sets", "auc_uncentered", "auc_centered", "difference",
        "within_median_uncentered", "cross_median_uncentered",
        "max_abs_score_difference", "invariant_to_per_state_shift", "implementation"]


def main() -> int:
    A = pd.read_csv(R / "reference_implementation_auc.tsv", sep="\t")
    A = A[A.sign_convention == "as returned"].set_index("method")
    PY = pd.read_csv(R / "pathway_scoring_methods.tsv", sep="\t").set_index("method")

    rows = []
    m = PY.loc["pathway mean (main)"]
    rows.append({"method": "pathway mean (main)", "n_sets": int(m.n_sets),
                 "auc_uncentered": m.auc_uncentered, "auc_centered": m.auc_centered,
                 "difference": m.difference,
                 "within_median_uncentered": m.within_median_uncentered,
                 "cross_median_uncentered": m.cross_median_uncentered,
                 "max_abs_score_difference": m.max_abs_score_difference,
                 "invariant_to_per_state_shift": bool(m.invariant_to_per_state_shift),
                 "implementation": "arithmetic mean, NumPy"})
    for label, key, impl in REF_ROWS:
        r = A.loc[key]
        rows.append({"method": label, "n_sets": int(r.n_sets),
                     "auc_uncentered": r.auc_uncentered, "auc_centered": r.auc_centered,
                     "difference": r.difference,
                     "within_median_uncentered": r.within_median_uncentered,
                     "cross_median_uncentered": r.cross_median_uncentered,
                     "max_abs_score_difference": r.max_abs_score_difference,
                     "invariant_to_per_state_shift": bool(r.invariant_to_per_state_shift),
                     "implementation": impl})
    T = pd.DataFrame(rows)[COLS]
    T.to_csv(R / "table_s4_rows.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(T.round(4).to_string(index=False))
    print(f"\n書き出し: {R / 'table_s4_rows.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
