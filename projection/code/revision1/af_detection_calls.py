# -*- coding: utf-8 -*-
"""AF3: 検出コールそのものを補足として出せる形にする。

Δ 行列の全遺伝子について、3 つのプロテオームでの検出可否と Group A 所属を書き出す。
査読者が再計算せずに選択を確認できるようにするのが目的。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import build_delta_matrix as BDM  # noqa: E402

OUT = HERE / "results" / "roundR"
FRAC = 0.5


def detected(name: str, species: str) -> set:
    d = BDM.read_int(name)
    s = pd.Series(1.0, index=d.index[d.notna().mean(axis=1) >= FRAC])
    h, _ = BDM.to_human(s, species)
    return set(h.index)


def main() -> int:
    ctx = detected("cat_prot_ctx", "cat")
    med = detected("cat_prot_med", "cat")
    mo = detected("mouse_prot", "mouse")
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    ga = set((HERE / "groupA_intersection.txt").read_text().split())
    gu = set((HERE / "groupA_union.txt").read_text().split())
    complete = set(D.dropna(axis=0, how="any").index)

    T = pd.DataFrame({
        "gene": sorted(D.index),
    })
    T["feline_cortex_proteome"] = ["yes" if g in ctx else "no" for g in T.gene]
    T["feline_medulla_proteome"] = ["yes" if g in med else "no" for g in T.gene]
    T["pod_treck_proteome"] = ["yes" if g in mo else "no" for g in T.gene]
    T["group"] = ["A" if g in ga else "B" for g in T.gene]
    T["group_a_union_variant"] = ["yes" if g in gu else "no" for g in T.gene]
    T["complete_in_16_states"] = ["yes" if g in complete else "no" for g in T.gene]
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "groupA_detection_calls.tsv", sep="\t", index=False)

    n = lambda c, v="yes": int((T[c] == v).sum())   # noqa: E731
    S = pd.DataFrame([
        {"step": "Genes in the Δ matrix", "genes": len(T)},
        {"step": "Detected in the feline cortical proteome", "genes": n("feline_cortex_proteome")},
        {"step": "Detected in the feline medullary proteome",
         "genes": n("feline_medulla_proteome")},
        {"step": "Detected in both feline proteomes",
         "genes": int(((T.feline_cortex_proteome == "yes") &
                       (T.feline_medulla_proteome == "yes")).sum())},
        {"step": "Detected in the Pod-TRECK proteome", "genes": n("pod_treck_proteome")},
        {"step": "Group A: detected in all three", "genes": int((T.group == "A").sum())},
        {"step": "Group A with finite values in all 16 states",
         "genes": int(((T.group == "A") & (T.complete_in_16_states == "yes")).sum())},
        {"step": "Group B: the remainder of the Δ matrix", "genes": int((T.group == "B").sum())},
    ])
    S.to_csv(OUT / "groupA_detection_summary.tsv", sep="\t", index=False)
    print(S.to_string(index=False))
    print(f"\n書き出し: {OUT / 'groupA_detection_calls.tsv'}（{len(T)} 行）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
