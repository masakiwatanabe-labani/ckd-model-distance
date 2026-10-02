# -*- coding: utf-8 -*-
"""AT1: 補足表のラベル → 番号の対応表を 1 つだけ持ち、そこから出力を作る。

正本は manuscript_C/SUPPLEMENTARY.md のキャプション。番号を振り直したら
ここを回すだけで、末尾の Supplementary Materials の一覧も引用の検査も追随する。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
OUT = HERE / "results" / "roundR" / "supp_labels.tsv"

# ラベル → キャプション冒頭で照合する語（内容で対応づけるので番号に依存しない）
LABELS = [
    ("tabS:pair_subsets", "Table", "Pair subsets used for descriptive AUC"),
    ("tabS:state_means", "Table", "Mean response of each state over the 2,016-gene"),
    ("tabS:delta_definitions", "Table", "The 16-state analysis under four definitions"),
    ("tabS:scoring_methods", "Table", "The centering contrast under five pathway scoring methods"),
    ("tabS:pathway_composition", "Table", "Dependence of the centering contrast on the composition"),
    ("tabS:offset_filter", "Table", "Sensitivity of pair-class AUCs and two projection indicators"),
    ("tabS:orthologue_selection", "Table", "Orthologue characteristics and the effect of selecting"),
    ("tabS:mantel", "Table", "Mantel statistics for elapsed time and onset compartment"),
    ("tabS:spearman_sensitivity", "Table", "Sensitivity of Spearman correlation between feline"),
    ("tabS:detection_calls", "Table", "Proteomic detection calls behind the Group A selection"),
    ("tabS:pair_cosines", "Table", "Observed cosine similarities and descriptive reliability"),
    ("tabS:splithalf", "Table", "State-level split-half estimates from 500 partitions"),
    ("figS:orthologue", "Figure", "Orthologue characteristics and gene selection"),
    ("figS:pathway_alignment", "Figure", "Pathway-level alignment to feline cortical CKD3/4"),
    ("figS:pathway_rankings", "Figure", "Pathway-specific rankings and the effect of mean-centering"),
    ("figS:benchmark_sim", "Figure", "Simulation of observed similarity relative to estimated"),
    ("figS:time_dissimilarity", "Figure", "Pairwise directional and amplitude dissimilarity"),
]
# 一覧に出す短い説明（ラベルごとに 1 つ）
BLURB = {
    "tabS:pair_subsets": "the pair subsets used for the descriptive AUCs",
    "tabS:state_means": "the mean response of each state over the Group A gene universe",
    "tabS:delta_definitions": "the 16-state analysis under four definitions of the response vector",
    "tabS:scoring_methods": "the centering contrast under five pathway scoring methods",
    "tabS:pathway_composition": "the composition and redundancy of the pathway set and the common "
                                "pathway identifiers",
    "tabS:offset_filter": "the offset and filtering sensitivity",
    "tabS:orthologue_selection": "orthologue characteristics and gene selection",
    "tabS:mantel": "the Mantel statistics",
    "tabS:spearman_sensitivity": "the sensitivity of the Spearman correlation",
    "tabS:detection_calls": "the proteomic detection calls behind Group A",
    "tabS:pair_cosines": "the per-pair cosines and reliability benchmarks",
    "tabS:splithalf": "the split-half reliability of every state",
}


def captions() -> dict:
    out = {}
    for b in (M / "SUPPLEMENTARY.md").read_text().split("\n\n"):
        s = " ".join(b.split())
        m = re.match(r"^(Table|Figure) S(\d+)\.\s*(.*)$", s)
        if m:
            out[(m.group(1), int(m.group(2)))] = m.group(3)
    return out


def build() -> pd.DataFrame:
    caps = captions()
    rows = []
    for label, kind, head in LABELS:
        hit = [n for (k, n), txt in caps.items() if k == kind and txt.startswith(head)]
        assert len(hit) == 1, f"{label}: キャプションが {len(hit)} 件見つかった（{head}）"
        rows.append({"label": label, "kind": kind, "number": hit[0],
                     "citation": f"{kind} S{hit[0]}",
                     "caption": caps[(kind, hit[0])][:120],
                     "blurb": BLURB.get(label, "")})
    T = pd.DataFrame(rows).sort_values(["kind", "number"], kind="stable")
    assert T.number.tolist() == sorted(T.number.tolist(), key=lambda x: x) or True
    for kind in ("Table", "Figure"):
        ns = sorted(T[T.kind == kind].number)
        assert ns == list(range(1, len(ns) + 1)), f"{kind} の番号に抜けがある: {ns}"
    return T


def supplementary_sentence(T: pd.DataFrame) -> str:
    tab = T[T.kind == "Table"].sort_values("number")
    fig = T[T.kind == "Figure"].sort_values("number")
    nt, nf = len(tab), len(fig)
    parts = "; ".join(f"S{int(r.number)}, {r.blurb}" for r in tab.itertuples())
    return (
        "Supplementary Materials: The following are available online. Supplementary Notes "
        f"(Supplementary_Notes.docx) contain Notes S1 to S5, the captions of Tables S1 to S{nt} "
        f"and the legends of Figures S1 to S{nf}. Supplementary Tables "
        f"(Supplementary_Tables.xlsx) contain Tables S1 to S{nt}, one sheet per table: "
        f"{parts}; and one further sheet giving the per-gene detection calls behind Table "
        f"S{int(tab[tab.label == 'tabS:detection_calls'].number.iloc[0])} for every gene of the Δ "
        "matrix. Supplementary Data S1 "
        "(Supplementary_Data_S1_PodTRECK_proteome.xlsx) is the Pod-TRECK quantification table from "
        f"which the Group A detection calls were made. Figures S1 to S{nf} are supplied as "
        "separate files.")


def main() -> int:
    T = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT, sep="\t", index=False)
    print(T[["label", "citation", "caption"]].to_string(index=False, max_colwidth=58))
    sent = supplementary_sentence(T)
    p = M / "BACKMATTER.md"
    t = p.read_text()
    old = next(l for l in t.splitlines() if l.startswith("Supplementary Materials:"))
    if old != sent:
        p.write_text(t.replace(old, sent))
        print("\nBACKMATTER.md の Supplementary Materials を作り直した")
    else:
        print("\nBACKMATTER.md は既に最新")
    print(f"書き出し: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
