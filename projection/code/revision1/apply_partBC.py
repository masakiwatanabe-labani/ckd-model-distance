"""Part B と Part C の本文反映。

B-1 Figure 3B は従来3コレクションの 74 セットに固定し、そのことをキャプションと
     Section 4.13 に明記する。
B-2 IRI 2 h のリード率 11% の内訳（47/48 が Reactome の小さな重複セット）を本文に添える。
B-3 経路間順位相関の 5-95 パーセンタイルを 422 セットで再計算した値に直す。
C   天井の過小評価が信頼性の低い側で大きいため、群間差 +0.219 は下限である旨を 3.2 に加える。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'

# ---- B-3: 5-95 パーセンタイルを 422 セットの値に ----
rp(R, "the median rank correlation between two pathways is +0.336 (5th-95th percentile -0.374 to "
      "+0.846)",
      "the median rank correlation between two pathways is +0.336 (5th-95th percentile -0.385 to "
      "+0.888)")

# ---- B-2: IRI 2 h のリード率の内訳 ----
rp(R, "and IRI 2 h, the lowest-aligned state of the panel at 0.193, leads 11%.",
      "and IRI 2 h, the lowest-aligned state of the panel at 0.193, leads 11%. That last figure is "
      "an artefact of set redundancy rather than a broad signal: 47 of the 48 sets IRI 2 h leads are "
      "Reactome sets, and the ten highest are proteasomal and ubiquitin-degradation sets of 30 to 33 "
      "genes that share most of their members. We note it and draw nothing from it.")

# ---- B-1: Figure 3B を 74 セットに固定することを明記 ----
rp(R, "cos θ to the reference axis for all 16 states across the 74 pathways with at least 50 Group A "
      "genes, ordered by mean mouse cos θ.",
      "cos θ to the reference axis for all 16 states across the 74 pathways of the three original "
      "collections (GO Biological Process, KEGG, Hallmark) with at least 50 Group A genes, ordered "
      "by mean mouse cos θ. The panel is deliberately held to those three collections; adding the "
      "104 Reactome sets that also pass the size criterion would widen it to 178 columns, most of "
      "them near-duplicates of one another (Section 4.13).")
rp(M, "No per-pathway test was performed and no pathway is described as significant.",
      "No per-pathway test was performed and no pathway is described as significant. Figure 3B, "
      "which displays cells rather than summarising them, is held to the 74 sets of the three "
      "original collections that pass the size criterion rather than to all 178: the Reactome sets "
      "that would be added overlap one another heavily, so the extra columns would not add "
      "independent information and would make the panel unreadable. Every number reported in "
      "Section 2.3 uses all 422 sets.")

# ---- Part C ----
rp(D, "The cross-species shortfall is continuous with the same-species one, and the evidence "
      "supports a difference of degree and nothing stronger.",
      "The cross-species shortfall is continuous with the same-species one, and the evidence "
      "supports a difference of degree and nothing stronger. That difference is itself a lower "
      "bound. The simulation of Section 4.10 shows the ceiling estimate sits below the value two "
      "identical states attain, and the shortfall of the estimate grows as reliability falls; the "
      "uncorrected ceilings of the cat-mouse pairs are lower than those of the same-species pairs "
      "(median 0.765 against 0.852 and 0.896 in Table 1), so their gaps are the more strongly "
      "understated of the two. The +0.219 we report is therefore smaller than the difference would "
      "be if both groups were measured equally well.")
print('Part B and C applied')
