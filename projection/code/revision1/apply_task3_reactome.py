"""Task 3 の本文反映（Reactome）。rev4 由来の Markdown に当てる。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'

rp(R, "every GO Biological Process, KEGG and Hallmark set with at least 30 Group A genes "
      "(193 sets: 118, 48 and 27 respectively; Reactome was not available locally). Each of the "
      "16 states was projected onto the reference axis using only that pathway's genes, giving "
      "3,088 cells (Figure 3B).",
      "every Reactome, GO Biological Process, KEGG and Hallmark set with at least 30 Group A "
      "genes (422 sets: 229, 118, 48 and 27 respectively; three names occur in two collections "
      "with different membership and are kept separate). Each of the 16 states was projected "
      "onto the reference axis using only that pathway's genes, giving 6,752 cells (Figure 3B).")
rp(R, "Pathway sizes here are median 44 (interquartile range 34-61), and 119 of the 193 sets fall "
      "below 50 genes.",
      "Pathway sizes here are median 45 (interquartile range 33-67), and 244 of the 422 sets fall "
      "below 50 genes.")
rp(R, "Thirty-four of the 2,316 mouse × pathway cells reach 0.9 or above",
      "Forty-five of the 5,064 mouse × pathway cells reach 0.9 or above")
rp(R, "Of the 888 cells meeting the size condition", "Of the 2,136 cells meeting the size condition")
rp(R, "Kendall's W across the 193 pathways is 0.362, and the median rank correlation between two "
      "pathways is +0.315",
      "Kendall's W across the 422 pathways is 0.320, and the median rank correlation between two "
      "pathways is +0.336")
rp(R, "IRI 12 mo leads 25%, IRI 28 d 20%, IRI 72 h 11%, and the remaining nine states share the "
      "rest. IRI 12 mo has a mid-ranked whole-panel cos θ of 0.442 yet leads the most pathways.",
      "IRI 28 d leads 22%, IRI 12 mo 21%, IRI 2 h 11%, and nine further states share the remainder. "
      "IRI 28 d has a mid-ranked whole-panel cos θ of 0.484 yet leads the most pathways, and IRI 2 h, "
      "the lowest-aligned state of the panel at 0.193, leads 11%.")
rp(R, "0.788 in the median (interquartile range 0.733-0.839), with 33% of pathways at or above the "
      "whole-panel value, 13 of 193 reaching 0.90 and 15 falling to 0.60 or below.",
      "0.760 in the median (interquartile range 0.707-0.808), with 23% of pathways at or above the "
      "whole-panel value, 16 of 422 reaching 0.90 and 25 falling to 0.60 or below.")
rp(R, "the same 87 pairs compared in the resulting 193-dimensional space",
      "the same 87 pairs compared in the resulting 422-dimensional space")
rp(R, "Averaging the centred values into 193 pathway means",
      "Averaging the centred values into 422 pathway means")
rp(R, "its variation is not a size artefact (rank correlation with gene count +0.035; Figure 4C): "
      "the lowest",
      "we do not read its variation as biology. Set size accounts for part of it - the rank "
      "correlation with gene count is +0.198 (Figure 4C), against +0.035 before Reactome was added - "
      "and what the earlier collections made look like a size-independent spread does not survive a "
      "wider set of collections. The per-set claims we do make are therefore gated on the "
      "size-matched null of each set rather than on the spread across sets. Reporting the ordering "
      "for context: the lowest")
rp(R, "pathway resolution separates better than aggregation in every space (0.788 vs 0.511; 0.788 "
      "vs 0.633; 0.754 vs 0.608)",
      "pathway resolution separates better than aggregation in every space (0.760 vs 0.511; 0.782 "
      "vs 0.600; 0.771 vs 0.626)")
rp(R, "| Pathway-level AUC, median | 0.788 | 0.788 | 0.754 |",
      "| Pathway-level AUC, median | 0.760 | 0.782 | 0.771 |")
rp(R, "| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.633 | 0.608 |",
      "| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.626 |")
rp(R, "119 of 193 pathways fall below 50 genes, where that range exceeds 0.54",
      "244 of 422 pathways fall below 50 genes, where that range exceeds 0.54")
rp(R, "within each of the 193 pathways", "within each of the 422 pathways")

rp(D, "with 193 sets and no per-set test", "with 422 sets and no per-set test")
rp(D, "the interquartile range across the 193 pathways is 0.733 to 0.839 and 15 of",
      "the interquartile range across the 422 pathways is 0.707 to 0.808 and 25 of")
rp(D, "Kendall's W across 193 pathways is 0.362", "Kendall's W across 422 pathways is 0.320")
rp(D, "and 119 of our 193 pathways fall below 50 genes", "and 244 of our 422 pathways fall below 50 genes")
rp(D, "a median AUC of 0.788 with 15 of 193 pathways at 0.60 or below",
      "a median AUC of 0.760 with 25 of 422 pathways at 0.60 or below")
rp(D, "Reactome was not included. ", "")

rp(M, "Gene sets were taken from MSigDB Hallmark [31], the Gene Ontology Biological Process "
      "collection [32] and KEGG [33], and restricted to those with at least 30 Group A genes, "
      "giving 193 sets: 118 GO Biological Process, 48 KEGG and 27 Hallmark. Reactome was not "
      "available locally and is not included.",
      "Gene sets were taken from Reactome, MSigDB Hallmark [31], the Gene Ontology Biological "
      "Process collection [32] and KEGG [33], and restricted to those with at least 30 Group A "
      "genes, giving 422 sets: 229 Reactome, 118 GO Biological Process, 48 KEGG and 27 Hallmark. "
      "Three set names occur in two collections with different membership and are kept as separate "
      "sets. The collections are included because each is independently curated and in general use, "
      "so that the pathway-level result does not rest on one curator's choice of sets; no collection "
      "was selected or dropped on the number or size of the sets it contributes. They are not "
      "independent of one another: the median Reactome set shares a maximum Jaccard index of 0.467 "
      "with its closest KEGG set and 49% of Reactome sets reach 0.5, so the 422 sets are not 422 "
      "independent tests, and we make no claim that any set is non-redundant. Adding Reactome moves "
      "the per-set results - the median area under the curve falls from 0.788 to 0.760 and its rank "
      "correlation with set size rises from +0.035 to +0.198 - while leaving the comparison between "
      "pathway resolution and aggregation unchanged (Section 2.3).")
rp(M, "the same 87 pairs in the 193-dimensional space", "the same 87 pairs in the 422-dimensional space")
rp(M, "a vector of 193 pathway means", "a vector of 422 pathway means")
print('Task 3 applied')
