"""図 3・4 の入力を Reactome 込みの結果に向け、§2.3 の四条件の値を 422 セットのものに直す。

図 3B だけは意図的に従来3コレクションの 74 セットに固定する（Section 4.13）。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R = 'manuscript/RESULTS.md'

rp(R, "Averaging the centred values into 422 pathway means raises both medians to 0.869 "
      "(0.791-0.934) and 0.882 (0.792-0.927), and the AUC falls to 0.511.",
      "Averaging the centred values into 422 pathway means raises both medians to 0.831 "
      "(0.682-0.902) and 0.832 (0.683-0.917), and the AUC falls to 0.511.")
rp(R, "Taking pathway means of the uncentred Δ leaves within-species pairs at 0.896 (0.829-0.941) "
      "and cross-species pairs at 0.489 (0.355-0.590), separating at AUC 0.899.",
      "Taking pathway means of the uncentred Δ leaves within-species pairs at 0.870 (0.699-0.937) "
      "and cross-species pairs at 0.462 (0.314-0.597), separating at AUC 0.868.")
rp('manuscript/DISCUSSION.md',
   "With centring, within- and cross-species pairs separate at AUC 0.511, which is chance, both "
   "classes having risen to a median of about 0.87. Without it, the same pathway means separate at "
   "0.899, better than the 0.813 of the gene-level comparison they replace. Cross-species pairs are "
   "the ones that move: their median goes from 0.448 to 0.882 when the mean is removed first, and "
   "to 0.489 when it is not.",
   "With centring, within- and cross-species pairs separate at AUC 0.511, which is chance, both "
   "classes having risen to a median of about 0.83. Without it, the same pathway means separate at "
   "0.868, better than the 0.813 of the gene-level comparison they replace. Cross-species pairs are "
   "the ones that move: their median goes from 0.448 to 0.832 when the mean is removed first, and "
   "to 0.462 when it is not.")
rp('manuscript/FRONTMATTER.md',
   "Third, averaging genes within pathways separates the species at AUC 0.899 or at 0.511 according "
   "to whether each state's mean change is removed first",
   "Third, averaging genes within pathways separates the species at AUC 0.868 or at 0.511 according "
   "to whether each state's mean change is removed first")
rp('manuscript/DISCUSSION.md',
   "AUC 0.899 therefore does not mean the uncentred comparison is closer to the truth.",
   "AUC 0.868 therefore does not mean the uncentred comparison is closer to the truth.")
print('Reactome figure values applied')
