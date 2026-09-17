#!/bin/zsh
# Part F-1. matched Group A（1,578 遺伝子）空間で Table 4 の量を計算する。
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python
export XSP_TAG="_groupA_matched"
export XSP_GENES="aa_groupA_matched.txt"
$V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
    --genes "$XSP_GENES" --time time_map.tsv --out "results/alpha_groupA_matched" > /dev/null
echo "  alpha ok"
$V cos_matrix.py > /dev/null;        echo "  cos_matrix ok"
$V reliability.py > /dev/null;       echo "  reliability ok"
$V mantel_s1.py > /dev/null;         echo "  mantel ok"
$V auc_permutation.py --pairs "results/control_analysis_groupA_matched/pairs.tsv" \
    --out "results/control_analysis_groupA_matched" > /dev/null; echo "  auc_perm ok"
$V pathway_separation.py > /dev/null; echo "  pathway ok"
echo "F-1 gene space done"
