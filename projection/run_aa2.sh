#!/bin/zsh
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python

echo "### ortholog_all: pathway_separation のみ再開"
XSP_TAG="_ortholog_all" XSP_GENES="aa_ortholog_all.txt" $V pathway_separation.py > /dev/null
echo "  ok"

echo "### groupB_matched: 全ステップ"
export XSP_TAG="_groupB_matched" XSP_GENES="aa_groupB_matched.txt"
$V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
    --genes "$XSP_GENES" --time time_map.tsv --out results/alpha_groupB_matched > /dev/null
echo "  alpha ok"
$V cos_matrix.py > /dev/null;            echo "  cos ok"
$V reliability.py > /dev/null;           echo "  reliability ok"
$V mantel_s1.py > /dev/null;             echo "  mantel ok"
$V auc_permutation.py --pairs results/control_analysis_groupB_matched/pairs.tsv \
    --out results/control_analysis_groupB_matched > /dev/null; echo "  auc ok"
$V pathway_separation.py > /dev/null;    echo "  pathway ok"
echo "### AA 完了"
