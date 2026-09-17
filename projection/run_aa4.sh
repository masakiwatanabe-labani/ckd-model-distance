#!/bin/zsh
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python
L=/private/tmp/claude-501/-Users-masaki/f7fba45f-f8e1-4d4a-b8f2-da5834e0fbe6/scratchpad/aa4.log
: > $L

echo "### groupB_matched の残り" >> $L
export XSP_TAG="_groupB_matched" XSP_GENES="aa_groupB_matched.txt"
$V reliability.py > /dev/null;        echo "  reliability ok" >> $L
$V mantel_s1.py > /dev/null;          echo "  mantel ok" >> $L
$V auc_permutation.py --pairs results/control_analysis_groupB_matched/pairs.tsv \
    --out results/control_analysis_groupB_matched > /dev/null; echo "  auc ok" >> $L
$V pathway_separation.py > /dev/null; echo "  pathway ok" >> $L
$V ceiling_validation.py > /dev/null; echo "  ceiling ok" >> $L

echo "### 天井検証（残り 2 空間）" >> $L
unset XSP_TAG XSP_GENES
$V ceiling_validation.py > /dev/null; echo "  groupA ok" >> $L
$V verify_numbers.py | tail -1 >> $L
XSP_TAG="_ortholog_all" XSP_GENES="aa_ortholog_all.txt" $V ceiling_validation.py > /dev/null
echo "  ortholog_all ok" >> $L

echo "### 集計" >> $L
$V aa_summary.py >> $L 2>&1
echo "### 全完了" >> $L
cat $L
