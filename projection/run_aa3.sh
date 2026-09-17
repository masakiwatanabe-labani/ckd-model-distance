#!/bin/zsh
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python
while pgrep -f run_aa2.sh > /dev/null; do sleep 10; done
echo "### run_aa2 完了。天井検証を 3 空間で実行"
$V ceiling_validation.py > /dev/null;                       echo "  groupA ok"
$V verify_numbers.py | tail -1
XSP_TAG="_ortholog_all"   XSP_GENES="aa_ortholog_all.txt"   $V ceiling_validation.py > /dev/null; echo "  ortholog_all ok"
XSP_TAG="_groupB_matched" XSP_GENES="aa_groupB_matched.txt" $V ceiling_validation.py > /dev/null; echo "  groupB_matched ok"
echo "### 集計"
$V aa_summary.py
echo "### 全完了"
