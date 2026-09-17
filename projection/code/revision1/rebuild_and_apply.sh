#!/bin/zsh
# rev4 を正本に Markdown を作り直し、今ラウンドの変更をすべて当て直す。
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python
$V code/revision1/rebuild_markdown_from_rev4.py
for s in restore_formulas apply_task3_reactome apply_task1_gap apply_task2_sim apply_task4_filter apply_task5_move apply_partBC apply_reactome_figures apply_round4 apply_round5 apply_round6 apply_round7 apply_round8; do
  [ -f code/revision1/$s.py ] && $V code/revision1/$s.py
done
echo "--- 完了 ---"
