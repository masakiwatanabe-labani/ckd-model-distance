#!/bin/zsh
# C案の Markdown を docx から作り直し、第C ラウンドの修正をすべて当てる。
# 途中でやり直せるように、この順で通しで走らせる。
set -e
cd "$(cd "$(dirname "$0")" && pwd)/../.."
V=../.venv/bin/python
R=code/revision1

echo "### 0. docx → Markdown"
rm -rf manuscript_C
$V $R/docx_to_markdown_C.py "$1"

echo "### 1. 対照解析（Part 2）を計算し本文へ"
$V $R/centering_controls.py
$V $R/apply_C_part2.py
$V $R/apply_C_part2b.py

echo "### 2. Part 0・3・4・5"
$V $R/apply_C_parts.py
$V $R/apply_C_part5b.py
$V $R/apply_C_refs_eq.py

echo "### 3. 文献の振り直しと仕上げ"
$V $R/renumber_refs_C.py
$V $R/apply_C_final_polish.py

echo "### 4. 図"
$V $R/make_figures_C.py > /dev/null 2>&1
echo "  図 8 点"

echo "### 5. 検査"
$V verify_numbers.py | tail -3
