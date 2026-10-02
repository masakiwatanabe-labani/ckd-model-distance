#!/bin/zsh
# BE. S4 の段落を設計別の値中心に直し、図に両設計を描き、提出一式を作り直す。
# 回答書は verify_numbers.py が書いた検査件数を読むので、最後に 2 周させる。
set -e
cd "$(cd "$(dirname "$0")" && pwd)/../.."
V=../.venv/bin/python
R=code/revision1

echo "### 1. シミュレーションの超過率を設計別に集計（BE1）"
$V $R/be_simulation_summary.py | tail -4

echo "\n### 2. 補足図を描き直し（S4 に両設計）、描画データを TSV に（BE2）"
$V $R/bd_figures.py | tail -14

echo "\n### 3. 図↔本文の検査"
$V $R/bd_check_figure_text.py | tail -2

echo "\n### 3b. 文献の追加（Crossref 照会）と番号の振り直し（BI）"
$V $R/bi_add_references.py | tail -4
$V $R/renumber_refs_C.py | tail -2
$V $R/ay_check_citations.py | tail -1
$V $R/bj_check_refnums.py | tail -1

echo "\n### 4. 本文 docx と構造検査"
$V $R/an_build_from_base.py | tail -3
$V $R/an_check_docx.py | tail -2

echo "\n### 5. 補足（図を埋めたノート・S1-S5 順の図 PDF・表）"
$V $R/build_supplement.py | tail -4

echo "\n### 6. 変更履歴版（check_tracked.py を通す）と青字版"
$V $R/at_package_tracked.py | tail -4
$V $R/al_check_tracked.py | tail -1
$V $R/bh_package_blue.py | tail -4
$V $R/bh_check_blue.py | tail -2
$V $R/al_metadata_audit.py | tail -1

echo "\n### 7. NUMBERS.md → 回答書 → マニフェスト（件数を合わせるため 2 周）"
$V $R/build_numbers_review.py | tail -1
for i in 1 2; do
  $V verify_numbers.py | tail -1
  $V $R/build_response_letter.py | tail -1
  $V $R/build_letters_docx.py > /dev/null
  $V $R/al_manifest.py > /dev/null
done
$V verify_numbers.py | tail -2
echo "\n--- BE 完了 ---"
