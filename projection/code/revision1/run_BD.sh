#!/bin/zsh
# BD. 補足図を現在の結果から描き直し、図↔本文の検査を通してから提出一式を作り直す。
# 回答書は verify_numbers.py が書いた検査件数を読むので、この順で 2 周する。
set -e
cd "$(cd "$(dirname "$0")" && pwd)/../.."
V=../.venv/bin/python
R=code/revision1

echo "### 1. 補足図を描き直し、描画データを TSV に落とす（BD2）"
$V $R/bd_figures.py | tail -14

echo "\n### 2. 図↔本文の検査（BD3）"
$V $R/bd_check_figure_text.py | tail -2

echo "\n### 3. 本文 docx を土台から作り直す"
$V $R/an_build_from_base.py | tail -3
$V $R/an_check_docx.py | tail -2

echo "\n### 4. 補足（図を埋め込んだノート・図 PDF・表）"
$V $R/build_supplement.py | tail -4

echo "\n### 5. 変更履歴版"
$V $R/at_package_tracked.py | tail -3

echo "\n### 6. 検査 → 回答書 → 検査（件数を合わせるため 2 周）"
for i in 1 2; do
  $V verify_numbers.py | tail -1
  $V $R/build_response_letter.py | tail -1
  $V $R/build_letters_docx.py > /dev/null
  $V $R/al_manifest.py > /dev/null
done
$V verify_numbers.py | tail -2
echo "\n--- BD 完了 ---"
