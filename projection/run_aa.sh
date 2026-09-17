#!/bin/zsh
# AA. 遺伝子空間の感度解析。既定空間（Group A）の再実行で回帰がないことを確かめてから、
# 全 1:1 orthologue と発現量マッチ Group B で同じ 4 量を再計算する。
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python

echo "########## 0. 既定空間の再実行（回帰確認）"
$V cos_matrix.py > /dev/null
$V reliability.py > /dev/null
$V mantel_s1.py > /dev/null
$V pathway_separation.py > /dev/null
$V verify_numbers.py | tail -3

for space in ortholog_all groupB_matched; do
  echo "########## ${space}"
  export XSP_TAG="_${space}"
  export XSP_GENES="aa_${space}.txt"
  $V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
      --genes "$XSP_GENES" --time time_map.tsv --out "results/alpha_${space}" > /dev/null
  $V cos_matrix.py > /dev/null
  $V reliability.py > /dev/null
  $V mantel_s1.py > /dev/null
  $V auc_permutation.py --pairs "results/control_analysis_${space}/pairs.tsv" \
      --out "results/control_analysis_${space}" > /dev/null
  $V pathway_separation.py > /dev/null
  echo "  done ${space}"
done
echo "########## AA 完了"
