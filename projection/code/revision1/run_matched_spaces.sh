#!/bin/zsh
# Part 4. 対応表から定義し直した matched Group A / matched Group B の 2 空間を計算する。
set -e
cd "$(cd "$(dirname "$0")" && pwd)"
while [ ! -f make_figures.py ] && [ "$PWD" != / ]; do cd ..; done
V=../.venv/bin/python
for space in groupA_matched groupB_matched2; do
  echo "########## ${space}"
  export XSP_TAG="_${space}"
  export XSP_GENES="aa_${space}.txt"
  $V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
      --genes "$XSP_GENES" --time time_map.tsv --out "results/alpha_${space}" > /dev/null
  $V cos_matrix.py > /dev/null
  $V reliability.py > /dev/null
  $V ceiling_validation.py > /dev/null
  $V mantel_s1.py > /dev/null
  $V auc_permutation.py --pairs "results/control_analysis_${space}/pairs.tsv" \
      --out "results/control_analysis_${space}" > /dev/null
  $V pathway_separation.py > /dev/null
  echo "  done ${space}"
done
echo "matched spaces done"
