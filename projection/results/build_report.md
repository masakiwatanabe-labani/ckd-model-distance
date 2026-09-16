# Step 0 build report

- Δ行列: 22345 遺伝子 x 16 状態
- 全状態で有限な遺伝子（complete case）: 9480
- core（全データセットで Ensembl ortholog 由来）: 12291
- added（どこかで大文字フォールバック）: 10054

## Group A

| 集合 | n（Δ行列内） | n（写像後の素の集合） |
|---|---|---|
| intersection: (cortex ∩ medulla) ∩ mouse prot | 2240 | 2252 |
| union: (cortex ∪ medulla) ∩ mouse prot | 2964 | 2982 |

中間量: cat_cortex_detected=2927, cat_medulla_detected=2652, podtreck_prot_detected=12301, cat_ctx_and_med=2393, cat_ctx_or_med=3186

## 状態あたりの有限値数

| 状態 | 有限値数 |
|---|---|
| cat_CKD12 | 15053 |
| cat_CKD34 | 15053 |
| cat_med_CKD12 | 15170 |
| cat_med_CKD34 | 15170 |
| mouse_5D | 15213 |
| mouse_2W | 15213 |
| mouse_3W | 15213 |
| IRI_2h | 13484 |
| IRI_4h | 13696 |
| IRI_24h | 14149 |
| IRI_48h | 14066 |
| IRI_72h | 14155 |
| IRI_7d | 14274 |
| IRI_14d | 14225 |
| IRI_28d | 14237 |
| IRI_12mo | 13580 |
