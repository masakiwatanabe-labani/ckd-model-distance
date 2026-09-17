# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34

Group A intersection, 全16状態で complete case の 1632 遺伝子。
基準軸は固定していない（全 120 ペアの cos を対称に計算）。

## クラス別 cos 分布

| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |
|---|---|---|---|---|---|---|
| within_dataset | 45 | 0.166 | 0.350 | 0.633 | 0.763 | 0.949 |
| same_species_diff_dataset | 27 | 0.085 | 0.317 | 0.589 | 0.721 | 0.807 |
| cross_species | 48 | 0.137 | 0.263 | 0.326 | 0.363 | 0.438 |
| within_cat | 6 | 0.723 | 0.741 | 0.764 | 0.809 | 0.876 |
| within_podtreck | 3 | 0.653 | 0.658 | 0.662 | 0.805 | 0.949 |
| within_iri | 36 | 0.166 | 0.318 | 0.575 | 0.736 | 0.894 |
| shared_control_pairs | 33 | 0.166 | 0.330 | 0.633 | 0.774 | 0.949 |
| no_shared_control_pairs | 87 | 0.085 | 0.281 | 0.358 | 0.579 | 0.824 |

## 対照と検定対象

| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |
|---|---|---|---|---|---|
| 検定対象 | cat_CKD12 x cat_CKD34 | 0.876 | 0.844 | 93 | 97 |
| 上振れ対照(中央値) | IRI 時点ペア 36組 | 0.575 | - | 40 | 50 |
| 上振れ対照(最大) | IRI 時点ペア 36組 | 0.894 | - | 96 | 97 |
| 下振れ対照(指示どおり) | cat_CKD34 x cat_med_CKD34 | 0.824 | 0.736 | 91 | 97 |
| 下振れ対照(進行軸) | cat_ctx_prog x cat_med_prog | 0.044 | -0.070 | 0 | 0 |

## ネコ内ペアの内訳

| ペア | 種別 | cos | Spearman rho |
|---|---|---|---|
| cat_CKD12 x cat_CKD34 | same_compartment_diff_stage | 0.876 | 0.844 |
| cat_CKD12 x cat_med_CKD12 | same_stage_diff_compartment | 0.723 | 0.747 |
| cat_CKD12 x cat_med_CKD34 | diff_compartment_diff_stage | 0.763 | 0.756 |
| cat_CKD34 x cat_med_CKD12 | diff_compartment_diff_stage | 0.734 | 0.758 |
| cat_CKD34 x cat_med_CKD34 | same_stage_diff_compartment | 0.824 | 0.736 |
| cat_med_CKD12 x cat_med_CKD34 | same_compartment_diff_stage | 0.764 | 0.699 |

## 分割対照（共有対照バイアスの除去）

対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として
計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも
同じ処理を施してある（どちらも対照 n=6 -> 3）。

- cat_CKD12 x cat_CKD34: 共有対照 0.876 -> 分割対照 中央値 0.735 [95% 0.543, 0.829]（20 通り）
- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 0.121 / 中央値 0.527 / 最大 0.835
- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の 第 86 百分位

分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は
全体に下がる。両側に同じ処理をしているので比較は成立するが、
分割後の絶対値を共有対照の値と直接比べてはいけない。
