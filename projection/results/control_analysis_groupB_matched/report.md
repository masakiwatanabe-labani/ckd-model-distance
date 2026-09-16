# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34

Group A intersection, 全16状態で complete case の 1618 遺伝子。
基準軸は固定していない（全 120 ペアの cos を対称に計算）。

## クラス別 cos 分布

| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |
|---|---|---|---|---|---|---|
| within_dataset | 45 | 0.161 | 0.354 | 0.623 | 0.757 | 0.945 |
| same_species_diff_dataset | 27 | 0.053 | 0.299 | 0.534 | 0.718 | 0.788 |
| cross_species | 48 | 0.124 | 0.260 | 0.327 | 0.385 | 0.453 |
| within_cat | 6 | 0.707 | 0.722 | 0.759 | 0.801 | 0.881 |
| within_podtreck | 3 | 0.623 | 0.628 | 0.632 | 0.789 | 0.945 |
| within_iri | 36 | 0.161 | 0.334 | 0.571 | 0.713 | 0.891 |
| shared_control_pairs | 33 | 0.161 | 0.347 | 0.621 | 0.762 | 0.945 |
| no_shared_control_pairs | 87 | 0.053 | 0.273 | 0.356 | 0.539 | 0.813 |

## 対照と検定対象

| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |
|---|---|---|---|---|---|
| 検定対象 | cat_CKD12 x cat_CKD34 | 0.881 | 0.842 | 93 | 97 |
| 上振れ対照(中央値) | IRI 時点ペア 36組 | 0.571 | - | 40 | 50 |
| 上振れ対照(最大) | IRI 時点ペア 36組 | 0.891 | - | 96 | 97 |
| 下振れ対照(指示どおり) | cat_CKD34 x cat_med_CKD34 | 0.813 | 0.714 | 91 | 97 |
| 下振れ対照(進行軸) | cat_ctx_prog x cat_med_prog | 0.082 | -0.025 | 0 | 0 |

## ネコ内ペアの内訳

| ペア | 種別 | cos | Spearman rho |
|---|---|---|---|
| cat_CKD12 x cat_CKD34 | same_compartment_diff_stage | 0.881 | 0.842 |
| cat_CKD12 x cat_med_CKD12 | same_stage_diff_compartment | 0.711 | 0.743 |
| cat_CKD12 x cat_med_CKD34 | diff_compartment_diff_stage | 0.766 | 0.742 |
| cat_CKD34 x cat_med_CKD12 | diff_compartment_diff_stage | 0.707 | 0.728 |
| cat_CKD34 x cat_med_CKD34 | same_stage_diff_compartment | 0.813 | 0.714 |
| cat_med_CKD12 x cat_med_CKD34 | same_compartment_diff_stage | 0.753 | 0.677 |

## 分割対照（共有対照バイアスの除去）

対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として
計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも
同じ処理を施してある（どちらも対照 n=6 -> 3）。

- cat_CKD12 x cat_CKD34: 共有対照 0.881 -> 分割対照 中央値 0.736 [95% 0.531, 0.833]（20 通り）
- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 0.111 / 中央値 0.519 / 最大 0.824
- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の 第 96 百分位

分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は
全体に下がる。両側に同じ処理をしているので比較は成立するが、
分割後の絶対値を共有対照の値と直接比べてはいけない。
