# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34

Group A intersection, 全16状態で complete case の 7897 遺伝子。
基準軸は固定していない（全 120 ペアの cos を対称に計算）。

## クラス別 cos 分布

| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |
|---|---|---|---|---|---|---|
| within_dataset | 45 | 0.181 | 0.372 | 0.680 | 0.789 | 0.954 |
| same_species_diff_dataset | 27 | 0.065 | 0.387 | 0.645 | 0.763 | 0.836 |
| cross_species | 48 | 0.179 | 0.337 | 0.377 | 0.456 | 0.572 |
| within_cat | 6 | 0.706 | 0.737 | 0.761 | 0.797 | 0.898 |
| within_podtreck | 3 | 0.692 | 0.694 | 0.696 | 0.825 | 0.954 |
| within_iri | 36 | 0.181 | 0.346 | 0.657 | 0.779 | 0.924 |
| shared_control_pairs | 33 | 0.181 | 0.361 | 0.668 | 0.820 | 0.954 |
| no_shared_control_pairs | 87 | 0.065 | 0.353 | 0.436 | 0.646 | 0.836 |

## 対照と検定対象

| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |
|---|---|---|---|---|---|
| 検定対象 | cat_CKD12 x cat_CKD34 | 0.898 | 0.871 | 93 | 97 |
| 上振れ対照(中央値) | IRI 時点ペア 36組 | 0.657 | - | 40 | 50 |
| 上振れ対照(最大) | IRI 時点ペア 36組 | 0.924 | - | 96 | 97 |
| 下振れ対照(指示どおり) | cat_CKD34 x cat_med_CKD34 | 0.809 | 0.728 | 78 | 81 |
| 下振れ対照(進行軸) | cat_ctx_prog x cat_med_prog | 0.109 | -0.022 | 0 | 0 |

## ネコ内ペアの内訳

| ペア | 種別 | cos | Spearman rho |
|---|---|---|---|
| cat_CKD12 x cat_CKD34 | same_compartment_diff_stage | 0.898 | 0.871 |
| cat_CKD12 x cat_med_CKD12 | same_stage_diff_compartment | 0.706 | 0.760 |
| cat_CKD12 x cat_med_CKD34 | diff_compartment_diff_stage | 0.760 | 0.755 |
| cat_CKD34 x cat_med_CKD12 | diff_compartment_diff_stage | 0.730 | 0.755 |
| cat_CKD34 x cat_med_CKD34 | same_stage_diff_compartment | 0.809 | 0.728 |
| cat_med_CKD12 x cat_med_CKD34 | same_compartment_diff_stage | 0.762 | 0.692 |

## 分割対照（共有対照バイアスの除去）

対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として
計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも
同じ処理を施してある（どちらも対照 n=6 -> 3）。

- cat_CKD12 x cat_CKD34: 共有対照 0.898 -> 分割対照 中央値 0.762 [95% 0.589, 0.843]（20 通り）
- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 0.127 / 中央値 0.585 / 最大 0.872
- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の 第 75 百分位

分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は
全体に下がる。両側に同じ処理をしているので比較は成立するが、
分割後の絶対値を共有対照の値と直接比べてはいけない。
