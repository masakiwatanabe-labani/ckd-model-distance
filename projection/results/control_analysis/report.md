# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34

Group A intersection, 全16状態で complete case の 2016 遺伝子。
基準軸は固定していない（全 120 ペアの cos を対称に計算）。

## クラス別 cos 分布

| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |
|---|---|---|---|---|---|---|
| within_dataset | 45 | 0.153 | 0.447 | 0.755 | 0.832 | 0.958 |
| same_species_diff_dataset | 27 | -0.033 | 0.502 | 0.707 | 0.840 | 0.884 |
| cross_species | 48 | 0.186 | 0.341 | 0.447 | 0.496 | 0.548 |
| within_cat | 6 | 0.759 | 0.791 | 0.813 | 0.851 | 0.907 |
| within_podtreck | 3 | 0.768 | 0.772 | 0.776 | 0.867 | 0.958 |
| within_iri | 36 | 0.153 | 0.410 | 0.703 | 0.813 | 0.950 |
| shared_control_pairs | 33 | 0.153 | 0.422 | 0.717 | 0.853 | 0.958 |
| no_shared_control_pairs | 87 | -0.033 | 0.375 | 0.494 | 0.711 | 0.884 |

## 対照と検定対象

| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |
|---|---|---|---|---|---|
| 検定対象 | cat_CKD12 x cat_CKD34 | 0.907 | 0.895 | 93 | 97 |
| 上振れ対照(中央値) | IRI 時点ペア 36組 | 0.703 | - | 40 | 50 |
| 上振れ対照(最大) | IRI 時点ペア 36組 | 0.950 | - | 96 | 97 |
| 下振れ対照(指示どおり) | cat_CKD34 x cat_med_CKD34 | 0.857 | 0.807 | 82 | 86 |
| 下振れ対照(進行軸) | cat_ctx_prog x cat_med_prog | -0.042 | -0.107 | 0 | 0 |

## ネコ内ペアの内訳

| ペア | 種別 | cos | Spearman rho |
|---|---|---|---|
| cat_CKD12 x cat_CKD34 | same_compartment_diff_stage | 0.907 | 0.895 |
| cat_CKD12 x cat_med_CKD12 | same_stage_diff_compartment | 0.759 | 0.822 |
| cat_CKD12 x cat_med_CKD34 | diff_compartment_diff_stage | 0.794 | 0.816 |
| cat_CKD34 x cat_med_CKD12 | diff_compartment_diff_stage | 0.790 | 0.831 |
| cat_CKD34 x cat_med_CKD34 | same_stage_diff_compartment | 0.857 | 0.807 |
| cat_med_CKD12 x cat_med_CKD34 | same_compartment_diff_stage | 0.832 | 0.774 |

## 分割対照（共有対照バイアスの除去）

対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として
計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも
同じ処理を施してある（どちらも対照 n=6 -> 3）。

- cat_CKD12 x cat_CKD34: 共有対照 0.907 -> 分割対照 中央値 0.772 [95% 0.553, 0.850]（20 通り）
- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 0.124 / 中央値 0.648 / 最大 0.923
- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の 第 64 百分位

分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は
全体に下がる。両側に同じ処理をしているので比較は成立するが、
分割後の絶対値を共有対照の値と直接比べてはいけない。
