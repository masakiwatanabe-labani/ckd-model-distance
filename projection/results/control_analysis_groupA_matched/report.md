# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34

Group A intersection, 全16状態で complete case の 1632 遺伝子。
基準軸は固定していない（全 120 ペアの cos を対称に計算）。

## クラス別 cos 分布

| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |
|---|---|---|---|---|---|---|
| within_dataset | 45 | 0.157 | 0.454 | 0.744 | 0.825 | 0.957 |
| same_species_diff_dataset | 27 | -0.037 | 0.500 | 0.715 | 0.845 | 0.886 |
| cross_species | 48 | 0.152 | 0.319 | 0.394 | 0.454 | 0.534 |
| within_cat | 6 | 0.744 | 0.776 | 0.794 | 0.835 | 0.901 |
| within_podtreck | 3 | 0.771 | 0.779 | 0.786 | 0.871 | 0.957 |
| within_iri | 36 | 0.157 | 0.414 | 0.708 | 0.818 | 0.949 |
| shared_control_pairs | 33 | 0.157 | 0.422 | 0.729 | 0.853 | 0.957 |
| no_shared_control_pairs | 87 | -0.037 | 0.338 | 0.451 | 0.716 | 0.886 |

## 対照と検定対象

| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |
|---|---|---|---|---|---|
| 検定対象 | cat_CKD12 x cat_CKD34 | 0.901 | 0.886 | 91 | 94 |
| 上振れ対照(中央値) | IRI 時点ペア 36組 | 0.708 | - | 40 | 50 |
| 上振れ対照(最大) | IRI 時点ペア 36組 | 0.949 | - | 96 | 97 |
| 下振れ対照(指示どおり) | cat_CKD34 x cat_med_CKD34 | 0.843 | 0.779 | 76 | 78 |
| 下振れ対照(進行軸) | cat_ctx_prog x cat_med_prog | 0.004 | -0.093 | 0 | 0 |

## ネコ内ペアの内訳

| ペア | 種別 | cos | Spearman rho |
|---|---|---|---|
| cat_CKD12 x cat_CKD34 | same_compartment_diff_stage | 0.901 | 0.886 |
| cat_CKD12 x cat_med_CKD12 | same_stage_diff_compartment | 0.744 | 0.819 |
| cat_CKD12 x cat_med_CKD34 | diff_compartment_diff_stage | 0.778 | 0.795 |
| cat_CKD34 x cat_med_CKD12 | diff_compartment_diff_stage | 0.776 | 0.820 |
| cat_CKD34 x cat_med_CKD34 | same_stage_diff_compartment | 0.843 | 0.779 |
| cat_med_CKD12 x cat_med_CKD34 | same_compartment_diff_stage | 0.811 | 0.747 |

## 分割対照（共有対照バイアスの除去）

対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として
計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも
同じ処理を施してある（どちらも対照 n=6 -> 3）。

- cat_CKD12 x cat_CKD34: 共有対照 0.901 -> 分割対照 中央値 0.751 [95% 0.551, 0.849]（20 通り）
- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 0.128 / 中央値 0.655 / 最大 0.921
- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の 第 64 百分位

分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は
全体に下がる。両側に同じ処理をしているので比較は成立するが、
分割後の絶対値を共有対照の値と直接比べてはいけない。
