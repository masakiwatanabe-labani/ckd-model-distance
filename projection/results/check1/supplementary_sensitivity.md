# Supplementary: Δ の作り方に対する感度分析（cat_CKD34 x mouse_2W）

主解析は S1（ortholog / FPKM>=1 / log2FC）に固定する。

| 変種 | 写像 | 発現フィルタ | Δ の定義 | n | Spearman rho | 主解析との差 | 主解析 |
|---|---|---|---|---|---|---|---|
| S1 | ortholog | on | log2FC | 8815 | +0.4893 | +0.0000 | **はい** |
| S2 | ortholog | on | z-scored | 8815 | +0.2575 | -0.2317 |  |
| S3 | ortholog | off | log2FC | 10499 | +0.3878 | -0.1014 |  |
| S4 | ortholog | off | z-scored | 10499 | +0.2263 | -0.2630 |  |
| S5 | naive | on | log2FC | 10483 | +0.4998 | +0.0105 |  |
| S6 | naive | on | z-scored | 10483 | +0.2715 | -0.2178 |  |
| S7 | naive | off | log2FC | 12438 | +0.3982 | -0.0911 |  |
| S8 | naive | off | z-scored | 12438 | +0.2377 | -0.2515 |  |

## 参照点

- 主解析 S1（ortholog / filter on / log2FC）: rho = +0.4893（旧稿 0.490）
- naive / filter off / log2FC:  rho = +0.3982（新稿 0.398）
- naive / filter off / z-scored: rho = +0.2377
- S1 との差: -0.0911 / -0.2515

## 要因ごとの主効果（他2要因で平均した差）

- delta: -0.1955
- fpkm_filter: -0.0670
- mapping: +0.0116
