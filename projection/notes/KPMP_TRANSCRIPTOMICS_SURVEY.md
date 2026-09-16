# M. KPMP トランスクリプトームの per-sample データ調査（解析なし）

調査日 2026-09-05。KPMP Atlas Repository の検索 API（Elastic App Search、
`atlas.kpmp.org/spatial-viewer/search/api/as/v1/engines/atlas-repository/search.json`、
公開検索キーはフロントエンドのバンドルに埋め込まれているもの）へのカタログ照会のみ。
**データのダウンロードも解析も行っていません。**

---

## 1. per-sample 発現行列は存在するか → **モダリティによる**

KPMP のアクセス階層は「open = 正規化済み／非正規化の発現行列と実験メタデータ」、
「controlled = 臨床・病理データと配列レベル（FASTQ/BAM）」。
つまり**発現行列は原則 open**（ただし利用規約への同意が必要。既に取得済みの
proteomics と同じ扱い）。

| モダリティ | open ファイル数 | 形式 | per-sample 発現行列 |
|---|---|---|---|
| Single-nucleus RNA-Seq | 552 | tsv mtx 244 / xlsx 272 / h5Seurat 4 | **あり** |
| Single-cell RNA-Seq | 656 | tsv mtx 326 / xlsx 326 / h5Seurat 4 | **あり** |
| 10x Multiome | 305 | h5 tsv mtx 152 / xlsx 153 | **あり** |
| **Regional Transcriptomics (LMD)** | **39** | csv 20 / xlsx 18 / txt 1 | **実質なし**（下記） |
| Bulk RNA-Seq | **0** | — | **なし** |

**Regional Transcriptomics（我々が使った LMD proteomics と同じ区画分離）は使えません。**
全 319 ファイルのうち 280 が controlled（fastq 150 / bam 130）で、open な 39 件の内訳は
per-participant の生カウント csv が 20 件、メタデータ xlsx が 18 件、
Atlas Explorer 用の集約 zip が 1 件です。内部の事前調査どおり、Atlas Explorer 版は
「区画特異性スコア」であって CKD vs 対照ではありません。

## 2. 群サイズと split-half の可否

open access の per-participant 発現行列を持つ**参加者数**:

| モダリティ | CKD | Healthy Reference | split-half の可否 |
|---|---|---|---|
| **Single-nucleus RNA-Seq** | **98** | **68** | **可能。半分でも 49 対 34** |
| Single-cell RNA-Seq | 169 | 60 | 可能 |
| **Regional Transcriptomics** | **8** | **1** | **不可能（対照が1名）** |

（参考: 全体の参加者登録数は AKI 92 / CKD 524 / DM-R 67、Healthy Reference の
付随研究は 125 名で 2023-06-30 に終了。上表は「open な発現行列が実際にある人数」。）

**snRNA-seq / scRNA-seq なら天井は十分推定できます。**対照 68 名は本研究で最大の
対照群（現行最大は ERCB GPL22945 の living donor 18 名）で、split-half の
不等長分割の問題も起きません。

**Regional Transcriptomics は Healthy Reference が1名で、対照群を構成できません。**
天井以前に Δ が定義できません。

## 3. ERCB GPL22945 との共通遺伝子数

**未計算です**（ダウンロードしていないため）。ただし上限は構造的に決まります。

- 現行の3種共通空間 = Group A ∩ ERCB = **1,690 遺伝子**
- ERCB GPL22945 は Brainarray ENTREZG カスタム CDF で **12,070 プローブセット**
- KPMP snRNA-seq は 10x の標準アノテーションで蛋白コード領域をほぼ網羅する

したがって律速は **Group A（2,016）と ERCB（12,070）であって KPMP ではありません。**
KPMP snRNA-seq を足しても3種共通空間は 1,690 前後から大きく変わらないと見込まれます
（参考: 既に計算済みの Group A ∩ KPMP proteomics は 1,967）。
確定値が要るなら 1 参加者分の行列を取得すれば足ります。

---

## 判断への含意

**ヒト側は ERCB 1コホートが上限、という前提は誤りでした。**
KPMP snRNA-seq（CKD 98 / Healthy Reference 68）が open access で存在し、
天井も推定できます。「単一ヒトコホート」という Limitations の書き方は、
このままでは不正確です。

**ただし、これは §2.5 を復活させません。**

§2.5 を止めた理由は**ヒトコホートの数ではなく、オルソログ写像の非対称性**でした
（cat→human と mouse→human で片側だけが写像を経る。両種とも perc_id 上位50%の
遺伝子に絞ると AUC 0.898 → 0.468 と消え、cat の perc_id が mouse 以上か否かで
効果量が +0.110 → +0.026 と動く）。これは**オルソログ表の性質であって
ヒトデータセットの性質ではない**ため、ヒトのコホートを増やしても解消しません。

§2.5 を成立させるには、依然として写像品質を揃えた設計（perc_id を揃えた遺伝子集合、
または3種同時アラインメントによるオルソログ群の再定義）が必要です。
KPMP snRNA-seq の追加は、その設計を採った場合の**独立検証コホート**としては
極めて有用ですが、それ単独では交絡を解きません。

## 追加で考慮すべき点（採用する場合）

- snRNA-seq は疑似バルク化が必要で、bulk の log2FC とは正規化が異なります。
  内部の解析メモも「ヒトを正面から比較するなら snRNA-seq の疑似バルクのほうが素直」
  としており、方向としては整合します。
- Healthy Reference には living kidney donor / kidney stone donor / pilot nephrectomy が
  混在します。ERCB の対照（living donor）と完全には揃いません。
- 取得量は snRNA-seq の発現行列 244 パッケージ。Mount Sinai の Data Ark ミラーは
  AIR·MS の DUA が別途必要です。KPMP 本体からの取得なら既存の同意で足ります。

## 出典

- KPMP Available Data（アクセス階層とモダリティ別の形式）https://www.kpmp.org/available-data
- KPMP Controlled access data https://www.kpmp.org/controlled-data
- KPMP Atlas Repository 検索 API（ファセット集計。上表の数値の出所）
- Atlas Explorer v1.0 Regional Transcriptomics RNA-seq data, doi:10.48698/t9fh-qn48
