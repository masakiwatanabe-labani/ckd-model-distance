# P. 先行研究の確認（Discussion 執筆の前提）

調査日 2026-09-05。各項目で「先行あり/なし」と差分を明記する。
**逐語引用が取れたものと、検索要約からの記述は区別して書く。**
「先行なし」は「検索で見つからなかった」であって不在の証明ではない。

---

## P1. 射影型指標 — **指標そのものは先行あり／反例の指摘は見つからず**

### 先行あり: Signature Projection Score（SPS）

US 7,519,519 B1 "Signature projection score"（発明者 Eynon BP, Natsoulis G, Jarnigan K、
出願人 Entelos Inc）。条件 c に対する SPS は

    SPS(c) = Σ_g (X_gc − R_g)(T_g − R_g) / S_g

で、T は処置効果、R は対照/参照状態、S は尺度。これは対照を原点として、
参照ベクトル (T − R) 方向への射影であり、**本稿の α と同じ族**である。

**振幅を含めることを明示的に利点として主張している（Google Patents から逐語取得）:**

> "This method is advantageous over the correlation coefficient method because it takes into
> account the amplitude of the expression changes as well as their direction; it is superior
> to Euclidian distance methods because it measures direction and magnitude relative to the
> untreated or control state"

正規化についても記載がある:

> "the sum over genes of the product of these two differences is computed, and normalized such
> that the highest score of any of the characterizing treatments when calculated in this way
> is 100."

**つまり正規化は「特徴づけに使った処置群の最大値」を 100 とするもので、
参照自身を上限とする正規化ではない。**新しい条件が 100 を超えうる構造だが、
取得したテキストにその議論は見当たらない（自動読み取りによる否定なので、
「見当たらない」であって「無い」ではない）。

同族特許 US 7,396,645（Cholestasis signature）、US 7,778,782（PPARα signatures）も
同じ SPS を用いる。両者の USPTO PDF は画像ベースで、本調査では本文を取得できなかった。

### 先行あり: CMap / LINCS

Connectivity Map の connectivity score は双方向の重み付き Kolmogorov–Smirnov 統計量
（ES）に基づき、CMap 2.0 の WCS、LINCS の tau へと発展した。**これらは順位ベースで、
参照集合に対して正規化される**ため、SPS のような非有界の射影とは設計が異なる。
再現性の批判（Lim & Pavlidis; Sci Rep 2021, s41598-021-97005-z）や
複数スコアの整合（Brief Bioinform 2021, bbab161）はあるが、
**いずれも「参照自身を超える値を取りうる」という問題設定ではない。**

### 先行あり: 射影型のクロス種フレームワーク

TransComp-R（Translatable Components Regression）は動物データで PCA を組み、
ヒトデータをその潜在空間へ**射影**して回帰する。骨関節症（Osteoarthritis and
Cartilage 2025, S1063458425011574 / bioRxiv 2025.02.23.639777）、炎症性腸疾患、
神経変性へ適用例がある。射影の値が similarity として読めるかという議論は見当たらない。

### 先行なし（見つからず）: 「参照自身を超えうる」ことの指摘

single-sample scoring の限界に関する批判文献は多数ある（eLife 2022, 71994;
Brief Bioinform 2025, bbaf684; Foroutan et al. singscore, BMC Bioinformatics 2018）が、
いずれも **スコアの分布安定性・サンプル構成依存・二重利用（double dipping）**の議論で、
**射影スコアが自己参照値を超えうるため similarity として解釈できない、という指摘は
見つからなかった。**

### → §2.1 の書き換え

「距離が方向と振幅を混ぜている」という書き方をやめる（L で確認したとおり
1 − Spearman はスケール不変で、そもそも方向のみの量。この表現は誤り）。
代わりに **SPS が振幅包含を意図的な利点として主張していることを引用し、
その帰結として「similarity として解釈できない」反例を示す**構成にする。

---

## P2. 経路集約と種差 — **論法の先行は明確にあり**

### 中心的な例: COVID-19 マウスモデル

Bishop CR, Dumenil T, Rawle DJ, Le TT, Yan K, Tang B, Hartel G, Suhrbier A.
*PLoS Pathog* 2022;18(9):e1010867. doi:10.1371/journal.ppat.1010867

> "Overlap of single copy orthologue differentially expressed genes (scoDEGs) between human
> and mouse studies was generally poor (≈15–35%)."

> "As immunity and immunopathology are the focus of most studies, these mouse models can thus
> be viewed as representative and relevant models of COVID-19."

Discussion でも同旨:

> "All three mouse models showed poor scoDEG overlap with four studies of infected human lung
> tissues. The overlaps ranged from 21–35% for up-regulated scoDEGs and 15–27% for
> down-regulated scoDEGs."
> "As dominant pathways in immunity and inflammation are the target of most COVID-19
> interventions, these mouse models can be viewed as providing representative and pertinent
> models for pre-clinical assessments of new interventions."

**「遺伝子レベルは不一致だが経路レベルは一致 → モデルは妥当」という論法そのもの。**

### 同種の論法の他例

- **敗血症/Pathprint**: 経路レベル解析が遺伝子レベルより創薬リードの生成に有効であり、
  マウス内毒素血症モデルで検証されたとする報告（PMC5974511）。
- **結核**: "Concordant and discordant gene expression patterns in mouse strains identify
  best-fit animal model for human tuberculosis"（PMC5608750）。
- **CAMO**（下記 P3/制約参照）も経路中心の congruence を提案する立場。

### → 我々の位置づけ（Discussion の独立段落）

我々は**同じ集約操作を実データで再現し、それが種差を捨てていることを示した**:
経路単位では 89% の経路が AUC ≥ 0.90 で猫とマウスを分離するのに、
各経路を1スコアに潰すと AUC 0.541 まで落ちる。
**経路レベルの一致は、種差が無いことの証拠ではなく、集約が種差を見えなくした結果でありうる。**
これは上記論法への直接の反証になる。ただし「経路レベル解析が有用でない」とは主張しない
（用途によっては集約が適切）。主張は「一致の根拠として集約を使うのは循環的」に限定する。

---

## P3. 時間軸 — **直接の先行は見つからず／隣接研究はあり**

### 隣接あり: 種間の時系列アラインメント

- **TimeMeter**（Jiang et al., *Nucleic Acids Res* 2020;48(9):e51）— DTW で発現時系列を
  整列し、progression advance score (PAS) で「どちらが時間的に進んでいるか」を測る。
  ヒト–マウス神経分化、マウス指再生–アホロートルに適用。
- **Genes2Genes**（Sumanaweera et al., *Nat Methods* 2024, s41592-024-02378-4）—
  ベイズ情報理論による pseudotime 軌道の遺伝子レベル整列。
  in vitro T 細胞が "match an immature in vivo state while lacking expression of genes
  associated with TNF signaling" と、**in vivo からの乖離を特定**した点は概念的に近い。
- **TemporalVAE**（*Nat Cell Biol* 2025）— 発生時系列のアトラス補助マッピング。

### 隣接あり: 複数時点を持つモデル–ヒト対応

Pandey et al., *Alzheimers Dement* 2024, doi:10.1002/alz.087565 /
PMC11710784 — 「distinct mouse models match to distinct human AD subtypes in
age-dependent manner」。モデルごとに複数時点を取り、どの時点がどのヒト亜型に
最も合うかを見ている。ただし**取得した本文からは、「後の時点ほどヒトに近づくか」を
明示的に問うた記述と、その否定的結論は確認できなかった。**

### 先行なし（見つからず）: ペア間距離と固定参照軸への角度の分離

以下はいずれも複数の検索で見つからなかった:

1. 「時点を選び直してもモデルは疾患方向に近づかない」に相当する明示的な報告
2. **ペア間の相対的近さ**（隣接時点は互いに似る）と**固定参照軸への絶対角度**
   （疾患方向へ向かっているか）を**分離**した報告
3. 時間との相関を方向成分と振幅成分に分解した報告

DTW/pseudotime 系はいずれも**2軌道間の対応付け**が目的で、
「軌道が外部の固定軸へ向かっているか」は問うていない。

### → §2.3 が本稿で最も新規性が高い

ただし Genes2Genes の in vitro/in vivo 乖離の知見は概念的に近いので必ず引用し、
「軌道間の対応付けではなく、外部の固定疾患軸への角度を測る点が異なる」と差分を明示する。

---

## 制約: 「単一の最良モデルは存在しない」— **既報が複数あるため新規性を主張しない**

必ず引用する:

| 文献 | 内容 |
|---|---|
| **Zong W, et al. CAMO. *PNAS* 2023;120(6):e2202584120. PMID 36730203** | 経路中心の congruence 定量。"visually identify molecular mechanisms and pathway subnetworks that are **best or least mimicked** by model organisms"。**経路ごとにモデル適合度が違うことは CAMO が既に示している。** |
| **Lee HJ, et al. *eLife* 2022;11:e70763**（マラリア） | "Mouse models of infection will **not recapitulate all features** of the human response"。P. yoelii 17XL が重症マラリアに最も近いが、モデルごとに再現する側面が異なる。 |
| **Pandey et al. *Alzheimers Dement* 2024**（AD 亜型） | 異なるマウスモデルが異なるヒト AD 分子亜型に対応。年齢依存。 |
| **TransComp-R**（OA: *Osteoarthritis Cartilage* 2025 ほか） | DMM と ACLR で translatable な PC が異なる＝モデルごとに翻訳可能な生物学が違う。 |
| Seok J, et al. *PNAS* 2013;110:3507 / Takao & Miyakawa *PNAS* 2015;112:1167 | マウスモデルの妥当性をめぐる有名な対立。CAMO はこの論争への応答として書かれている。 |

**我々の貢献は手続き面に限定して書く:**

1. **Kendall の W = 0.362** — 順位の入れ替わりを1つの数値で定量した
2. **ランダム遺伝子集合の帯**（n=30 で 0.652、n=50 で 0.540）— 経路スコアがサイズだけで
   どれだけ動くかを示し、主張できる範囲を明示的に線引きした
3. **減衰の天井** — 観測値が測定限界か否かを分離した

### 訂正（CAMO 本文 PMC9963430 を取得後）

**当初「CAMO は測定信頼性の天井を明示していない」と書いたが、これは誤りだった。
CAMO は天井を明示的に実装している。**本文（自動読み取り）より:

> "assess the 'ceiling' of (i.e., the highest possible) concordance given the cohort,
> sample size, and experimental design"

intracohort c-score（同一コホートを分割して同種内の一致度を測る）を天井として用い、
"median c-score = 0.15 for HI1 vs. MI vs. 0.43 for HI1 vs. HI2" のように
種間の値を同種内の天井と並べて報告している。**天井という発想は我々の独自性ではない。**

残る差分は次の2点に限られる。

1. **経路サイズ由来のノイズ床の定量。** CAMO は問題を認識してはいる
   （"For pathway-specific analysis, however, variabilities of c-scores are much increased
   and results are less reliable."）が、対処は遺伝子集合サイズの下限・上限の設定
   （5-200 / 3-500）であり、**サイズの関数としてスコアがどれだけ偶然で動くかは
   定量していない。**我々はランダム遺伝子集合の 95% 帯（n=30 で 0.652、n=50 で 0.540）を
   出し、主張できる範囲を線引きしている。
2. **方向と振幅の分解。** CAMO の c-score / d-score は F 値ベースで、遺伝子ごとの
   符号の一致／不一致を 0-1 に正規化したもの。**ベクトルの向き（cos θ）と長さ（ν）の
   分解ではない。**

加えて **Kendall の W による順位不安定性の定量**も CAMO には該当記述がない。

天井については CAMO を先行として明示的に引用し、我々の split-half + Spearman-Brown は
**同じ発想の別実装**として位置づける。独自性は主張しない。

---

## 本文への反映（このあと実施）

- **§2.1**: SPS の逐語引用 → その帰結として反例、という構成に書き換え（P1）
- **Discussion 段落**: 集約論法への反証を独立段落に（P2）。PLoS Pathog を名指しで引用
- **Discussion**: §2.3 を最も新規性の高い部分として書き、Genes2Genes / TimeMeter との差分を明示（P3）
- **Discussion**: 「単一の最良モデルは存在しない」は既報として CAMO / eLife マラリア /
  AD 亜型 / TransComp-R を引用し、我々は手続き面の貢献のみ主張（制約）

## 未確認・限界

- US 7,396,645 と US 7,778,782 の本文（USPTO PDF が画像ベースで取得できず）
- CAMO 本文 → PMC9963430 から取得済み。上記の訂正を反映
- Pandey et al. の本文（abstract のみ確認）
- 「先行なし」は網羅的検索の結果であり、不在の証明ではない
