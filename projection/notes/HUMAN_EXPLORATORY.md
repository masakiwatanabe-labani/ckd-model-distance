# K. ヒト外挿（探索的） — 本採用の可否を決めるための報告

実行: `python human_extrapolation.py` → `results/human_exploratory/`
**主解析には組み込んでいません。** `delta_matrix.tsv` も `NUMBERS.md` も変更していません。

---

## 確認事項 1: split-half による天井が推定可能か

| ヒト状態 | プラットフォーム | 症例 n | 対照 n | n_max | r_half | 信頼性 r (SB) | 推定可能 |
|---|---|---|---|---|---|---|---|
| hs_ERCB_allCKD_A | GPL22945 | 39 | 18 | 9 | 0.888 | **0.941** | はい |
| hs_ERCB_RPGN_A | GPL22945 | 21 | 18 | 9 | 0.918 | **0.957** | はい |
| hs_ERCB_DN_A | GPL22945 | 7 | 18 | 3 | 0.763 | **0.865** | はい |
| hs_ERCB_allCKD_B | GPL24120 | 130 | 5 | 2 | 0.027 | **0.052** | 推定はできるが使えない |
| hs_ERCB_DN_B | GPL24120 | 10 | 5 | 2 | 0.356 | **0.525** | 同上 |
| hs_ERCB_HT_B | GPL24120 | 20 | 5 | 2 | −0.140 | **0.000** | 同上 |
| hs_KPMP_TI / hs_KPMP_G | KPMP-MS | — | — | 0 | — | — | **いいえ** |

**結論:**

- **GPL22945 は使えます。** 対照が living donor 18 検体あり、これは本研究で最も
  検体数の多い対照群です。信頼性 0.865–0.957 は、マウス側の多くの状態より高い水準です。
- **GPL24120 は使えません。** 対照が LD 3 + TN 2 = 5 検体しかなく、半分に割ると 2 検体に
  なるため、Δ が対照のサンプリング誤差に埋もれます（HT_B では r_half が負）。
  症例側が 130 検体あっても対照側が律速です。以降すべて除外しました。
- **KPMP は原理的に推定不可能です。** `DataLake_DEPs.txt` は事前計算済みコントラスト
  （Gene_name / Comparison / LogFC）のみで per-sample データがありません。
  加えてプロテオーム層であり、RNA 層の Δ 行列とは測定量が違います（本研究は層を混ぜない方針）。
  **制約どおり、KPMP の cos は主張に使いません。** 参考値のみ:
  cat × human 中央 0.384、mouse × human 中央 0.171。

## 確認事項 2: プラットフォーム差（マイクロアレイ vs RNA-seq）

**見積もれませんでした。** ERCB は2つのアレイ版に同じ疾患 (DN) を持つので、種・疾患・組織を
固定してプラットフォームだけ変えたペアが作れる設計でしたが、その相手が GPL24120 側であり、
信頼性が足りません。

| ペア | 内容 | cos | 天井 (SB) | 判定 |
|---|---|---|---|---|
| DN_A × DN_B | 同一疾患・同一組織・別アレイ版 | 0.776 | 0.674 | 天井超過。解釈不可 |
| allCKD_A × allCKD_B | 全CKD・別アレイ版 | 0.609 | 0.220 | 天井超過。解釈不可 |
| GPL22945 内のヒト状態ペア | 同一アレイ版・別疾患 | 0.930 | 0.910 | 同一プラットフォーム内の参照値 |

**ただし、この不確かさは判断を妨げません。** cat × human と mouse × human は**どちらも
RNA-seq 対 マイクロアレイ**であり、プラットフォーム差は2クラスに対称にかかります。
したがって「cat × human と mouse × human のどちらが大きいか」の比較はプラットフォームで
交絡しません。**交絡するのは cat × human と cat × mouse の絶対値比較のほう**で、
これは行いません。

## 確認事項 3: 3種共通空間の遺伝子数

| 集合 | n |
|---|---|
| Group A intersection（cat+mouse、全16状態 complete case） | 2,016 |
| ERCB の遺伝子（Entrez CDF） | 12,070 |
| **3種共通（Group A ∩ ERCB）** | **1,690** |
| KPMP のタンパク | 8,646 |
| （参考）Group A ∩ KPMP | 1,967 |

**十分です。** 1,690 遺伝子は主解析の 2,016 の 84% にあたり、この部分空間で
cat × mouse の cos 中央値は 0.437（主解析の全 2,016 遺伝子では 0.4475）と、
空間を絞ったことによる歪みはほぼありません。

---

## 結果（GPL22945 の3状態のみ）

| クラス | ペア数 | cos 中央 | cos 範囲 | 天井 中央 | 天井までの最小余裕 | 天井超過 |
|---|---|---|---|---|---|---|
| **cat × human** | 12 | **0.465** | 0.422–0.543 | 0.872 | 0.302 | 0 |
| **mouse × human** | 36 | **0.390** | 0.022–0.505 | 0.927 | 0.460 | 0 |
| cat × mouse（比較用、同一部分空間） | 48 | 0.437 | 0.171–0.543 | 0.862 | 0.307 | 0 |

**期待どおりの方向でした（cat × human > mouse × human）。**
どのクラスも天井から十分離れており、測定限界の記述ではありません。

状態レベル網羅並べ替え（§2.3 と同じ手続き。非ヒト16状態の種ラベルを入れ替え、
C(16,4)=1,820 通りを全列挙）:

| 部分集合 | 状態 | ペア | AUC | 帰無中央 | 厳密 p |
|---|---|---|---|---|---|
| 全マウス状態 | 16 | 48 | **0.8981** | 0.5012 | **0.0016** |
| IRI 2h を除く | 15 | 45 | 0.8889 | 0.5051 | 0.0022 |
| 急性期 (IRI 2h/4h) を除く | 14 | 42 | 0.8778 | 0.5028 | 0.0030 |
| 7 d 以上のみ | 11 | 33 | 0.9127 | 0.5000 | 0.0030 |
| 急性期 + Pod-TRECK D5 を除く | 13 | 39 | 0.8642 | 0.5000 | 0.0042 |

**急性期の時点を落としても結論は変わりません**（AUC 0.86–0.91、p ≤ 0.0042）。
「マウスパネルに急性状態が混ざっているから離れて見えるだけ」という説明は成り立ちません。

### ただし §2.3 とは性質が違う

状態ごとの cos 中央値（ヒト3状態に対して）:

- 猫 4 状態: 0.450 / 0.461 / 0.464 / 0.502 — 狭く固まる
- マウス最良: IRI_24h **0.457**、mouse_2W 0.440、IRI_7d 0.412
- マウス最低: IRI_2h 0.040、IRI_4h 0.241

**マウスの最良状態 (0.457) は猫の最低状態 (0.450) と重なります。**
§2.3 の種差は「異種は 0.548 を超えず、同種は 0.884 まで届く」という分離でしたが、
ここにあるのは分離ではなく**分布のシフト**です。主張の強さを §2.3 と同じにはできません。

---

## 判断材料のまとめ

**採用に有利:**

- 天井が推定でき（信頼性 0.865–0.957、マウス側より高い）、全ペアが天井から十分離れている
- 遺伝子数 1,690 は十分で、部分空間による歪みがない
- 核心の比較（cat × human vs mouse × human）はプラットフォーム差に対して対称
- 状態レベルの厳密検定で p = 0.0016、急性期を落としても p ≤ 0.0042
- 期待した方向

**採用に不利:**

- 使えるヒト状態は **1 コホート・1 プラットフォームの 3 状態のみ**
- **分布は重なる**（マウス最良 0.457 ≒ 猫最低 0.450）。§2.3 のような分離ではない
- KPMP は使えないので、ヒト側の独立検証がない
- プラットフォーム差の絶対量が不明（対称性で回避しているだけ）
- 種 × データセットの交絡は §2.3 と同じ構造で残る（データセットレベルなら p ≥ 1/3）

**当初の推奨は「§2.5 として採用」でしたが、下記の検査で撤回します。**

---

# 追加検査の結果 — 採用を撤回します

`python human_checks.py` → `results/human_exploratory/`

## (i) オルソログ写像の非対称性 — **交絡を排除できませんでした**

cat × human はネコ側だけが、mouse × human はマウス側だけがオルソログ写像を経ます。
写像の質は種間で非対称です:

- cat→human perc_id 中央 **89.0**、mouse→human **92.4**。マウスのほうが高い遺伝子が 55.3%

写像品質で層別すると、**種差は写像品質の関数として動き、最も写像が良い層では逆転します。**

| 層 | n | cat × human | mouse × human | cat × mouse | 差（cat−mouse） |
|---|---|---|---|---|---|
| 全遺伝子 | 1,690 | 0.465 | 0.390 | 0.437 | **+0.075** |
| **両種とも perc_id 上位50%** | 500 | 0.375 | **0.397** | 0.363 | **−0.021（逆転）** |
| 両種とも perc_id 下位50% | 500 | 0.497 | 0.429 | 0.462 | +0.068 |
| **cat の perc_id ≥ mouse** | 639 | 0.485 | 0.376 | 0.451 | **+0.110** |
| **cat の perc_id < mouse** | 792 | 0.465 | 0.439 | 0.461 | **+0.026** |

AUC と厳密 p でも同じ:

| 層 | n | AUC | p | 同サイズ乱択 AUC 中央 [5%点] | サイズで説明できるか |
|---|---|---|---|---|---|
| 全遺伝子 | 1,690 | 0.8981 | 0.0016 | — | — |
| 両種とも perc_id 上位50% | 500 | **0.4676** | 0.6011 | 0.8519 [0.5507] | **できない** |
| 両種とも perc_id 下位50% | 500 | 0.9398 | 0.0005 | 0.8588 [0.5716] | できる（典型的） |
| cat の perc_id ≥ mouse | 639 | 0.8981 | 0.0022 | 0.8773 [0.6502] | できる（典型的） |
| cat の perc_id < mouse | 792 | **0.6736** | 0.1324 | 0.8796 [0.7036] | **できない** |

**決定的なのは最後の2行です。** cat × mouse の cos は両層でほぼ同じ（0.451 と 0.461）で
信号量は変わらないのに、対ヒトの差だけが **+0.110 → +0.026** と約4分の1になります。
つまり「ネコがヒトに近い」の大きさは、**どちらの種の写像が良いかで決まっています。**
両種とも写像が良い遺伝子に絞ると差は消え（AUC 0.468、p=0.60）、
これは同サイズのランダム部分集合の分布（中央 0.852、5%点 0.551）から明らかに外れており、
遺伝子数を減らしたことによる検出力低下では説明できません。

**したがって「写像品質の非対称性」という代替仮説は棄却できず、むしろ支持されます。**
§2.3 では同じ検査で写像品質が棄却できた（同一性四分位が単調でない、高信頼度1:1に絞っても不変）
のに対し、ヒト軸では逆の結果になりました。§2.3 と §2.5 で結論が違うのは、
§2.3 の cat × mouse は両種とも写像を経る対称な比較であるのに対し、
§2.5 の cat × human と mouse × human は**片側だけが写像を経る非対称な比較**だからです。

## (ii) 状態の信頼性 — 交絡していません

信頼性 × 対ヒト cos の順位相関 rho = **+0.112**（16状態）。
猫4状態の信頼性中央 **0.838** に対しマウス12状態は **0.945** で、
マウスのほうが信頼性が高いのに対ヒト cos は低い。観測と逆向きなので、この交絡は排除できます。

## (iii) ヒト状態の入れ子と共有対照 — 部分的な問題

ヒト3状態は同じ LD 18 検体を対照とし、さらに allCKD_A は DN_A / RPGN_A を包含します。
分割対照でヒト内 cos は 0.930→0.877、0.982→0.932、0.884→0.844 と下がりますが、
ネコ／マウスの状態とヒトの状態は検体を共有しないため、cross-species ペアには効きません。

ただしヒト状態を1つずつ使うと結果が大きく振れます:

| 使ったヒト状態 | AUC | p |
|---|---|---|
| hs_ERCB_DN_A のみ | 1.0000 | 0.0005 |
| hs_ERCB_allCKD_A のみ | 0.8750 | 0.0148 |
| **hs_ERCB_RPGN_A のみ** | **0.7500** | **0.0852** |

RPGN（ANCA関連血管炎）単独では分離しません。ヒト側の疾患の選び方に依存します。

---

# 写像品質を揃えた再検定 — **成立せず（2026-09-05）**

`python human_mapping_matched.py` → `results/human_exploratory/mapping_matched_test.tsv`

層別ではなく**遺伝子ごとの caliper マッチング**で組み直した。
|perc_id(cat→human) − perc_id(mouse→human)| <= delta を満たす遺伝子だけを残すので、
集合内の全遺伝子で両種の写像品質が同程度になる（層別では「高cat/低mouse」と
「低cat/高mouse」が層内で相殺されてしまう）。

| caliper | n | pid_cat | pid_mouse | cat × human | mouse × human | 差 | cat × mouse | AUC | p |
|---|---|---|---|---|---|---|---|---|---|
| なし | 1,431 | 85.4 | 89.6 | 0.472 | 0.418 | **+0.054** | 0.465 | 0.833 | 0.0055 |
| ≤ 10 | 1,054 | 89.9 | 90.0 | 0.477 | 0.437 | +0.040 | 0.460 | 0.745 | 0.054 |
| ≤ 5 | 812 | 91.6 | 91.2 | 0.491 | 0.464 | +0.027 | 0.458 | 0.660 | 0.154 |
| **≤ 2** | 466 | **93.2** | **93.0** | 0.421 | **0.478** | **−0.058** | 0.457 | **0.331** | **0.874** |
| ≤ 5 かつ perc_id 上位半分 | 601 | 95.4 | 94.9 | 0.447 | 0.408 | +0.039 | 0.408 | 0.704 | 0.086 |

**交絡の用量反応がそのまま出ています。** caliper を狭める（写像品質を揃える）ほど差は
単調に縮み、最も揃った caliper ≤ 2（pid_cat 93.2 対 pid_mouse 93.0）で**符号が反転**します
（mouse × human 0.478 > cat × human 0.421）。

検出力低下ではありません。同サイズのランダム部分集合 500 回の AUC 分布は
90% 区間 [0.537, 0.961]（中央 0.806）で、**観測 0.331 はその下側に外れています。**
また cat × mouse はどの部分集合でも 0.457–0.465 とほぼ一定で、信号量自体は落ちていません。
動いているのは cat と mouse の「対ヒトの順位」だけです。

**したがって、写像品質を揃えた設計でも §2.5 は成立しません。**
「ネコがヒトに近い」という観測は、ネコとマウスのオルソログ写像品質の非対称性で
説明されます。KPMP snRNA-seq（M で存在を確認）を足しても、これはヒトデータセットの
性質ではなくオルソログ表の性質なので解消しません。

---

# 最終判断: **Limitations 記載で確定**

採用を支持していた材料（AUC 0.898、p=0.0016、急性期を落としても頑健、天井から十分離れている）は
いずれも有効ですが、**(i) の写像非対称性を排除できない以上、種の効果として報告できません。**
「ネコがヒトに近い」の効果量は、ネコとマウスのどちらの写像が良いかで +0.110 から +0.026 まで動き、
両種とも写像が良い遺伝子では逆転します。

§2.3（cat × mouse）は両側が写像を経る対称な比較なのでこの問題を受けませんが、
§2.5（cat × human vs mouse × human）は片側だけが写像を経る非対称な比較であり、
構造的にこの交絡を持ちます。**遺伝子ごとに写像品質を揃えた設計**（例: cat と mouse の
perc_id が同程度の遺伝子だけを使う、あるいは3種同時アラインメントで
オルソログ群を定義し直す）を組まない限り、この主張は立ちません。

## 本文 Limitations 用の英文（そのまま使える版）

> **Extrapolation to human CKD.** We examined whether the feline states sit closer to
> human CKD than the mouse states do, using the ERCB tubulointerstitial transcriptome
> (GSE104954) and KPMP. Of the two ERCB array platforms, only GPL22945 has a control
> group large enough to estimate a reliability ceiling (18 living-donor biopsies; the
> GPL24120 control group of five gives split-half reliabilities of 0.00-0.53). The KPMP
> regional proteomics release provides only pre-computed contrasts, so no ceiling can be
> estimated for it, and it is a proteome layer.
>
> On the three usable human states the feline states were closer to human CKD than the
> mouse states (median cos θ 0.465 versus 0.390; AUC 0.898, exact p = 0.0016 by
> state-level permutation), and the result was robust to excluding the acute IRI time
> points. We nonetheless do not report it, because it is confounded with an asymmetry in
> orthologue mapping that this design cannot remove: cos θ for cat × human depends on cat
> to human orthology alone, and cos θ for mouse × human on mouse to human orthology
> alone, and the two are not of equal quality (median percent identity 89.0 versus 92.4).
> Matching genes pairwise on mapping quality reverses the sign of the difference
> (|Δ percent identity| ≤ 2, n = 466, cat 93.2 versus mouse 93.0: median cos θ 0.421 for
> cat × human versus 0.478 for mouse × human, AUC 0.331), and the difference shrinks
> monotonically as the matching caliper is tightened from 10 to 5 to 2 (+0.040, +0.027,
> −0.058). The reversal is not a loss of power: an AUC of 0.331 falls below the 90%
> interval of equally sized random gene subsets ([0.537, 0.961]), and cos θ between cat
> and mouse states is unchanged across these subsets (0.457-0.465).
>
> This is a limitation of the comparison, not of the human data available. KPMP releases
> open-access per-participant single-nucleus expression matrices for 98 CKD and 68
> healthy-reference participants, so additional human cohorts could be added; but the
> asymmetry is a property of the orthologue tables rather than of any human dataset, and
> adding cohorts would not remove it. A design in which both species are mapped to human
> through a common procedure of equal quality — for example a three-species simultaneous
> alignment — would be required before this comparison could be made.

（対応する日本語の要点は `FINDINGS.md` の Limitations に記載）

## Limitations に書く内容（案・初版）

- ヒト外挿を ERCB (GSE104954) と KPMP regional proteomics で探索的に検討した。
- KPMP は事前計算済みコントラストのみで per-sample データがなく信頼性の天井を推定できず、
  かつプロテオーム層であるため、本研究の RNA 層の枠組みには組み込めない。
- ERCB は2つのアレイ版のうち GPL22945 のみが天井を推定できる対照群サイズ（LD 18検体）を持つ。
  GPL24120 は対照 5 検体で信頼性が 0.00–0.53 と推定不能水準。
- 使える3状態で cat × human 0.465 対 mouse × human 0.390 と期待方向の差が出たが、
  この差はオルソログ写像品質の種間非対称性（cat→human 89.0 対 mouse→human 92.4）と
  交絡しており、両種とも写像品質が高い遺伝子に限ると差は消える（AUC 0.468）。
  片側のみが写像を経る比較では、この交絡を現在のデータで分離できない。
- したがってヒトへの外挿は本研究の主張に含めない。
