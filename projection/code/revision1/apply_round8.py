"""第8ラウンド（最終校正）の本文反映。新しい解析も新しい数値も作らない。

A  Table 4 を 1 セル 1 値に作り直し、§2.6 を表と一致させる
B  Figure 2 のキャプションをパネルの中身に合わせる
C  §3.2・§3.8・§4.10・§4.12 を §2.2 の方針にそろえる
D  中心化効果の個体再標本化を Results（§2.3）で初出にし、Methods に定義を書く
E  数式の位置、重複文、Table S1 の凡例、作業用注記
"""
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'
C, F, S = 'manuscript/CONCLUSIONS.md', 'manuscript/FRONTMATTER.md', 'manuscript/SUPPLEMENTARY.md'
BM = 'manuscript/BACKMATTER.md'

# =====================================================================
# A-3. Table 4 を 1 セル 1 値の版に差し替える（build_table4.py が結果ファイルから生成）。
# =====================================================================
def _table4():
    new = (_P('results') / 'round8' / 'table4.md').read_text().strip()
    p = _P(R)
    t = p.read_text()
    head = "**Table 4.**"
    i = t.index(head)
    j = t.index("\n\n", i) + 2                      # キャプションの直後
    k = t.index("\n\n", j)                          # 表ブロックの終わり
    old = t[j:k]
    assert old.lstrip().startswith("|"), "Table 4 の表ブロックが見つからない"
    p.write_text(t[:j] + new + t[k:])
    print(f'A-3: Table 4 を {len(new.splitlines()) - 2} 行 x 4 列の 1 セル 1 値版に差し替え')


_table4()

rp(R, "**Table 4.** *Principal results recomputed in four gene spaces. Group A is the main "
      "analysis. The two matched columns are the two sides of one saved 1:1 matching on "
      "control-group expression rank, so they hold the same 1,632 genes each and are balanced to "
      "a standardised mean difference of −0.017, against +1.236 between the unmatched sets; the "
      "384 Group A genes with no partner are the high-expression ones (Section 4.14).",
      "**Table 4.** *Principal results recomputed in four gene spaces. Group A is the main "
      "analysis. The two matched columns are the two sides of one saved 1:1 matching on "
      "control-group expression rank, so they hold the same 1,632 genes each and are balanced to "
      "a standardised mean difference of −0.017, against +1.236 between the unmatched sets; the "
      "384 Group A genes with no partner are the high-expression ones (Section 4.14). Each cell "
      "holds one value, so that α, cos θ and ν and the point estimate and its interval are on "
      "separate rows.")

# =====================================================================
# A-1. §2.6 の残った旧値を表に合わせる。
# =====================================================================
rp(R, "And pathway resolution separates better than aggregation in every space (0.760 vs 0.511; "
      "0.782 vs 0.600; 0.771 vs 0.626; 0.779 vs 0.542).",
      "And pathway resolution separates better than aggregation in every space (0.760 vs 0.511; "
      "0.782 vs 0.600; 0.780 vs 0.624; 0.781 vs 0.529).")
rp(R, "First, α exceeds the reference’s own value of 1.000 only in Group A, and there by a point "
      "estimate whose interval includes 1.000 (1.020, 0.965–1.079, against 0.796 and 0.930).",
      "First, α exceeds the reference’s own value of 1.000 only in Group A, and there by a point "
      "estimate whose interval includes 1.000 (1.020, 0.965–1.079, against 0.796, 0.987 and "
      "0.947).")
print('A-1: §2.6 の旧値を Table 4 に合わせた')

# =====================================================================
# B. Figure 2 のキャプションをパネルの中身に合わせる。
# 図は第6ラウンドで B = Spearman-Brown、C = 無補正に入れ替えたが、キャプションが
# 旧配置のままだった。丸囲みは B に付く。D は補正後の差を描いている。
# =====================================================================
rp(R, "(B, C) Observed cos θ for all 120 state pairs against the reliability-derived benchmark "
      "for that pair, under the uncorrected and the Spearman–Brown- corrected benchmark. The "
      "diagonal is the benchmark and the shaded region above it is above the estimated benchmark "
      "under the independent-error model. Cross-species pairs lie far below their benchmarks "
      "under both estimators. The four pairs above the benchmark in (C) are circled; each shares "
      "a cohort or its control samples.",
      "(B, C) Observed cos θ for all 120 state pairs against the reliability-derived benchmark "
      "for that pair: (B) the Spearman–Brown benchmark, which is the primary form, and (C) the "
      "uncorrected one. The diagonal is the benchmark and the shaded region above it is above the "
      "estimated benchmark under the independent-error model; the arrow marks the smallest gap a "
      "cross-species pair leaves, 0.297 in (B) and 0.171 in (C). Cross-species pairs lie far "
      "below their benchmarks under both. The four pairs above the Spearman–Brown benchmark are "
      "circled in (B); each shares a cohort or its control samples, which we report without "
      "reading a cause from it (Section 2.2). Nine pairs lie above the uncorrected benchmark in "
      "(C), which sits lower.")
rp(R, "(D) The same comparison with animals resampled instead of genes: for each pair the "
      "observed cos θ and the uncorrected benchmark are recomputed from the same 1,000 resamples "
      "of the animals",
      "(D) The same comparison with animals resampled instead of genes: for each pair the "
      "observed cos θ and the Spearman–Brown benchmark are recomputed from the same 1,000 "
      "resamples of the animals")
rp(R, "A range containing zero does not mean the pair reached its benchmark.",
      "The three horizontal lines are the group medians, 0.196, 0.487 and 0.338, and they are "
      "what the section reports; the ranges convey spread.")
print('B: Figure 2 のキャプションをパネルに合わせた')

# =====================================================================
# C. §3.2・§3.8・§4.10・§4.12 を §2.2 の方針（判定を撤回し、効果量で述べる）にそろえる。
# =====================================================================
rp(D, "co-vary, is less conclusive: the difference stays above zero in 26 of the 48 cross-species "
      "pairs. The split follows the split-half reliability of the feline state — the 24 pairs "
      "built on cortical CKD3/4 or medullary CKD1/2 (0.841, 0.851) all hold, while 22 of the 24 "
      "built on cortical CKD1/2 or medullary CKD3/4 (0.568, 0.616) cannot be resolved. Those 24 "
      "pairs share two feline states between them, so they are not 24 independent replications, "
      "and the smallest of their lower bounds is 0.003.",
      "co-vary, is less conclusive, and we report it as an effect size rather than as a verdict "
      "on each pair, because the resampled range is not calibrated for that use (Section 4.11). "
      "The size of the difference follows the split-half reliability of the feline state: the 24 "
      "pairs built on cortical CKD3/4 or medullary CKD1/2 (0.841, 0.851) have a median of 0.487, "
      "the 24 built on cortical CKD1/2 or medullary CKD3/4 (0.568, 0.616) a median of 0.338, and "
      "the same-species pairs 0.196. Those 24 pairs share two feline states between them, so they "
      "are not 24 independent replications.")

rp(D, "Under animal-level resampling the distance between an observed value and its benchmark "
      "could be resolved only where the feline state is reliably measured: 22 of the 24 "
      "cross-species pairs built on the two less reliable feline states have an interval "
      "containing zero, and we make no claim for those pairs. Three things set that range "
      "together — the reliability of the feline state, three disease samples in most mouse "
      "states, and the conditioning of every interval on the replicates in which a split-half "
      "estimate exists (Section 3.2) — and we do not attribute it to any one of them.",
      "Under animal-level resampling the distance between an observed value and its benchmark is "
      "interpretable only where the feline state is reliably measured: the pairs built on the two "
      "less reliable feline states carry a smaller median difference and a spread that reaches "
      "below zero, and we make no claim for them. We do not report how many pairs were resolved, "
      "because the resampled range covers the quantity it estimates in about half of simulated "
      "datasets and its lower limit exceeds zero even when the true difference is zero (Section "
      "4.11). Three things set that spread together — the reliability of the feline state, three "
      "disease samples in most mouse states, and the conditioning of every range on the "
      "replicates in which a split-half estimate exists (Section 3.2) — and we do not attribute "
      "it to any one of them.")

rp(M, "The benchmark assumes the two states’ errors are independent. States sharing a cohort or "
      "control samples violate this, and the excess of the observed value over the benchmark is "
      "reported as a lower bound on the resulting systematic error.",
      "The benchmark assumes the two states’ errors are independent, and states sharing a cohort "
      "or control samples violate that assumption. We do not convert an excess of the observed "
      "value over the benchmark into a bound on the resulting error: the simulation above shows "
      "that two states with identical true responses exceed the corrected benchmark in 78.7% of "
      "replicates with independent errors, so an excess is consistent with near-identical "
      "responses as well as with dependent ones. Which pairs exceed their benchmark is reported "
      "as an observation in Section 2.2.")

rp(M, "Separation was scored as the area under the curve, AUC, of within-species against "
      "cross-species pair values. Because the 120 pairs are generated by 16 states and are not "
      "independent, p values were obtained by permuting the species labels attached to the "
      "states, enumerating all 164=1,820 assignments exactly. A stricter level, permuting which "
      "dataset is designated feline, has a floor of p = 1/3 with three datasets; this is reported "
      "in the Discussion as a limitation rather than as a result.",
      "Separation was scored as the area under the curve, AUC, of within-species against "
      "cross-species pair values, and is reported as a description of the panel. The unit "
      "classified is the pair of states, not the state. No p value is attached to it. Because the "
      "120 pairs are generated by 16 states and are not independent, the only permutation "
      "available relabels the states, enumerating all C(16,4) = 1,820 assignments of the feline "
      "label; that null contains assignments in which a single dataset is split between species, "
      "which the design cannot produce, so it is not the null of a hypothesis this design can "
      "test. The stricter permutation, over which dataset is designated feline, has a floor of "
      "p = 1/3 with three datasets. Both are tabulated in Table S4 and neither is used for "
      "inference anywhere in the article.")
print('C: §3.2・§3.8・§4.10・§4.12 を §2.2 の方針にそろえた')

# =====================================================================
# D. 中心化効果の個体再標本化は §2.3 が初出。§3.3 はそれを参照する形に整理し、
#    Methods に回数・パーセンタイル・保持した構造の定義を書く。
# =====================================================================
rp(D, "Under an animal-level resampling that preserves the design, the uncentred value exceeds "
      "the centred one in 98.6% of replicates, so the direction is not an artefact of which cats "
      "and mice were sampled; the size of the gap is less certain, the drop running from 0.037 to "
      "0.438 across resamples (Section 2.3).",
      "Section 2.3 reports what an animal-level resampling does to this comparison: the direction "
      "reproduces in 98.6% of replicates while the size of the gap does not hold as tightly, so "
      "the effect is not an artefact of which cats and mice were sampled but its magnitude is "
      "estimated loosely.")

rp(M, "What is reported. The quantity is the paired difference",
      "Resampling the centring comparison. The same machinery is used once more in Section 2.3, "
      "on the two aggregated AUCs rather than on a pair of states. Animals are drawn with "
      "replacement 1,000 times under the pooling described above, so that states sharing a "
      "control group receive the same resampled controls and the feline cortical and medullary "
      "samples, being the same cats, are drawn together; all sixteen Δ vectors are rebuilt from "
      "each draw; the 422 pathway means are formed with and without removing each state's mean "
      "change across genes; and the AUC separating the 39 within-species from the 48 "
      "cross-species pairs that share no controls is computed for each. Reported are the "
      "proportion of replicates in which the uncentred AUC exceeds the centred one, and the "
      "median and the 2.5th and 97.5th percentiles of each AUC and of their difference over the "
      "1,000 replicates (`results/round7/centering_auc_animal_bootstrap.tsv`). No replicate was "
      "undefined here, because no split-half estimate is required.\n\n"
      "What is reported. The quantity is the paired difference")

# =====================================================================
# E. 仕上げ。§2.1 の 2 つの数式を LaTeX に戻し、文の直後に置く。
# =====================================================================
rp(R, "α=a⋅va⋅a,cosθ=a⋅v‖a‖‖v‖,ν=‖v‖‖a‖,R⊥=‖v-αa‖‖v‖",
      r"$$\alpha = \frac{a \cdot v}{a \cdot a}, \qquad \cos\theta = \frac{a \cdot v}"
      r"{\|a\|\,\|v\|}, \qquad \nu = \frac{\|v\|}{\|a\|}, \qquad "
      r"R_\perp = \frac{\|v - \alpha a\|}{\|v\|}$$")
rp(R, "α=νcosθ,  R⊥=1-cos2θ=sinθ.  1",
      r"$$\alpha = \nu \cos\theta, \qquad R_\perp = \sqrt{1 - \cos^2\theta} = \sin\theta "
      r"\tag{1}$$")
print('D と E-1: 再標本化の定義と §2.1 の数式を整えた')

# E の続き。Table 1 のカウント列も 1 セル 1 値にする（Table 4 と同じ原則）。
rp(R, "| Pair class | n pairs | cos θ median | cos θ max | Benchmark, Spearman–Brown: median "
      "[min] | Benchmark, uncorrected: median [min] | Pairs above benchmark (S–B / uncorr.) |",
      "| Pair class | n pairs | cos θ median | cos θ max | Benchmark, Spearman–Brown: median "
      "[min] | Benchmark, uncorrected: median [min] | Pairs above the Spearman–Brown benchmark | "
      "Pairs above the uncorrected benchmark |")
rp(R, "| Within dataset | 45 | 0.755 | 0.958 | 0.920 | 0.852 | 4 / 8 |",
      "| Within dataset | 45 | 0.755 | 0.958 | 0.920 | 0.852 | 4 | 8 |")
rp(R, "| Same species, different dataset | 27 | 0.707 | 0.884 | 0.945 | 0.896 [0.781] | 0 / 1 |",
      "| Same species, different dataset | 27 | 0.707 | 0.884 | 0.945 | 0.896 [0.781] | 0 | 1 |")
rp(R, "| Cross-species (cat vs mouse) | 48 | 0.448 | 0.548 | 0.862 [0.765] | 0.765 [0.620] | 0 / 0 |",
      "| Cross-species (cat vs mouse) | 48 | 0.448 | 0.548 | 0.862 [0.765] | 0.765 [0.620] | 0 | 0 |")
print('E: Table 1 のカウント列を 2 列に分けた')
