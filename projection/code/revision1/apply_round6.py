"""第6ラウンドの本文反映。

Part 1-2 / 1-4: 条件 (iv)(v) で Spearman-Brown 補正済みベンチマークが較正されている
（真の cos < 1 の 120 条件すべてで超過ゼロ、同一応答で 78.7% 超過）ことが分かったので、
主たる基準を SB 側に切り替える。無補正形は保守的な参照として併記する。
Part 3-1: 補足要素のキャプションと参照の対応を整える。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'
C, F, S = 'manuscript/CONCLUSIONS.md', 'manuscript/FRONTMATTER.md', 'manuscript/SUPPLEMENTARY.md'
BM = 'manuscript/BACKMATTER.md'

# ---- §2.2: 主たる基準を SB に切り替える ----
rp(R, "We report it both uncorrected and with the Spearman–Brown correction, and we make the "
      "primary claim on the distance between the observed value and the benchmark rather than on "
      "any ratio of the two (Figure 2B,C; corrected values are given in Table S1 and are not used "
      "in the text).",
      "We report it both with the Spearman–Brown correction and uncorrected, and make the primary "
      "claim on the corrected form, because simulation of the whole procedure places that form at "
      "the value two states with identical true responses attain while the uncorrected form sits "
      "below it (Section 4.10, Figure S1): the uncorrected benchmark is estimated from half "
      "samples and compared with a full-sample observation, and the correction is the "
      "extrapolation that reconciles the two. The claim is made on the distance between the "
      "observed value and the benchmark rather than on any ratio of the two (Figure 2B,C; the "
      "attenuation-corrected cos θ values are given in Table S1 and are not used in the text).")

rp(R, "No cross-species pair approaches its benchmark. Under the uncorrected — that is, the most "
      "conservative — benchmark, the lowest benchmark among the 48 cross-species pairs is 0.620, "
      "above the highest observed cross-species value of 0.548; the closest any cross-species "
      "pair comes to its own benchmark is a gap of 0.171. None of the 48 exceeds either benchmark "
      "(Table 1).",
      "No cross-species pair approaches its benchmark. The lowest Spearman–Brown benchmark among "
      "the 48 cross-species pairs is 0.765, above the highest observed cross-species value of "
      "0.548, and the closest any cross-species pair comes to its own benchmark is a gap of "
      "0.297. Under the uncorrected form, which sits lower, the corresponding figures are 0.620 "
      "and 0.171. None of the 48 exceeds either benchmark (Table 1).")
print('§2.2: 基準の位置づけを SB 主に切り替え')

# ---- Part 1-4: 個体ブートストラップの本数と中央値を SB 基準に置き換える ----
rp(R, "The difference stays above zero in 26 of the 48 cross-species pairs, and which pairs those "
      "are follows the split-half reliability of the feline state in the pair. The 24 pairs whose "
      "feline state is cortical CKD3/4 or medullary CKD1/2 (reliabilities 0.841 and 0.851) all "
      "keep the difference above zero, with a median of 0.403 and lower bounds running from 0.003 "
      "to 0.347. The largest tail proportion among them is 0.026 — the plain share of defined "
      "replicates at or below zero (Section 4.11), not a test p value — and it exceeds 0.025 for "
      "the pair whose lower bound is 0.003 because the percentile interpolates between order "
      "statistics that straddle zero. Of the 24 whose feline state is cortical CKD1/2 (0.568) or "
      "medullary CKD3/4 (0.616, from five animals), 22 have an interval containing zero: the "
      "difference could not be resolved there, which is not the same as showing those pairs reach "
      "their benchmark. A leave-one-animal-out jackknife, a much milder perturbation, puts the "
      "difference below zero in 1 of the 48. Among the 39 same-species pairs sharing no controls, "
      "24 keep the difference above zero and 15 could not be resolved.",
      "Against the Spearman–Brown benchmark the difference stays above zero in 32 of the 48 "
      "cross-species pairs, and which pairs those are follows the split-half reliability of the "
      "feline state in the pair. The 24 pairs whose feline state is cortical CKD3/4 or medullary "
      "CKD1/2 (reliabilities 0.841 and 0.851) all keep the difference above zero, with a median "
      "of 0.487 and lower bounds running from 0.142 to 0.501. The largest tail proportion among "
      "them is 0.020 — the plain share of defined replicates at or below zero (Section 4.11), not "
      "a test p value. Of the 24 whose feline state is cortical CKD1/2 (0.568) or medullary "
      "CKD3/4 (0.616, from five animals), 16 have a range containing zero: the difference could "
      "not be resolved there, which is not the same as showing those pairs reach their benchmark. "
      "Among the 39 same-species pairs sharing no controls, 34 keep the difference above zero and "
      "5 could not be resolved. Under the uncorrected benchmark the same counts are 26 of 48, 24 "
      "of 24 with a median of 0.403 and a smallest lower bound of 0.003, 2 of 24 and 24 of 39; "
      "the ordering of the groups is the same and only the level moves. A leave-one-animal-out "
      "jackknife on the uncorrected form, a much milder perturbation, puts the difference below "
      "zero in 1 of the 48.")

rp(R, "Group C, the 24 pairs built on the two less reliably measured feline states, has 22 "
      "intervals containing zero and a median difference of 0.215 against 0.119 for same-species "
      "pairs; we do not offer it as evidence in either direction, and no claim below rests on it.",
      "Group C, the 24 pairs built on the two less reliably measured feline states, has 16 ranges "
      "containing zero and a median difference of 0.338 against 0.196 for same-species pairs; we "
      "do not offer it as evidence in either direction, and no claim below rests on it.")

rp(R, "Taken as a yes-or-no question, falling short of the benchmark is not a cross-species "
      "phenomenon at all: the interval excludes zero for 24 of the 39 same-species pairs (62%) "
      "and for 26 of the 48 cat-mouse pairs (54%). The distributions differ in magnitude. Under "
      "the uncorrected benchmark the median difference is 0.119 (interquartile range 0.065-0.291) "
      "for same-species pairs, 0.403 (0.350-0.469) for cat-mouse pairs with a reliably measured "
      "feline state and 0.215 (0.189-0.299) for the rest; under the Spearman-Brown benchmark the "
      "three are 0.196 (0.117-0.360), 0.487 (0.437-0.561) and 0.338 (0.322-0.434).",
      "Taken as a yes-or-no question, falling short of the benchmark is not a cross-species "
      "phenomenon at all: the range excludes zero for 34 of the 39 same-species pairs (87%) and "
      "for 32 of the 48 cat-mouse pairs (67%). The distributions differ in magnitude. Under the "
      "Spearman-Brown benchmark the median difference is 0.196 (interquartile range 0.117-0.360) "
      "for same-species pairs, 0.487 (0.437-0.561) for cat-mouse pairs with a reliably measured "
      "feline state and 0.338 (0.323-0.434) for the rest; under the uncorrected benchmark the "
      "three are 0.119 (0.065-0.291), 0.403 (0.350-0.469) and 0.215 (0.189-0.299).")

rp(R, "The observed gap between the two class medians is +0.219 under the uncorrected benchmark "
      "and +0.239 under the Spearman-Brown one.",
      "The observed gap between the two class medians is +0.239 under the Spearman-Brown "
      "benchmark and +0.219 under the uncorrected one.")
print('Part 1-4: 個体ブートストラップの数値を SB 主に')

# 第5ラウンドで足した「SB に直すとこう動く」という補足段落は、SB が主になった以上
# 重複する。較正の根拠だけを残し、数の繰り返しを削る。
rp(R, "The Spearman-Brown form is the calibrated one — in simulation it is never exceeded by a "
      "pair whose true responses differ, while the uncorrected form is low because it compares a "
      "half-sample estimate with a full-sample observation (Section 4.10) — so every figure given "
      "here on the uncorrected benchmark understates the shortfall against the calibrated one. "
      "Repeating the resampling against the corrected benchmark leaves the group-B statement "
      "unchanged and less dependent on its smallest value: all 24 lower limits stay above zero, "
      "and the smallest is 0.142 rather than 0.003. It also moves the other two groups, to 8 of "
      "24 in group C and 34 of the 39 same-species pairs, and to 32 of the 48 cat-mouse pairs "
      "overall, which does not change the reading that falling short is not a cross-species "
      "phenomenon (`results/revision1/pair_uncertainty_both_ceilings.tsv`).",
      "The two forms order the three groups the same way and differ only in level, so the choice "
      "between them does not bear on the comparison; it bears on where zero sits "
      "(`results/revision1/pair_uncertainty_both_ceilings.tsv`).")
print('§2.2 の重複段落を整理')

# ---- Table 1: 補正済みの列を先に置く ----
rp(R, "| Pair class | n pairs | cos θ median | cos θ max | Benchmark, uncorrected: median [min] | "
      "Benchmark, Spearman–Brown: median [min] | Pairs above benchmark (uncorr. / S–B) |",
      "| Pair class | n pairs | cos θ median | cos θ max | Benchmark, Spearman–Brown: median "
      "[min] | Benchmark, uncorrected: median [min] | Pairs above benchmark (S–B / uncorr.) |")
rp(R, "| Within dataset | 45 | 0.755 | 0.958 | 0.852 | 0.920 | 8 / 4 |",
      "| Within dataset | 45 | 0.755 | 0.958 | 0.920 | 0.852 | 4 / 8 |")
rp(R, "| Same species, different dataset | 27 | 0.707 | 0.884 | 0.896 [0.781] | 0.945 | 1 / 0 |",
      "| Same species, different dataset | 27 | 0.707 | 0.884 | 0.945 | 0.896 [0.781] | 0 / 1 |")
rp(R, "| Cross-species (cat vs mouse) | 48 | 0.448 | 0.548 | 0.765 [0.620] | 0.862 [0.765] | 0 / 0 |",
      "| Cross-species (cat vs mouse) | 48 | 0.448 | 0.548 | 0.862 [0.765] | 0.765 [0.620] | 0 / 0 |")
rp(R, "**Table 1.** *Observed cos θ by pair class, with reliability-derived benchmarks.*",
      "**Table 1.** *Observed cos θ by pair class, with reliability-derived benchmarks. The "
      "Spearman–Brown form is given first because simulation places it at the value two states "
      "with identical true responses attain, while the uncorrected form sits below that (Section "
      "4.10).*")

# ---- Table 4 と §2.6: 基準の行も SB に ----
rp(R, "| Lowest uncorrected benchmark, cat–mouse | 0.620 | 0.641 | 0.614 | 0.627 |",
      "| Lowest Spearman–Brown benchmark, cat–mouse | 0.765 | 0.780 | 0.760 | 0.770 |")
rp(R, "No cat–mouse pair reaches its benchmark in any space — the highest observed value stays "
      "below the lowest benchmark (0.548 vs 0.620; 0.572 vs 0.641; 0.453 vs 0.614; 0.530 vs "
      "0.627) and none of the 48 pairs exceeds it.",
      "No cat–mouse pair reaches its benchmark in any space — the highest observed value stays "
      "below the lowest Spearman–Brown benchmark (0.548 vs 0.765; 0.572 vs 0.780; 0.453 vs 0.760; "
      "0.530 vs 0.770) and none of the 48 pairs exceeds it, under either form of the benchmark.")

# ---- §3.2: 群間差を SB 主に ----
rp(D, "it happens about as often within species as across them - the interval excludes zero for "
      "62% of same-species pairs and 54% of cat-mouse pairs - and what separates the two classes "
      "is how large the shortfall is. The median difference is 0.119 within species against 0.403 "
      "and 0.215 in the two cross-species groups, and permuting which four states are feline "
      "places the observed difference in medians, +0.219, against a null median of -0.012 (exact "
      "p = 0.0011)",
      "it happens about as often within species as across them - the range excludes zero for 87% "
      "of same-species pairs and 67% of cat-mouse pairs - and what separates the two classes is "
      "how large the shortfall is. The median difference is 0.196 within species against 0.487 "
      "and 0.338 in the two cross-species groups, and permuting which four states are feline "
      "places the observed difference in medians, +0.239, against a null median of -0.012 (exact "
      "p = 0.0005)")
print('Table 1 / Table 4 / §3.2 を SB 主に')

# ---- Abstract と Conclusions を SB 基準の数値に ----
rp(F, "which can be resolved only for the 24 pairs with a reliably measured feline state (median "
      "0.403 against 0.119 within species).",
      "which can be resolved only for the 24 pairs with a reliably measured feline state (median "
      "0.487 against 0.196 within species).")
rp(C, "under animal-level resampling the shortfall was resolved for the 24 cat–mouse pairs built "
      "on the two reliably measured feline states and not for 22 of the others.",
      "under animal-level resampling the shortfall was resolved for the 24 cat–mouse pairs built "
      "on the two reliably measured feline states and not for 16 of the others.")
print('Abstract と Conclusions を SB 基準に')

# ---- Part 2 / Part 3-1: Table S1 の出典を、ペアごとに両方の基準を持つファイルに差し替える ----
rp(S, "**Table S1.** Attenuation-corrected cos θ and the observed-to-benchmark ratio for all 120 "
      "state pairs. These values are not used in the main text "
      "(results/reliability/pairs_corrected.tsv).",
      "**Table S1.** Per-pair benchmarks and attenuation-corrected cos θ for all 120 state pairs: "
      "the observed cos θ, the Spearman–Brown and uncorrected benchmarks, the corrected cos θ and "
      "its interval, and the pair's class and shared-control status. The corrected cos θ and the "
      "observed-to-benchmark ratio are not used in the main text "
      "(results/reliability/pairs_ceilings.tsv).")
rp(BM, "Table S1, the attenuation-corrected cos θ and the observed-to-benchmark ratio for all 120 "
       "state pairs, which are not used in the main text "
       "(results/reliability/pairs_corrected.tsv).",
       "Table S1, the per-pair benchmarks in both forms together with the attenuation-corrected "
       "cos θ and the observed-to-benchmark ratio for all 120 state pairs, the last two of which "
       "are not used in the main text (results/reliability/pairs_ceilings.tsv).")
rp(BM, "Figure S1, the simulation separating the two reasons the reliability-derived benchmark is "
       "not an upper bound (results/round5/ceiling_decomposition_sim.tsv).",
       "Figure S1, the simulation showing which of the two benchmark forms is an upper bound and "
       "why (results/round5/ceiling_decomposition_sim.tsv).")
print('Part 2 / 3-1: Table S1 の出典と背表紙の一覧を更新')
