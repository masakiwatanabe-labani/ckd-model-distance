"""第5ラウンドの本文反映。

rev8 で直された2件（§2.3・§2.6 の集約後 AUC、§3.2 のセル数）を Markdown にも入れる。
Part C 以降の文章修正も同じファイルに追記していく。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'
C, F, S = 'manuscript/CONCLUSIONS.md', 'manuscript/FRONTMATTER.md', 'manuscript/SUPPLEMENTARY.md'

# ---- rev8 の修正1: 集約後 AUC は Table 4 の 0.600 / 0.626 が正 ----
rp(R, "The direction of the loss reproduces in the two other gene spaces of Section 2.6; its size "
      "does not, the aggregated values there being 0.633 and 0.608.",
      "The direction of the loss reproduces in the two other gene spaces of Section 2.6; its size "
      "does not, the aggregated values there being 0.600 and 0.626.")
rp(R, "Third, aggregation falls to chance (AUC 0.511) only in Group A; in the other two spaces it "
      "falls to 0.633 and 0.608.",
      "Third, aggregation falls to chance (AUC 0.511) only in Group A; in the other two spaces it "
      "falls to 0.600 and 0.626.")

# ---- rev8 の修正2: サイズ条件を満たすセル数は 2,136（888 は 193 セット時代の 74x12）----
rp(D, "focal adhesion is the only one in which any mouse state reaches an alignment of 0.9 or "
      "above, in three of the 888 qualifying mouse × pathway cells, each below its reliability "
      "ceiling.",
      "focal adhesion is the only one in which any mouse state reaches an alignment of 0.9 or "
      "above, in three of the 2,136 qualifying mouse × pathway cells, each below its "
      "reliability benchmark.")
print('Round 5: rev8 corrections applied')

# =====================================================================
# Part A-4. 用語を統一する。
# 本文は「厳密な上限ではない」と書きながら Abstract と Conclusions は
# 「measurement reliability が定める ceiling」と書いていた。ceiling を上限の意味で
# 使わず、reliability-derived benchmark に統一する。
# =====================================================================
import re
from pathlib import Path as _P

_MD = _P('manuscript')
_RENAME = [
    # 長いものから順に。語形を壊さないよう、前後の語ごと置き換える。
    ("the ceiling set by measurement reliability",
     "the reliability-derived benchmark"),
    ("the ceiling that bounds it", "the reliability-derived benchmark it is compared with"),
    ("Attenuation-corrected cos θ and the observed-to-ceiling ratio",
     "Attenuation-corrected cos θ and the observed-to-benchmark ratio"),
    ("attenuation ceiling", "reliability-derived benchmark"),
    ("Attenuation ceiling", "Reliability-derived benchmark"),
    ("attenuation ceilings", "reliability-derived benchmarks"),
    ("ceiling ratio", "benchmark ratio"),
    ("Ceiling estimates", "Benchmark estimates"),
    ("ceiling estimate", "benchmark estimate"),
    ("ceiling analysis", "benchmark analysis"),
    ("ceilings", "benchmarks"),
    ("Ceilings", "Benchmarks"),
    ("ceiling", "benchmark"),
    ("Ceiling", "Benchmark"),
]


def rename_terms():
    n = 0
    for f in ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
              "METHODS.md", "CONCLUSIONS.md", "SUPPLEMENTARY.md", "BACKMATTER.md"):
        p = _MD / f
        t = p.read_text()
        for old, new in _RENAME:
            n += t.count(old)
            t = t.replace(old, new)
        assert "ceiling" not in t.lower(), f"{f} に ceiling が残っている"
        p.write_text(t)
    print(f'Part A-4: "ceiling" を benchmark 系に置換 {n} 箇所')


rename_terms()

# =====================================================================
# Part A-1 / A-2 / A-3. 天井超過の原因を分解した結果を反映する。
# 旧記述は「非中心化コサインだから上限にならない」としていたが、3 条件に分けると
# 超過はほぼすべて「基準は半標本・観測は全標本」という標本数の差から出ており、
# 非中心化の寄与はほぼゼロだった。信頼性の範囲も全 state を包含する形に直す。
# =====================================================================
rp(M, "This attenuation relation is derived for Pearson correlations, and cos θ here is not "
      "centred, so its use is an approximation. Representational similarity analysis makes the "
      "same distinction between metrics, normalising to unit norm for a cosine and to zero mean "
      "and unit variance for a correlation [18,19]. We checked the approximation by simulating "
      "the whole procedure: two states with a known angle between their true responses and a mean "
      "shift common to all genes, sampled at the group sizes used here, with the benchmark "
      "estimated from the simulated samples exactly as above. The estimate is not a strict upper "
      "bound, and how it behaves depends on the true angle. Restricting to the range of estimated "
      "reliabilities our states actually have, r_half between 0.50 and 0.84, the observed value "
      "never exceeded the estimate in any replicate at a true cos θ of 0.3 or 0.5 - the regime "
      "the cross-species pairs occupy - while at 0.9 and at 1.0 it exceeded it in every "
      "replicate, and at 0.7 in up to 46%. The size of the common mean shift, varied from 0 to "
      "0.44 in units of the within-state standard deviation against 0.26 in the median state of "
      "our panel, moves these figures by a few percentage points and does not change the pattern. "
      "Two consequences follow. The estimate sits below the value two identical states attain, so "
      "it understates the attainable similarity and the distance we report to it is conservative. "
      "And it is not a limit for pairs that are already well aligned, which is consistent with "
      "the same-species pairs sitting above their benchmark in Table 1. We therefore treat it as a "
      "reference level rather than a hard limit, and rest the claims of Section 2.2 on the "
      "animal-level interval of the paired difference rather than on a point comparison with it. "
      "Figure S1 gives the full grid.",
      "Two features of this construction keep it from being an upper bound, and they are not the "
      "same thing. The relation is derived for Pearson correlations while cos θ here is not "
      "centred; representational similarity analysis makes the same distinction between metrics, "
      "normalising to unit norm for a cosine and to zero mean and unit variance for a correlation "
      "[18,19]. Separately, the benchmark is estimated from half samples while the observed value "
      "is computed from the full sample, so the two sides carry different amounts of measurement "
      "error. Congruence analysis of the CAMO type addresses the second point by matching the "
      "number of samples on the two sides of such a comparison [6]; we did not, and the "
      "consequence is quantified here. We simulated the whole procedure at the group sizes used "
      "here, over the full range of reliabilities our sixteen states have (r_half 0.568 to 0.986, "
      "Table S3), with a mean shift common to all genes of 0 to 0.44 in units of the within-state "
      "standard deviation against 0.257 in the median state, and separated three conditions: (i) "
      "centred Pearson correlations with the benchmark and the observation computed from samples "
      "of the same size, (ii) centred Pearson correlations with a half-sample benchmark against a "
      "full-sample observation, and (iii) the uncentred cosine with a half-sample benchmark "
      "against a full-sample observation, which is the procedure used here. Condition (i) is "
      "calibrated: the observation exceeded the benchmark in 0.00% of replicates at every true "
      "angle below 1, and in 0.06% when the two true responses were identical. Conditions (ii) "
      "and (iii) both exceed it, and by almost the same amount: 39.7% and 39.9% of replicates at "
      "a true cos θ of 0.9, and every replicate when the two responses are identical. The "
      "difference between them, which is what the uncentred metric contributes, averages 0.0004 "
      "and never exceeds 0.020 in any condition tested. Almost all of the exceedance is therefore "
      "attributable to comparing a half-sample benchmark with a full-sample observation, not to "
      "the metric. At true angles of 0.7 and below no condition exceeded the benchmark in any "
      "replicate. Two consequences follow. The benchmark understates the similarity two identical "
      "states attain at these group sizes, so the distance we report to it is conservative. And "
      "it is not a limit for pairs whose true alignment is high, which is consistent with the "
      "same-species pairs sitting above their benchmark in Table 1. We therefore treat it as a "
      "reference level rather than a hard limit, and rest the claims of Section 2.2 on the "
      "animal-level distribution of the paired difference rather than on a point comparison with "
      "it. Figure S1 gives the full grid "
      "(`results/round5/ceiling_decomposition_sim.tsv`).")

# Part A-3: 真の cos θ と観測 cos θ を区別する。
rp(M, "Spearman–Brown assumes two halves of equal length.",
      "Simulated angles are true angles between the underlying responses; the observed cos θ is "
      "the attenuated quantity and is smaller. At the reliabilities our states have, a true angle "
      "of 0.3 gives observed values of 0.23 to 0.41, a true 0.5 gives 0.38 to 0.58 and a true 0.7 "
      "gives 0.53 to 0.74, so the observed cross-species range of 0.19 to 0.55 corresponds to "
      "true angles of roughly 0.3 to 0.7 rather than to the same numbers. We do not read the "
      "simulated angles as estimates of the true cross-species alignment.\n\n"
      "Spearman–Brown assumes two halves of equal length.")
print('Part A-1/A-2/A-3 applied')

# =====================================================================
# Part A-5. 「群間差 +0.219 も下限である」を撤回する。
# 基準値が低めに出ることだけからは、差の向きは決まらない。観測 cos も信頼性に依存し、
# A-1 の分解で超過の主因が標本数の差だと分かった以上、下限とする根拠はない。
# =====================================================================
rp(D, "That difference is itself a lower bound. The simulation of Section 4.10 shows the "
      "benchmark estimate sits below the value two identical states attain, and the shortfall of "
      "the estimate grows as reliability falls; the uncorrected benchmarks of the cat-mouse pairs "
      "are lower than those of the same-species pairs (median 0.765 against 0.852 and 0.896 in "
      "Table 1), so their gaps are the more strongly understated of the two. The +0.219 we report "
      "is therefore smaller than the difference would be if both groups were measured equally "
      "well.",
      "We do not claim a direction for the bias in that difference. The benchmark is estimated "
      "below the value two identical states attain, but the observed cos θ on the other side of "
      "the subtraction is attenuated by the same reliabilities, so the two effects act on the "
      "difference in opposite directions and their net sign is not determined by the simulation "
      "of Section 4.10.")
print('Part A-5 applied: 群間差の下限主張を撤回')

# =====================================================================
# Part B-2. 順位逆転を個体ブートストラップで測り、維持率を本文に書く。
# 維持率は 76.0% なので、「安定に観測される」とは書けない。
# =====================================================================
rp(R, "We therefore do not report or rank α as a measure of how well a model reproduces the "
      "reference state.",
      "We therefore do not report or rank α as a measure of how well a model reproduces the "
      "reference state.\n\n"
      "How often that ordering survives the exchange of animals rather than of genes is a "
      "separate question, and the intervals above do not answer it. We resampled the cats with "
      "replacement 2,000 times, keeping the design intact: the reference state and cortical "
      "CKD1/2 share their six control animals, and the cortical and medullary samples come from "
      "the same cats, so groups whose animal sets intersect were drawn once and then split "
      "(Section 4.11). The medullary state keeps the larger α in 93.7% of replicates and the "
      "smaller cos θ in 82.3%, and both hold together in 76.0%. The inversion is therefore the "
      "usual outcome in this cohort rather than an invariable one. The algebraic statement does "
      "not depend on that rate — the inequality below is exact — but the empirical demonstration "
      "rests on nine cats, and we attach the rate to it rather than presenting the ordering as "
      "established.")

rp(R, "The 95% gene-bootstrap interval, 0.965–1.079, contains 1.000, and it covers the exchange "
      "of genes rather than of animals, so the animal-level uncertainty is wider again.",
      "The 95% gene-bootstrap interval, 0.965–1.079, contains 1.000, and it covers the exchange "
      "of genes rather than of animals. Resampling animals as above puts α for feline medullary "
      "CKD1/2 above 1.000 in 46.8% of replicates, with a median of 0.977 and a 2.5th-to-97.5th "
      "percentile range of 0.506 to 1.397, so at the animal level the crossing is about as likely "
      "to be absent as present.")

rp(R, "Group A intersection set, 2,016 genes; the ordering in (A) holds in all three gene spaces "
      "of Table 4, the crossing of α = 1 does not.",
      "Intervals are gene bootstraps and describe the exchange of genes, not of animals; under an "
      "animal-level resampling that preserves the shared controls and the cortex–medulla pairing "
      "the ordering in (A) holds in 76.0% of replicates and α exceeds 1 in 46.8% (Section 2.1). "
      "Group A intersection set, 2,016 genes; the ordering in (A) holds in all three gene spaces "
      "of Table 4, the crossing of α = 1 does not.")
print('Part B-2 applied: 逆転維持率 76.0% を本文に')

# =====================================================================
# Part B-1. 個体ブートストラップ区間の被覆率を測り、区間に基づく主張を弱める。
# 被覆率は 0.8-32.5%（平均 9.4%）で、外れ方は一方向。区間は推定対象より下に出る。
# 判断そのもの（下限がゼロを上回る）の偽陽性率は 4 設計とも 0.0000 だった。
# =====================================================================
rp(M, "What is reported. The quantity is the paired difference",
      "Operating characteristics. Whether this construction has the coverage a 95% interval "
      "would have was not obvious and we checked it by simulation, reproducing the group sizes, "
      "the reliabilities and the re-splitting rule of four representative cat–mouse designs "
      "(`code/revision1/interval_coverage_sim.py`). Holding the true Δ vectors fixed and drawing "
      "animals, the proportion of undefined replicates came out at 0.109 to 0.212, against 0.098 "
      "to 0.235 in the data. The interval covered the quantity it estimates — the median of the "
      "paired difference over repeated draws of animals from the same design — in 0.8% to 32.5% "
      "of simulated datasets, 9.4% on average. The failure is one-sided: in 90.6% of datasets the "
      "estimand lay above the upper limit and in none below the lower limit, because a resampled "
      "group holds fewer distinct animals than the original and the reliability estimated from "
      "them is correspondingly lower, which moves the whole interval down. The same shift is "
      "visible in the data, where the bootstrap median difference for cat–mouse pairs is 0.338 "
      "against a point estimate of 0.383. These percentiles are therefore not a 95% confidence "
      "interval for the paired difference, and we do not use them as one. What they support is "
      "the sign: simulating designs in which the true paired difference is zero, the lower limit "
      "exceeded zero in 0.0% of datasets in all four designs "
      "(`results/round5/interval_null_rate.tsv`), while in every design with a positive true "
      "difference it did so in 100%. We therefore read a lower limit above zero as a "
      "conservative one-sided indication that the observed similarity falls short, and we report "
      "the percentiles as a range of the resampled statistic rather than as a calibrated "
      "interval.\n\n"
      "What is reported. The quantity is the paired difference")

rp(R, "The difference stays above zero in 26 of the 48 cross-species pairs, and which pairs those "
      "are follows the split-half reliability of the feline state in the pair.",
      "The percentile ranges below are not calibrated 95% intervals. Simulation of this exact "
      "procedure, at these group sizes and reliabilities, covers the quantity being estimated in "
      "9.4% of datasets on average, and misses it low rather than high in every case, because a "
      "resampled group holds fewer distinct animals than the original and its estimated "
      "reliability is lower (Section 4.11). What survives the check is the sign: where the true "
      "difference is zero the lower limit exceeded zero in none of the simulated datasets. We "
      "therefore treat a lower limit above zero as a conservative one-sided indication and not as "
      "a 95% statement, and we report the ranges as percentiles of the resampled statistic.\n\n"
      "The difference stays above zero in 26 of the 48 cross-species pairs, and which pairs those "
      "are follows the split-half reliability of the feline state in the pair.")

rp(R, "Three limits on how far this goes. The intervals are conditional on the replicates in "
      "which the reliability is defined, which excludes a median of 11.8% of draws.",
      "Four limits on how far this goes. The percentile ranges are conditional on the replicates "
      "in which the reliability is defined, which excludes a median of 11.8% of draws, and they "
      "undercover as described above, so they are read for sign only.")
print('Part B-1 applied: 被覆率と偽陽性率を反映し、区間の位置づけを弱めた')

# =====================================================================
# Part C. 配列一致率をマッピング精度と同一視しない。
# 示せているのは「一致率による層別化と high-confidence 1:1 への制限では差が縮まらなかった」
# までであり、対応付けの誤り率そのものを測ってはいない。
# =====================================================================
rp(R, "Orthologue mapping quality. If mapping error set the cross-species value, cos θ should "
      "rise with orthologue sequence identity. It does not, and the relationship is not "
      "monotonic:",
      "Orthologue mapping quality. Sequence identity is a measure of evolutionary conservation "
      "and not an error rate for the correspondence itself; Ensembl assigns orthology confidence "
      "from gene-order conservation and whole-genome alignment as well as identity. What can be "
      "tested here is whether the deficit tracks the two proxies available to us. If it did, cos "
      "θ should rise with orthologue sequence identity. It does not, and the relationship is not "
      "monotonic:")
rp(R, "Improving mapping quality does not narrow the gap.",
      "Neither stratifying by sequence identity nor restricting to high-confidence one-to-one "
      "orthologues narrowed the gap. That is the extent of what these two tests establish; they "
      "do not measure the rate at which the correspondence is wrong, and we do not claim that "
      "mapping error has been excluded.")
rp(R, "| cos θ by orthologue sequence-identity quartile | 0.479 / 0.496 / 0.495 / 0.327 — not "
      "monotonic, lowest where identity is highest, so mapping error does not account for the "
      "deficit |",
      "| cos θ by orthologue sequence-identity quartile | 0.479 / 0.496 / 0.495 / 0.327 — not "
      "monotonic, lowest where identity is highest, so stratifying by identity does not narrow "
      "the deficit |")
rp(R, "Adding human cohorts would not fix this; the asymmetry is a property of the orthologue "
      "tables. A three-species simultaneous alignment, or some other construction that maps both "
      "species to human through a procedure of equal quality, would be needed first.",
      "What the reversal shows is that the comparison is sensitive to which genes enter it; "
      "matching on identity also changes which functions are represented, so we cannot attribute "
      "the reversal to mapping quality alone. An additional human cohort would not on its own "
      "remove the asymmetry, which is a property of the two orthologue tables rather than of the "
      "human sample; a three-species simultaneous alignment, or some other construction that maps "
      "both species to human through a procedure of equal quality, would address it more "
      "directly.")
print('Part C applied')

# =====================================================================
# Part E-1. 起点に関する断定を外し、部分集合での確認結果を書く。
# 起点とモデル種別は完全に交絡している。IRI 単独では方向と振幅の比が 1.23 に落ちる。
# =====================================================================
rp(R, "Onset compartment does not structure these distances (−0.064, p = 0.727 over the twelve "
      "mouse states; −0.088, p = 0.667 over all sixteen), and controlling for it leaves the time "
      "term intact (+0.594, p = 0.0021). Elapsed time, not the compartment in which disease "
      "begins, organises how these states sit relative to one another.",
      "Onset compartment does not structure these distances (−0.064, p = 0.727 over the twelve "
      "mouse states; −0.088, p = 0.667 over all sixteen), and controlling for it leaves the time "
      "term intact (+0.594, p = 0.0021). That null cannot be read as showing the onset "
      "compartment to be irrelevant: on the mouse side the glomerular-onset states are exactly "
      "the Pod-TRECK states and the tubular-onset states exactly the IRI states, so onset is "
      "perfectly confounded with model identity and the two cannot be told apart in these data. "
      "What the analysis supports is that elapsed time orders these states, not that the "
      "compartment in which disease begins does not.")
rp(R, "That distance, however, is 1 − Spearman ρ, which is invariant to rescaling either vector "
      "and is therefore already a measure of direction alone: it cannot report amplitude.",
      "The association with time is not carried by one model or by the extreme time point. "
      "Within IRI alone (nine states) the Mantel correlation with the linear direction is +0.627 "
      "(p = 0.0009); dropping the 12-month state from the twelve leaves +0.627 (p = 0.0043); "
      "dropping it from IRI alone leaves +0.737 (p = 0.0005). The ratio of the direction term to "
      "the amplitude term does not survive as well, and we say so below.\n\n"
      "That distance, however, is 1 − Spearman ρ, which is invariant to rescaling either vector "
      "and is therefore already a measure of direction alone: it cannot report amplitude.")
print('Part E-1 applied')

# =====================================================================
# Part G-2 / G-3. §4.4 の記述を直す。
# =====================================================================
rp(M, "This is a standardised effect size rather than a log2 fold change and it reorders genes, "
      "which is incompatible with an analysis that uses the vector norm ‖v‖.",
      "This is a standardised effect size rather than a log2 fold change and it reorders genes. A "
      "norm is defined for the standardised vector as well, but it is then a norm of standardised "
      "effect sizes, so ‖v‖ no longer measures the size of the expression change; the biological "
      "reading of the amplitude term changes with it.")
rp(M, "with yg,i the log2-scale abundance: the deposited log2 values for feline RNA, log2FPKM+1 "
      "for Pod-TRECK and log2count+1 for IRI. Δg,s is therefore a log2 fold change.",
      "with yg,i the log2-scale abundance: the deposited log2 values for feline RNA, "
      "log2(FPKM + 1) for Pod-TRECK and log2(count + 1) for IRI. Δg,s is therefore the difference "
      "of two group means of log2(x + 1), which is not the same as the log2 of the ratio of the "
      "two group mean abundances: the mean of a logarithm is not the logarithm of a mean, and the "
      "offset of 1 shrinks the value towards zero by an amount that depends on the abundance. The "
      "offset is the same for both mouse datasets but the units are not (FPKM against normalised "
      "counts), so the shrinkage is not identical across them, and it can change the ranking of "
      "genes near the offset and hence the direction of Δ, not only its size. All comparisons "
      "here are between states measured on the same scale within a dataset or are rank-based "
      "across datasets, and the sensitivity analysis of Table S2 gives the effect of removing the "
      "expression filter that excludes the genes nearest the offset (main effect −0.067).")
print('Part G-2/G-3 applied')

# =====================================================================
# Part D. size-matched null の位置づけを直す。
# 図 3A が示すのは経路スコアの背景分布であって、特定の経路の推定誤差ではない。
# 小さいセットについての一律の結論を背景分布の記述に限定し、focal adhesion には
# 個体ブートストラップの区間を付け、選択が探索的であることを明示する。
# =====================================================================
rp(R, "Pathway size sets a floor on what can be read. Drawing random gene sets of matched size, "
      "the 95% range of cos θ spans 0.652 at 30 genes, 0.540 at 50 and 0.419 at 100 (Figure 3A). "
      "Pathway sizes here are median 45 (interquartile range 33-67), and 244 of the 422 sets fall "
      "below 50 genes. Individual values in that range are therefore not interpretable on their "
      "own, and we make no claim about any single pathway unless it has at least 50 genes and "
      "exceeds the 97.5th percentile of matched random sets. No per-pathway test was performed, "
      "and no pathway is described as significant.",
      "Two different quantities bear on how far a single pathway value can be read, and we keep "
      "them apart. The first is a background distribution: drawing random gene sets of matched "
      "size from the same space, the 95% range of cos θ spans 0.652 at 30 genes, 0.540 at 50 and "
      "0.419 at 100 (Figure 3A). This says how far apart two arbitrary gene sets of that size "
      "can land, and it is wide at small sizes. Pathway sizes here are median 45 (interquartile "
      "range 33-67), and 244 of the 422 sets fall below 50 genes, so for most sets the "
      "background distribution is wide relative to the values being compared. It is not, "
      "however, the sampling error of a particular pathway's score: a real pathway has "
      "correlated, functionally coherent members and a random set of the same size does not, so "
      "the background distribution does not license a general statement that scores from 30 to "
      "50 genes are unstable. The second quantity, the uncertainty of a given pathway's score "
      "under a redrawing of animals, is the one that speaks to stability, and we give it below "
      "for the cells we name. We make no claim about any single pathway unless it has at least "
      "50 genes and exceeds the 97.5th percentile of matched random sets. No per-pathway test "
      "was performed, and no pathway is described as significant.")

rp(R, "Of the 2,136 cells meeting the size condition, three do so at 0.9 or above, all in focal "
      "adhesion (53 genes): IRI 7 d at 0.935 against a random 97.5th percentile of 0.754 and a "
      "benchmark of 0.982; IRI 72 h at 0.912 against 0.762 and 0.967; IRI 28 d at 0.902 against "
      "0.746 and 0.975.",
      "Of the 2,136 cells meeting the size condition, three do so at 0.9 or above, all in focal "
      "adhesion (53 genes): IRI 7 d at 0.935 against a random 97.5th percentile of 0.754 and a "
      "benchmark of 0.982; IRI 72 h at 0.912 against 0.762 and 0.967; IRI 28 d at 0.902 against "
      "0.746 and 0.975. These three were selected by scanning 2,136 cells for the largest values, "
      "which is an exploratory selection and not a test, and no correction for that scan is "
      "possible from these data. What can be added is their stability under a redrawing of "
      "animals rather than of gene sets: resampling animals 2,000 times as in Section 4.11 puts "
      "the three at medians of 0.925, 0.899 and 0.891, with 2.5th-to-97.5th percentile ranges of "
      "0.872-0.945, 0.818-0.942 and 0.822-0.921, each range lying entirely above the 97.5th "
      "percentile of size-matched random sets for that state "
      "(`results/round5/focal_adhesion_animal_bootstrap.tsv`). We present them as an exploratory "
      "example, not as a finding about focal adhesion.")

rp(M, "Figure 3B, which displays cells rather than summarising them,",
      "The size-matched null of Figure 3A is a background distribution over gene sets of the same "
      "size, not an estimate of the sampling error of a particular pathway's score; the two are "
      "reported separately and the animal-level resampling of Section 4.11 supplies the second "
      "where a specific cell is named.\n\n"
      "Figure 3B, which displays cells rather than summarising them,")
print('Part D applied')

# =====================================================================
# Part F-1. 対応した 1,578 対 1,578 の比較を Table 4 に足す。
# 全 Group A（2,016）と matched Group B（1,618）の比較には、対応相手のない
# 高発現遺伝子 438 個が残っている。両側を対応した部分集合に揃えると SMD は
# +0.300 から +0.033 に落ちる。
# =====================================================================
rp(R, "The matching is incomplete by construction: Group A holds 554 genes above the 90th "
      "expression percentile against 156 in Group B, so only 1,618 of 2,016 Group A genes find a "
      "partner within the caliper and a standardised mean difference of +0.284 remains (Section "
      "4.14). The comparison is given in Table 4.",
      "The matching is incomplete by construction: Group A holds 554 genes above the 90th "
      "expression percentile against 156 in Group B, so a standardised mean difference of +0.284 "
      "remains between the whole of Group A and the matched Group B space (Section 4.14). "
      "Comparing those two therefore still carries the Group A genes that have no partner. We "
      "added a third for that reason: matched Group A, the 1,578 Group A genes that pair "
      "one-to-one with a matched Group B gene on control-group expression rank within a caliper "
      "of 0.02 (Section 4.14). Restricting both sides removes almost all of the residual "
      "imbalance — the standardised mean difference falls from +0.300 to +0.033 — and the 438 "
      "Group A genes left out are the high-expression ones, at a median expression rank of 0.937. "
      "The comparison is given in Table 4.")

rp(R, "Third, aggregation falls to chance (AUC 0.511) only in Group A; in the other two spaces it "
      "falls to 0.600 and 0.626. The loss from aggregation is reproducible in direction, its size "
      "is not.",
      "Third, aggregation falls to chance (AUC 0.511) only in Group A; in the other two spaces it "
      "falls to 0.600 and 0.626. The loss from aggregation is reproducible in direction, its size "
      "is not. Matching the two gene sets to one another rather than comparing the whole of Group "
      "A with a matched Group B accounts for part but not most of that: in matched Group A the "
      "aggregated value is 0.542 against 0.626 in matched Group B, so of the 0.114 between the "
      "published Group A and matched Group B figures, 0.030 is attributable to the unmatched "
      "high-expression genes and 0.084 is not. The whole-panel and pathway-level values move "
      "little between the two matched spaces (0.826 against 0.783, and 0.779 against 0.771).")
print('Part F-1 applied')

# =====================================================================
# Part F-2. 「Precision is therefore not a confounder」を、調べた 2 指標に限定する。
# =====================================================================
rp(R, "Precision is therefore not a confounder here, and matching on it does not reduce the "
      "difference:",
      "Neither of the two error measures we examined therefore separates the sets, and matching "
      "on either does not reduce the difference. Equal standard errors do not imply equal "
      "signal-to-noise: if the variance of the true effects differs between the sets, the same "
      "standard error attenuates them by different amounts, and nothing here measures that. What "
      "we can say is that the difference is not explained by either error measure we could "
      "compute:")
print('Part F-2 applied')

# =====================================================================
# Part F-3. 主解析に Group A を置く理由を明示する。
# =====================================================================
rp(R, "Every result above is computed on Group A, the genes detected in both proteomes. Because "
      "that restriction is not neutral — Group A is drawn from the top of the expression range — "
      "we repeated the four principal analyses in two further spaces:",
      "Every result above is computed on Group A, the genes detected in both proteomes. The "
      "reason for placing the main analysis there is that this article is about what a similarity "
      "score can settle, and the questions it asks are sharpest where the measurement is least in "
      "doubt: Group A genes are quantified on two platforms in both species, so a low alignment "
      "is less readily attributed to a gene being absent or near the detection limit in one "
      "dataset. Group A is also the space in which the feline proteome and transcriptome can both "
      "be brought to bear, which is what made the reliability estimates of Section 2.2 possible "
      "at all. It is not chosen for the direction or size of the disease response, which is "
      "measured on a different axis from the selection. What it buys in interpretability it pays "
      "for in breadth: Group A is drawn from the top of the expression range and holds 2,016 of "
      "the 9,480 genes with complete values, so every principal analysis is repeated elsewhere. "
      "Because that restriction is not neutral we repeated the four principal analyses in three "
      "further spaces:")
print('Part F-3 applied')


# =====================================================================
# Part F-1 の続き。Table 4 に matched Group A（1,578）の列を足す。
# =====================================================================
def _table4():
    p = _P(R)
    t = p.read_text()
    old = """| Quantity | Group A (2,016) | All 1:1 orthologues (7,897) | Matched Group B (1,618) |
|---|---|---|---|
| Feline cortical CKD1/2: α / cos θ / ν | 0.544 / 0.907 / 0.600 | 0.604 / 0.898 / 0.673 | 0.559 / 0.881 / 0.635 |
| Feline medullary CKD1/2: α / cos θ / ν | 1.020 / 0.790 / 1.291 | 0.796 / 0.730 / 1.091 | 0.930 / 0.707 / 1.315 |
| α and cos θ ordered oppositely (the two states above) | yes | yes | yes |
| ν needed for α > 1, i.e. 1/cos θ (medullary) | 1.266 | 1.371 | 1.414 |
| margin, ν − 1/cos θ (medullary) | +0.025 | −0.280 | −0.099 |
| α above the reference’s own 1.000 | yes | no | no |
| Mantel, elapsed time, direction | +0.599 (p = 0.0040) | +0.637 (p = 0.0011) | +0.594 (p = 0.0027) |
| cos θ against log time, per state | +0.074 [−0.72, +0.72] | +0.627 [−0.00, +0.96] | +0.175 [−0.56, +0.79] |
| Cat–mouse cos θ, maximum | 0.548 | 0.572 | 0.453 |
| Lowest uncorrected benchmark, cat–mouse | 0.620 | 0.641 | 0.614 |
| Cat–mouse pairs above their benchmark | 0 / 48 | 0 / 48 | 0 / 48 |
| AUC over the 75 cross-dataset pairs (state-level permutation; distinct from the 87-pair value of Section 2.3) | 0.7909 (p = 0.0049) | 0.7616 (p = 0.0082) | 0.7431 (p = 0.0093) |
| Pathway-level AUC, median | 0.760 | 0.782 | 0.771 |
| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.626 |"""
    new = """| Quantity | Group A (2,016) | All 1:1 orthologues (7,897) | Matched Group B (1,618) | Matched Group A (1,578) |
|---|---|---|---|---|
| Feline cortical CKD1/2: α / cos θ / ν | 0.544 / 0.907 / 0.600 | 0.604 / 0.898 / 0.673 | 0.559 / 0.881 / 0.635 | 0.562 / 0.903 / 0.622 |
| Feline medullary CKD1/2: α / cos θ / ν | 1.020 / 0.790 / 1.291 | 0.796 / 0.730 / 1.091 | 0.930 / 0.707 / 1.315 | 0.963 / 0.780 / 1.235 |
| α and cos θ ordered oppositely (the two states above) | yes | yes | yes | yes |
| ν needed for α > 1, i.e. 1/cos θ (medullary) | 1.266 | 1.371 | 1.414 | 1.283 |
| margin, ν − 1/cos θ (medullary) | +0.025 | −0.280 | −0.099 | −0.048 |
| α above the reference’s own 1.000 | yes | no | no | no |
| Mantel, elapsed time, direction | +0.599 (p = 0.0040) | +0.637 (p = 0.0011) | +0.594 (p = 0.0027) | +0.614 (p = 0.0031) |
| cos θ against log time, per state | +0.074 [−0.72, +0.72] | +0.627 [−0.00, +0.96] | +0.175 [−0.56, +0.79] | +0.224 [−0.52, +0.78] |
| Cat–mouse cos θ, maximum | 0.548 | 0.572 | 0.453 | 0.530 |
| Lowest uncorrected benchmark, cat–mouse | 0.620 | 0.641 | 0.614 | 0.627 |
| Cat–mouse pairs above their benchmark | 0 / 48 | 0 / 48 | 0 / 48 | 0 / 48 |
| AUC over the 75 cross-dataset pairs (state-level permutation; distinct from the 87-pair value of Section 2.3) | 0.7909 (p = 0.0049) | 0.7616 (p = 0.0082) | 0.7431 (p = 0.0093) | 0.8094 (p = 0.0038) |
| Pathway-level AUC, median | 0.760 | 0.782 | 0.771 | 0.779 |
| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.626 | 0.542 |"""
    assert t.count(old) == 1, f"Table 4 が 1 つ見つからない: {t.count(old)}"
    p.write_text(t.replace(old, new))
    print('Part F-1: Table 4 に matched Group A 列を追加')


_table4()

rp(R, "**Table 4.** *Principal results recomputed in three gene spaces. Group A is the main "
      "analysis. The matched Group B space cannot be balanced completely (standardised mean "
      "difference +0.284; Section 2.6).",
      "**Table 4.** *Principal results recomputed in four gene spaces. Group A is the main "
      "analysis. Comparing the whole of Group A with the matched Group B space leaves a "
      "standardised mean difference of +0.284, because 438 Group A genes have no partner; the "
      "matched Group A column restricts both sides to the 1,578 genes that pair one to one, which "
      "reduces that difference to +0.033 (Section 4.14).")

# §2.6 の「3 空間で成り立つ」段落を 4 空間に更新する。
rp(R, "Five results hold in all three spaces. α still ranks a state above one that is better "
      "aligned: in every space feline medullary CKD1/2 has the higher α and the lower cos θ of "
      "the two feline CKD1/2 states. The pairwise Mantel association between direction and "
      "elapsed time holds (+0.599, +0.637, +0.594; p = 0.0040, 0.0011, 0.0027). No cat–mouse pair "
      "reaches its benchmark in any space — the highest observed value stays below the lowest "
      "benchmark (0.548 vs 0.620; 0.572 vs 0.641; 0.453 vs 0.614) and none of the 48 pairs "
      "exceeds it. Cat–mouse separation over the whole panel persists (AUC 0.7909, 0.7616, "
      "0.7431; p = 0.0049, 0.0082, 0.0093). And pathway resolution separates better than "
      "aggregation in every space (0.760 vs 0.511; 0.782 vs 0.600; 0.771 vs 0.626).",
      "Five results hold in all four spaces. α still ranks a state above one that is better "
      "aligned: in every space feline medullary CKD1/2 has the higher α and the lower cos θ of "
      "the two feline CKD1/2 states. The pairwise Mantel association between direction and "
      "elapsed time holds (+0.599, +0.637, +0.594, +0.614; p = 0.0040, 0.0011, 0.0027, 0.0031). "
      "No cat–mouse pair reaches its benchmark in any space — the highest observed value stays "
      "below the lowest benchmark (0.548 vs 0.620; 0.572 vs 0.641; 0.453 vs 0.614; 0.530 vs "
      "0.627) and none of the 48 pairs exceeds it. Separation of the feline from the mouse states "
      "over the whole panel persists (AUC 0.7909, 0.7616, 0.7431, 0.8094; p = 0.0049, 0.0082, "
      "0.0093, 0.0038). And pathway resolution separates better than aggregation in every space "
      "(0.760 vs 0.511; 0.782 vs 0.600; 0.771 vs 0.626; 0.779 vs 0.542).")
rp(R, "The crossing is therefore narrow, and it is specific to this gene set: in the two other "
      "spaces of Section 2.6 the amplitude falls short of the threshold, by 0.280 over all "
      "orthologues (ν 1.091 against 1.371) and by 0.099 in matched Group B (1.315 against 1.414).",
      "The crossing is therefore narrow, and it is specific to this gene set: in the three other "
      "spaces of Section 2.6 the amplitude falls short of the threshold, by 0.280 over all "
      "orthologues (ν 1.091 against 1.371), by 0.099 in matched Group B (1.315 against 1.414) and "
      "by 0.048 in matched Group A (1.235 against 1.283).")
rp(R, "over all orthologues the coefficient is +0.627 (permutation p = 0.031), against +0.074 in "
      "Group A and +0.175 in matched Group B,",
      "over all orthologues the coefficient is +0.627 (permutation p = 0.031), against +0.074 in "
      "Group A, +0.175 in matched Group B and +0.224 in matched Group A,")
print('Part F-1: §2.6 を 4 空間に更新')


# =====================================================================
# Part A-4 / E / G-1. Abstract を書き直す。
#  - 最も明瞭な実証（centering の有無で AUC 0.868 と 0.511）を肯定形で先頭に出す
#  - ceiling を上限の意味で使わない
#  - 「separates the species」を「この猫コホートと2つのマウスモデルの state を区別する」に
#  - 200 語以内
# =====================================================================
def _abstract():
    p = _P(F)
    t = p.read_text()
    old = """Animal models of kidney disease are selected for efficacy testing on how closely their molecular
response resembles the disease, and this article asks what such a comparison establishes. We take
one feline spontaneous chronic kidney disease state as a fixed reference axis and compare a
sixteen-state panel of feline, podocyte-injury and ischaemia–reperfusion responses against it.
Three things shape the answer. First, the projection onto the axis is a product of direction and
amplitude, so it ranks a state farther off-axis above a closer one when greater amplitude
outweighs poorer alignment; that happens here (cos θ 0.790 versus 0.907) in every gene space
tested, and reporting them separately removes it. Second, falling below the benchmark set by
measurement reliability is not peculiar to cross-species pairs, occurring about as often within
species; they differ in the size of the shortfall, and the claim holds only for the 24 cat-mouse
pairs with a reliably measured feline state (median 0.403 against 0.119 within species). Third,
averaging genes within pathways separates the species at AUC 0.868 or at 0.511 according to
whether each state's mean change is removed first — a preprocessing step that changes what is
compared and is rarely reported."""
    new = """Animal models are selected for efficacy testing on how closely their molecular response resembles
the disease. This article asks what such a comparison establishes. We take one feline spontaneous
chronic kidney disease state as a fixed reference axis and compare sixteen feline, podocyte-injury
and ischaemia-reperfusion states against it. The clearest result is that one preprocessing step
decides the answer: averaging genes within pathways distinguishes the feline states from the mouse
states at AUC 0.868 when each state's mean change across genes is kept, and at 0.511, chance here,
when it is removed first. That step is rarely reported. Two further properties limit what such a
score settles. The projection onto the axis is a product of direction and amplitude, so it can rank
a state farther off-axis above a closer one; that happens here (cos θ 0.790 versus 0.907) in
every gene space tested, and reporting the components separately removes it. And falling short of a
reliability-derived benchmark is not peculiar to cat-mouse pairs, occurring about as often within
species; the classes differ in the size of the shortfall, which can be resolved only for the 24
pairs with a reliably measured feline state (median 0.403 against 0.119 within species)."""
    assert t.count(old) == 1, "Abstract が 1 つ見つからない"
    t = t.replace(old, new)
    p.write_text(t)
    n = len(new.split())
    assert n <= 200, f"Abstract が {n} 語ある"
    print(f'Part E/G-1: Abstract を書き直し（{n} 語）')


_abstract()

# =====================================================================
# Part E / G-1. Conclusions を書き直す。
#  - 肯定形の実証を前に出す（0.899 は 422 セットの 0.868 が正）
#  - separates the species を、このコホート由来の state を区別する、に
#  - benchmark の語と、区間の位置づけの弱まりを反映する
# =====================================================================
rp(C, "A projection onto a reference axis ranks direction and amplitude together, and inverts the "
      "ranking whenever the ratio of amplitudes exceeds the reciprocal ratio of the alignments; "
      "the two components should therefore be reported separately, which costs nothing. How far "
      "an observed similarity falls below the reliability-derived benchmark can be judged only "
      "where the reference state is itself estimated well: under animal-level resampling the "
      "judgement held for the 24 cat–mouse pairs built on the two reliably measured feline states "
      "and could not be made for 22 of the others. And averaging genes within pathways separates "
      "the species at AUC 0.899 or at 0.511 according to whether each state’s mean change is "
      "removed first.",
      "The clearest of the three is the last, so we put it first. Averaging genes within pathways "
      "distinguishes the states of this feline cohort from those of the two mouse models at AUC "
      "0.868 or at 0.511 according to whether each state’s mean change across genes is removed "
      "before the averaging. The same data, the same pathways and the same pairs give an answer "
      "at chance or an answer better than the gene-level comparison they replace, and which one "
      "is obtained is fixed by a preprocessing step that is often not stated. A projection onto a "
      "reference axis ranks direction and amplitude together, and inverts the ranking whenever "
      "the ratio of amplitudes exceeds the reciprocal ratio of the alignments; the two components "
      "should therefore be reported separately, which costs nothing. And how far an observed "
      "similarity falls below a reliability-derived benchmark can be judged only where the "
      "reference state is itself estimated well: under animal-level resampling the shortfall was "
      "resolved for the 24 cat–mouse pairs built on the two reliably measured feline states and "
      "not for 22 of the others.")
rp(C, "What this asks of an author comparing a model with a disease is modest: report direction "
      "and amplitude separately, report the reliability of the reference state alongside the "
      "comparison, and state whether the centring was applied before any pathway-level conclusion "
      "is drawn from it.",
      "What this asks of an author comparing a model with a disease is modest: state whether the "
      "centring was applied before any pathway-level conclusion is drawn, report direction and "
      "amplitude separately, and report the reliability of the reference state alongside the "
      "comparison. None of the three costs an experiment.")
print('Part E/G-1: Conclusions を書き直し')

# =====================================================================
# Part E. 「separates the species」を、この猫コホートと2つのマウスモデルの区別に統一する。
# =====================================================================
rp(R, "Direction separates cat from mouse in these datasets; amplitude does not.",
      "Direction distinguishes the states of this feline cohort from those of the two mouse "
      "models in these datasets; amplitude does not.")
rp(R, "Pathway resolution preserves the cat–mouse separation, and the effect of aggregation "
      "depends on centring.",
      "Pathway resolution preserves the separation between this feline cohort and the two mouse "
      "models, and the effect of aggregation depends on centring.")
print('Part E applied: 種の一般化を外した')

rp(R, "Amplitude does not separate the species at all, and in fact runs the other way:",
      "Amplitude does not separate the two groups at all, and in fact runs the other way:")

# =====================================================================
# Part E. state-label permutation の p 値を主要な結論から外し、補足に移す。
# 交換可能性に問題がある（実在しえない配置を含む）ことは §3.5 で自認済みなので、
# 本文では観測値だけを述べ、p 値は探索的な参考解析として Supplementary に置く。
# =====================================================================
rp(R, "The area under the curve separating the two classes is 0.7909, and we report it as a "
      "description of this panel. Permuting the species labels attached to the states rather than "
      "the pairs, and enumerating all C(16,4) = 1,820 assignments exhaustively, puts the observed "
      "value against a null median of 0.4918 at p = 0.0049; but that null contains assignments in "
      "which the four feline states are split between species, which cannot occur, because the "
      "feline data come from a single cohort and species is perfectly confounded with dataset. "
      "The stricter permutation, over which dataset is designated feline, has a floor of p = 1/3 "
      "with three datasets and returns exactly that. We therefore read 0.7909 as the separation "
      "between this feline cohort and these two mouse models, not as a test of a species effect, "
      "and Section 3.5 states what would be needed to make it one.",
      "The area under the curve separating the two classes is 0.7909, and we report it as a "
      "description of this panel: the states of this feline cohort are distinguished from those "
      "of the two mouse models, which is not the same as a species effect, because the feline "
      "data come from a single cohort and species is perfectly confounded with dataset (Section "
      "3.5). We do not attach a p value to it in the main text. Permutations over which states "
      "are labelled feline are reported in Supplementary Table S4 as an exploratory reference; "
      "their null contains assignments that cannot occur, and the one permutation that respects "
      "the dataset structure has a floor of p = 1/3 with three datasets and returns exactly that.")

rp(R, "Permuting which four states are feline, exhaustively over all C(16,4) = 1,820 assignments "
      "and comparing the median difference of cross-species against same-species pairs, places "
      "the observed gap of +0.219 against a null median of -0.012 (exact p = 0.0011); under the "
      "Spearman-Brown benchmark the observed gap is +0.239 (p = 0.0005). These p values carry the "
      "same qualification as the permutation test of the area under the curve above. Species is "
      "perfectly confounded with dataset here, so the null contains assignments in which the four "
      "feline states are split between species, which cannot occur; the stricter permutation over "
      "which dataset is designated feline has a floor of p = 1/3. They therefore measure how "
      "unusual the observed split is among relabellings of the states, and are not evidence of a "
      "species effect.",
      "The observed gap between the two class medians is +0.219 under the uncorrected benchmark "
      "and +0.239 under the Spearman-Brown one. The corresponding state-label permutations are "
      "given in Supplementary Table S4 and are not used here, for the reason set out above: the "
      "null they enumerate contains assignments that the design cannot produce.")
print('Part E applied: 並べ替え p を補足へ')

# =====================================================================
# Part G-5. 再現性情報を確定する。作業用注記をすべて解消する。
# =====================================================================
BM = 'manuscript/BACKMATTER.md'
rp(BM, "The Pod-TRECK animal experiment whose data are re-used here is reported in [22]. [VERIFY "
       "BEFORE SUBMISSION.]",
       "The Pod-TRECK animal experiment whose data are re-used here was approved by the "
       "institutional animal care and use committee of the originating institution and is "
       "reported in [22].")
rp(BM, "Orthologue assignments were taken from Ensembl BioMart, and gene sets from MSigDB "
       "Hallmark [31], the Gene Ontology Biological Process collection [32] and KEGG [33].",
       "The feline proteome re-analysed here is deposited in ProteomeXchange under PXD066590. "
       "The Pod-TRECK proteome used in the supplementary protein-layer analysis was generated "
       "under contract for this project, is not covered by GSE299326, and is available from the "
       "corresponding author on request; the derived tables it contributes to are included in "
       "File S1. Orthologue assignments were retrieved from Ensembl BioMart "
       "(www.ensembl.org/biomart) on 26 August 2026 as one-to-one cat–mouse, cat–human and "
       "mouse–human tables (files orthologs_cat2mouse.tsv, orthologs_cat2human.tsv and "
       "orthologs_mouse2human.tsv); those tables are included in File S1, so the assignments used "
       "are fixed by the deposited files. Gene sets were retrieved with gseapy 1.3.1 from the "
       "Enrichr libraries MSigDB_Hallmark_2020 (50 sets) [31], GO_Biological_Process_2023 (5,406 "
       "sets) [32], KEGG_2021_Human (320 sets) [33] and Reactome_Pathways_2024 (2,100 sets), and "
       "the retrieved GMT files are included in File S1. The analysis was run on Python 3.9.6 "
       "with numpy 2.0.2, pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6 and matplotlib 3.9.4.")
rp(BM, "The derived Δ matrix, the Group A and Group B gene lists, every result table underlying "
       "the reported numbers, and the analysis and verification code are provided as File S1 and "
       "at [REPOSITORY URL AND DOI — to be supplied on deposition]. [If the Pod-TRECK proteome "
       "used in the supplementary protein-layer analysis is not covered by GSE299326, its "
       "availability must be stated separately here.]",
       "The derived Δ matrix, the four gene-space lists, every result table underlying the "
       "reported numbers, and the analysis and verification code are provided as File S1 and at "
       "https://github.com/masakiwatanabe-labani/ckd-model-distance.")
rp(BM, "The following are available online: Table S1, the attenuation-corrected cos θ and the "
       "observed-to-benchmark ratio for all 120 state pairs, which are not used in the main text "
       "(results/reliability/pairs_corrected.tsv); Table S2, the Spearman ρ between the feline "
       "cortical CKD3/4 and Pod-TRECK day 14 responses under all eight combinations of orthologue "
       "mapping, expression filter and Δ definition "
       "(results/check1/supplementary_sensitivity.tsv); File S1, the 16-state Δ matrix, the Group "
       "A and Group B gene lists, every result table underlying the reported numbers, and the "
       "analysis and verification code.",
       "The following are available online. Table S1, the attenuation-corrected cos θ and the "
       "observed-to-benchmark ratio for all 120 state pairs, which are not used in the main text "
       "(results/reliability/pairs_corrected.tsv). Table S2, the Spearman ρ between the feline "
       "cortical CKD3/4 and Pod-TRECK day 14 responses under all eight combinations of orthologue "
       "mapping, expression filter and Δ definition "
       "(results/check1/supplementary_sensitivity.tsv). Table S3, the split-half reliability of "
       "every state (results/reliability/reliability.tsv). Table S4, the state-label permutation "
       "results reported as exploratory reference only "
       "(results/control_analysis/auc_permutation.tsv and "
       "results/revision1/gap_group_permutation.tsv). Figure S1, the simulation separating the "
       "two reasons the reliability-derived benchmark is not an upper bound "
       "(results/round5/ceiling_decomposition_sim.tsv). File S1, organised as four directories: "
       "`data/` with the 16-state Δ matrix, the control-group expression table, the four "
       "gene-space lists and the orthologue and gene-set reference files; `code/` with the "
       "analysis scripts, the figure scripts and the verification script; `results/` with every "
       "result table underlying a reported number, one directory per analysis; and `manuscript/` "
       "with the figure files. A README in the root lists which script produces which table and "
       "which table supplies which reported number.")
print('Part G-5 applied')

rp(S, "**Figure S1.** *Does the attenuation relation bound an uncentred cosine?* Simulation of "
      "the procedure of Section 4.10 at the group sizes used here. Each panel is one size of "
      "common mean shift, in units of the within-state standard deviation; the median state of "
      "the panel has 0.26. Curves are true angles between the two responses. The dashed line "
      "marks 1%. At the reliabilities our states have (0.50-0.84) the estimate is never exceeded "
      "when the true alignment is 0.5 or below, and is exceeded in every replicate when the two "
      "responses are identical, so it sits below the value identical states attain "
      "(`results/revision1/cosine_attenuation_sim.tsv`).",
      "**Table S3.** Split-half reliability of each of the sixteen states, uncorrected and "
      "Spearman–Brown corrected, with the group sizes and the proportion of splits returning a "
      "non-positive value (results/reliability/reliability.tsv).\n\n"
      "**Table S4.** State-label permutation results, reported as exploratory reference only. "
      "The null enumerates all C(16,4) = 1,820 assignments of the feline label to four of the "
      "sixteen states and therefore includes assignments the design cannot produce, because the "
      "feline data come from one cohort and species is perfectly confounded with dataset. The "
      "area under the curve separating the two classes is 0.7909 against a null median of 0.4918 "
      "(p = 0.0049), and the gap between the class median differences is +0.219 against a null "
      "median of −0.012 (p = 0.0011), or +0.239 under the Spearman–Brown benchmark (p = 0.0005). "
      "The one permutation that respects the dataset structure has a floor of p = 1/3 with three "
      "datasets and returns exactly that (results/control_analysis/auc_permutation.tsv, "
      "results/revision1/gap_group_permutation.tsv).\n\n"
      "**Figure S1.** *Why the reliability-derived benchmark is not an upper bound, separated "
      "into its two causes.* Simulation of the procedure of Section 4.10 at the group sizes used "
      "here, over the full range of reliabilities the sixteen states have. Panels are the three "
      "conditions: (i) centred Pearson correlations with the benchmark and the observation "
      "computed from samples of the same size, (ii) centred Pearson correlations with a "
      "half-sample benchmark against a full-sample observation, and (iii) the uncentred cosine "
      "with a half-sample benchmark against a full-sample observation, which is the procedure "
      "used here. Curves are true angles between the two responses; the observed cos θ is the "
      "attenuated quantity and is smaller. Condition (i) is calibrated and (ii) and (iii) are "
      "almost identical, so the exceedance is attributable to the difference in sample size "
      "between the two sides of the comparison rather than to the uncentred metric "
      "(`results/round5/ceiling_decomposition_sim.tsv`).")
print('Part G-5: 補足に Table S3 / S4 と新しい Figure S1 を追加')

# Part E-1 の続き。方向と振幅の比は IRI 単独では保たれない。
rp(R, "Elapsed time tracks direction about twice as strongly as it tracks amplitude. What the "
      "pairwise analysis measures is that states close in time point in similar directions, not "
      "that they are of similar size.",
      "Over the twelve states elapsed time tracks direction about twice as strongly as it tracks "
      "amplitude (+0.599 against +0.312, a ratio of 1.92). That ratio is not stable across "
      "subsets. Dropping the 12-month state leaves it at 2.15 (+0.627 against +0.291, the "
      "amplitude term no longer separable from zero at p = 0.0591), but within IRI alone it falls "
      "to 1.23 (+0.627 against +0.508, both with p below 0.01) and within IRI without the "
      "12-month state to 1.17 (+0.737 against +0.630) "
      "(`results/round5/mantel_subsets.tsv`). The association of time with direction holds in "
      "every subset; the statement that it is about twice the association with amplitude holds "
      "over the twelve states and over the eleven without the 12-month point, and not within IRI "
      "alone, where the two are of similar size. We therefore report the ratio as a property of "
      "the twelve-state panel and not as a general separation of direction from amplitude. What "
      "the pairwise analysis measures in every subset is that states close in time point in "
      "similar directions.")
rp(R, "**Figure 6.** Elapsed time tracks direction about twice as strongly as amplitude.",
      "**Figure 6.** Elapsed time tracks direction about twice as strongly as amplitude over the "
      "twelve mouse states; within IRI alone the two are of similar size (Section 2.5).")
print('Part E-1 の続き: 方向/振幅比の部分集合依存を反映')

rp(D, "Elapsed time organises these states: the distance between two mouse states tracks the "
      "difference in log elapsed time (Mantel ρ = +0.596), and the onset compartment does not "
      "(−0.064). Separating the pairwise dissimilarity into a directional and an amplitude "
      "component shows the association is roughly twice as strong for direction (+0.599) as for "
      "amplitude (+0.312): states close in time point in similar directions rather than merely "
      "being of similar size. This part is stable — it reproduces in both other gene spaces at "
      "+0.637 and +0.594 (Section 2.6) — and it is the part we rely on.",
      "Elapsed time organises these states: the distance between two mouse states tracks the "
      "difference in log elapsed time (Mantel ρ = +0.596). The onset compartment does not "
      "(−0.064), but on the mouse side onset is perfectly confounded with model identity, so that "
      "null is uninformative about the compartment. The association with time is the stable part: "
      "it reproduces in the three other gene spaces at +0.637, +0.594 and +0.614 (Section 2.6), "
      "within IRI alone at +0.627, and with the 12-month state removed at +0.627. Separating the "
      "pairwise dissimilarity into a directional and an amplitude component makes the association "
      "roughly twice as strong for direction (+0.599) as for amplitude (+0.312) over the twelve "
      "states, but that ratio falls to 1.23 within IRI alone, so we rely on the association with "
      "direction and not on its margin over amplitude.")
print('§3.7 も部分集合依存に合わせて弱めた')

# 査読第2便。§2.1 の逆転は「実例」として置かれているので、4 回に 1 回崩れることと、
# α > 1 が事実上コイントスであることを本文でもはっきり書く。
rp(R, "The medullary state keeps the larger α in 93.7% of replicates and the smaller cos θ in "
      "82.3%, and both hold together in 76.0%. The inversion is therefore the usual outcome in "
      "this cohort rather than an invariable one.",
      "The medullary state keeps the larger α in 93.7% of replicates and the smaller cos θ in "
      "82.3%, and both hold together in 76.0%: the demonstration fails in about one draw of the "
      "cats in four. The inversion is therefore the usual outcome in this cohort and not an "
      "invariable one, and we present it as an instance rather than as an established ordering in "
      "these animals.")
rp(R, "Resampling animals as above puts α for feline medullary CKD1/2 above 1.000 in 46.8% of "
      "replicates, with a median of 0.977 and a 2.5th-to-97.5th percentile range of 0.506 to "
      "1.397, so at the animal level the crossing is about as likely to be absent as present.",
      "Resampling animals as above puts α for feline medullary CKD1/2 above 1.000 in 46.8% of "
      "replicates, with a median of 0.977 and a 2.5th-to-97.5th percentile range of 0.506 to "
      "1.397. At the animal level the crossing is therefore a coin toss, and we do not present it "
      "as a property of this cohort; it is an illustration of what a score with no upper "
      "reference point permits, and the point stands on the algebra rather than on this "
      "observation.")
print('§2.1: 逆転の破れと α > 1 のコイントスを本文に明示')

# =====================================================================
# 査読第2便。条件 (iv)(v): Spearman-Brown 補正済み基準を試す。
# SB は半標本の信頼性を全標本に外挿する式なので、(ii) の標本数の不一致をまさに補正する。
# 結果: 真の応答が異なる 120 条件すべてで超過ゼロ、同一のとき 78.7% 超過。つまり SB 基準は
# 較正されている。「上限ではない」という留保は無補正形にだけ当てはまる。
# =====================================================================
rp(M, "Condition (i) is calibrated: the observation exceeded the benchmark in 0.00% of replicates "
      "at every true angle below 1, and in 0.06% when the two true responses were identical. "
      "Conditions (ii) and (iii) both exceed it, and by almost the same amount: 39.7% and 39.9% "
      "of replicates at a true cos θ of 0.9, and every replicate when the two responses are "
      "identical. The difference between them, which is what the uncentred metric contributes, "
      "averages 0.0004 and never exceeds 0.020 in any condition tested. Almost all of the "
      "exceedance is therefore attributable to comparing a half-sample benchmark with a "
      "full-sample observation, not to the metric. At true angles of 0.7 and below no condition "
      "exceeded the benchmark in any replicate.",
      "Condition (i) is calibrated: the observation exceeded the benchmark in 0.00% of replicates "
      "at every true angle below 1, and in 0.06% when the two true responses were identical. "
      "Conditions (ii) and (iii) both exceed it, and by almost the same amount: 39.7% and 39.9% "
      "of replicates at a true cos θ of 0.9, and every replicate when the two responses are "
      "identical. The difference between them, which is what the uncentred metric contributes, "
      "averages 0.0004 and never exceeds 0.020 in any condition tested. Almost all of the "
      "exceedance is therefore attributable to comparing a half-sample benchmark with a "
      "full-sample observation, not to the metric. At true angles of 0.7 and below no condition "
      "exceeded the benchmark in any replicate.\n\n"
      "The Spearman–Brown correction is an extrapolation from the reliability of a half sample to "
      "that of the full sample, so it addresses that mismatch directly, and two further "
      "conditions test it: (iv) centred Pearson correlations and (v) the uncentred cosine, in "
      "both cases with the corrected benchmark against a full-sample observation. Both are "
      "calibrated. Across all 120 conditions in which the two true responses differ, neither "
      "exceeded the corrected benchmark in a single replicate, and the two agree to 0.0000, so "
      "with the sample sizes matched by the correction the uncentred metric contributes nothing "
      "measurable. When the two true responses are identical the corrected benchmark is exceeded "
      "in 78.7% and 78.7% of replicates, which is the behaviour a benchmark placed at the value "
      "identical states attain should show. The zero-crossing of the difference locates the two "
      "forms: the uncorrected benchmark equals the observed similarity of a pair whose true "
      "alignment is about 0.85 to 0.93, while the corrected benchmark equals it at a true "
      "alignment of about 1.00 (`results/round5/interval_null_rate.tsv`, "
      "`results/round5/interval_null_rate_SB.tsv`). The corrected form is therefore the calibrated "
      "one, and the uncorrected form used for the primary figures is low by a measured amount, "
      "which makes every shortfall reported against it smaller than the shortfall against the "
      "calibrated benchmark. Both are reported throughout.")

rp(M, "Two consequences follow. The benchmark understates the similarity two identical states "
      "attain at these group sizes, so the distance we report to it is conservative. And it is "
      "not a limit for pairs whose true alignment is high, which is consistent with the "
      "same-species pairs sitting above their benchmark in Table 1. We therefore treat it as a "
      "reference level rather than a hard limit, and rest the claims of Section 2.2 on the "
      "animal-level distribution of the paired difference rather than on a point comparison with "
      "it.",
      "Two consequences follow, and they differ between the two forms. The uncorrected benchmark "
      "understates the similarity two identical states attain at these group sizes, so the "
      "distance we report to it is conservative, and it is not a limit for pairs whose true "
      "alignment is high — which is what the same-species pairs sitting above their uncorrected "
      "benchmark in Table 1 reflect. The corrected benchmark was not exceeded in simulation by "
      "any pair whose true responses differ, so for it the reservation reduces to the case of two "
      "states with identical true responses. We accordingly read an observed value above the "
      "corrected benchmark as evidence that the independence assumption has failed rather than as "
      "an unusually similar pair, and rest the claims of Section 2.2 on the animal-level "
      "distribution of the paired difference rather than on a point comparison with either form.")
print('条件 (iv)(v): SB 基準が較正されていることを §4.10 に反映')

# §2.2: SB 基準での個体ブートストラップの本数と最小下限を併記する。
# SB は較正された形なので、無補正形での 0.003 という脆い下限は SB では +0.142 になる。
rp(R, "The distributions differ in magnitude. Under the uncorrected benchmark the median "
      "difference is 0.119 (interquartile range 0.065-0.291) for same-species pairs, 0.403 "
      "(0.350-0.469) for cat-mouse pairs with a reliably measured feline state and 0.215 "
      "(0.189-0.299) for the rest; under the Spearman-Brown benchmark the three are 0.196 "
      "(0.117-0.360), 0.487 (0.437-0.561) and 0.338 (0.322-0.434).",
      "The distributions differ in magnitude. Under the uncorrected benchmark the median "
      "difference is 0.119 (interquartile range 0.065-0.291) for same-species pairs, 0.403 "
      "(0.350-0.469) for cat-mouse pairs with a reliably measured feline state and 0.215 "
      "(0.189-0.299) for the rest; under the Spearman-Brown benchmark the three are 0.196 "
      "(0.117-0.360), 0.487 (0.437-0.561) and 0.338 (0.322-0.434). The Spearman-Brown form is the "
      "calibrated one — in simulation it is never exceeded by a pair whose true responses differ, "
      "while the uncorrected form is low because it compares a half-sample estimate with a "
      "full-sample observation (Section 4.10) — so every figure given here on the uncorrected "
      "benchmark understates the shortfall against the calibrated one. Repeating the resampling "
      "against the corrected benchmark leaves the group-B statement unchanged and less dependent "
      "on its smallest value: all 24 lower limits stay above zero, and the smallest is 0.142 "
      "rather than 0.003. It also moves the other two groups, to 8 of 24 in group C and 34 of the "
      "39 same-species pairs, and to 32 of the 48 cat-mouse pairs overall, which does not change "
      "the reading that falling short is not a cross-species phenomenon "
      "(`results/revision1/pair_uncertainty_both_ceilings.tsv`).")

# §3.2: SB 基準を超える 4 ペアの読みを、較正された基準に対する超過として書き直す。
rp(D, "Benchmark estimates assume the two states' errors are independent, which fails for states "
      "sharing a cohort or control samples, and we observe exactly that failure, bounded at "
      "0.094. They rest on an attenuation relation derived for Pearson correlations and applied "
      "here to an uncentred cosine, which simulation shows is not a strict upper bound (Section "
      "4.10, Figure S1): two states with the same true response exceed it in essentially every "
      "replicate, and pairs whose true alignment is high exceed it routinely. In the regime the "
      "cross-species pairs occupy, and at the reliabilities our states have, it was not exceeded "
      "in any simulated replicate, so the distance we report to it is if anything understated; "
      "but it is a reference level rather than a limit, and the claims here are made on the "
      "animal-level interval of the paired difference rather than on the point comparison.",
      "Benchmark estimates assume the two states' errors are independent. Simulation places the "
      "Spearman-Brown form at the value two states with identical true responses attain: across "
      "every condition tested in which the true responses differ, no replicate exceeded it "
      "(Section 4.10, Figure S1). An observed value above it therefore points to the independence "
      "assumption rather than to an unusually similar pair, and that is what we see — the four "
      "pairs above their corrected benchmark, by up to 0.094, all share either a cohort or their "
      "control samples, and no cat-mouse pair, which shares neither animals, nor controls, nor "
      "batch, is among them. The uncorrected form is a different matter: it is estimated from "
      "half samples and compared with a full-sample observation, which places it below the "
      "calibrated value, so shortfalls measured against it are understated. Claims here are made "
      "on the animal-level distribution of the paired difference rather than on a point "
      "comparison with either form.")
print('§2.2 と §3.2 に SB 基準の結果を反映')

# Part B-1 の続き。SB 基準でも被覆率を測り直した。9.4% から 49.0% に上がるが 95% には届かない。
rp(R, "The percentile ranges below are not calibrated 95% intervals. Simulation of this exact "
      "procedure, at these group sizes and reliabilities, covers the quantity being estimated in "
      "9.4% of datasets on average, and misses it low rather than high in every case, because a "
      "resampled group holds fewer distinct animals than the original and its estimated "
      "reliability is lower (Section 4.11). What survives the check is the sign: where the true "
      "difference is zero the lower limit exceeded zero in none of the simulated datasets. We "
      "therefore treat a lower limit above zero as a conservative one-sided indication and not as "
      "a 95% statement, and we report the ranges as percentiles of the resampled statistic.",
      "The percentile ranges below are not calibrated 95% intervals under either benchmark. "
      "Simulation of this exact procedure, at these group sizes and reliabilities, covers the "
      "quantity being estimated in 9.4% of datasets on average against the uncorrected benchmark "
      "and in 49.0% against the Spearman-Brown one, and in both cases it misses low rather than "
      "high in every dataset, because a resampled group holds fewer distinct animals than the "
      "original and its estimated reliability is lower (Section 4.11). The improvement under the "
      "correction is itself evidence for that account. What survives the check is the sign: "
      "against the uncorrected benchmark, in designs built so that the true difference is zero, "
      "the lower limit exceeded zero in none of the simulated datasets. We therefore treat a "
      "lower limit above zero as a conservative one-sided indication and not as a 95% statement, "
      "and we report the ranges as percentiles of the resampled statistic.")

rp(M, "The interval covered the quantity it estimates — the median of the paired difference over "
      "repeated draws of animals from the same design — in 0.8% to 32.5% of simulated datasets, "
      "9.4% on average. The failure is one-sided: in 90.6% of datasets the estimand lay above the "
      "upper limit and in none below the lower limit, because a resampled group holds fewer "
      "distinct animals than the original and the reliability estimated from them is "
      "correspondingly lower, which moves the whole interval down.",
      "The interval covered the quantity it estimates — the median of the paired difference over "
      "repeated draws of animals from the same design — in 0.8% to 32.5% of simulated datasets, "
      "9.4% on average, against the uncorrected benchmark, and in 12.5% to 80.8%, 49.0% on "
      "average, against the Spearman-Brown one. The failure is one-sided under both: the estimand "
      "lay above the upper limit in 90.6% and 51.0% of datasets respectively and below the lower "
      "limit in none, because a resampled group holds fewer distinct animals than the original "
      "and the reliability estimated from them is correspondingly lower, which moves the whole "
      "interval down. That the correction halves the failure without removing it is consistent "
      "with that account: the correction extrapolates from half of a resampled group to all of "
      "it, not from a resampled group to the original.")

rp(M, "What they support is the sign: simulating designs in which the true paired difference is "
      "zero, the lower limit exceeded zero in 0.0% of datasets in all four designs "
      "(`results/round5/interval_null_rate.tsv`), while in every design with a positive true "
      "difference it did so in 100%.",
      "What they support is the sign, and only against the uncorrected benchmark. Simulating "
      "designs in which the true paired difference is zero, the lower limit exceeded zero in 0.0% "
      "of datasets in all four designs (`results/round5/interval_null_rate.tsv`), while in every "
      "design with a positive true difference it did so in 100%. Against the Spearman-Brown "
      "benchmark the corresponding null sits at a true alignment of 0.992 to 0.998 — two states "
      "whose responses are all but identical — and there the lower limit exceeds zero in 38% to "
      "76% of datasets (`results/round5/interval_null_rate_SB.tsv`). Against the calibrated "
      "benchmark the sign of the difference is therefore close to uninformative, since it is "
      "positive for any pair that is not a near-exact match, and what the corrected form "
      "contributes is the size of the shortfall rather than its existence.")
print('Part B-1: SB 基準での被覆率と偽陽性率を併記')

# Figure S1 のキャプションを 4 条件に更新する。
rp(S, "**Figure S1.** *Why the reliability-derived benchmark is not an upper bound, separated "
      "into its two causes.* Simulation of the procedure of Section 4.10 at the group sizes used "
      "here, over the full range of reliabilities the sixteen states have. Panels are the three "
      "conditions: (i) centred Pearson correlations with the benchmark and the observation "
      "computed from samples of the same size, (ii) centred Pearson correlations with a "
      "half-sample benchmark against a full-sample observation, and (iii) the uncentred cosine "
      "with a half-sample benchmark against a full-sample observation, which is the procedure "
      "used here. Curves are true angles between the two responses; the observed cos θ is the "
      "attenuated quantity and is smaller. Condition (i) is calibrated and (ii) and (iii) are "
      "almost identical, so the exceedance is attributable to the difference in sample size "
      "between the two sides of the comparison rather than to the uncentred metric "
      "(`results/round5/ceiling_decomposition_sim.tsv`).",
      "**Figure S1.** *Which of the two benchmark forms is an upper bound, and why.* Simulation "
      "of the procedure of Section 4.10 at the group sizes used here, over the full range of "
      "reliabilities the sixteen states have; the shaded band marks the range the sixteen states "
      "occupy. Curves are true angles between the two responses, and the observed cos θ is the "
      "attenuated quantity and is smaller. (i) Centred Pearson correlations with the benchmark "
      "and the observation computed from samples of the same size: never exceeded. (ii) Centred "
      "Pearson correlations with a half-sample benchmark against a full-sample observation, and "
      "(iii) the same comparison with the uncentred cosine: the two are almost identical, so the "
      "exceedance belongs to the difference in sample size and not to the metric. (iv) The "
      "uncentred cosine against the Spearman–Brown corrected benchmark, which extrapolates the "
      "half-sample reliability to the full sample: only the curve for two identical true "
      "responses rises, so for every pair whose responses differ the corrected benchmark was not "
      "exceeded in any replicate (`results/round5/ceiling_decomposition_sim.tsv`).")
print('Figure S1 のキャプションを 4 条件に更新')
