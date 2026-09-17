"""第7ラウンドの本文反映。

Part 0  文献を 34 件体系に（Ito et al. を [26] として挿入、旧 26-33 を 27-34 へ）
Part 1  区間に基づく判定を効果量の記述に置き換え、独立性の断定を削除する
Part 2  交換可能性を満たさない p 値を本文と Table 4 から外す
Part 3  §4.4 の Δ の名称と尺度の記述を正す
Part 4  マッチ空間の遺伝子数と SMD を 1 つの対応表にそろえる
Part 5  AUC がペアの分離であることを明示し、0.511 の読みを弱める
Part 6  再現性情報を確定する
"""
import re
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from edit_helper import rp

R, D, M = 'manuscript/RESULTS.md', 'manuscript/DISCUSSION.md', 'manuscript/METHODS.md'
C, F, S = 'manuscript/CONCLUSIONS.md', 'manuscript/FRONTMATTER.md', 'manuscript/SUPPLEMENTARY.md'
BM, I, REF = 'manuscript/BACKMATTER.md', 'manuscript/INTRODUCTION.md', 'manuscript/REFERENCES.md'
MD = _P('manuscript')

# =====================================================================
# Part 0. 文献体系を 34 件にそろえる。
# 本文側は rev9 で Ito et al. が [26] に入っており、Markdown 側は 33 件のままだった。
# 受け渡しのたびに +1 シフトを当てるのをやめ、Markdown 側を 34 件体系にする。
# =====================================================================
ITO = ("26. Ito, K.; Watanabe, M.; Kawamoto, K.; Kono, Y.; Sasaki, T.; Nakagawa, R.; Kawashima, "
       "Y.; Sasaki, N. Integrated proteomic–transcriptomic profiling identifies NR3C2 "
       "up-regulation and altered filtered-plasma-protein handling in chronic podocyte injury. "
       "FEBS J. 2026, accepted for publication.")
BODY = ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
        "CONCLUSIONS.md", "SUPPLEMENTARY.md", "BACKMATTER.md")


def shift_references():
    ref = MD / REF.split('/')[-1]
    t = ref.read_text()
    if "Ito, K.; Watanabe, M." in t:
        print('Part 0: 文献はすでに 34 件体系')
        return
    # 旧 33 -> 34 の順に後ろから繰り上げる（先に小さい番号を動かすと衝突する）
    for n in range(33, 25, -1):
        pat = re.compile(rf"^{n}\.\s", re.M)
        assert len(pat.findall(t)) == 1, f"文献 [{n}] が 1 件でない"
        t = pat.sub(f"{n + 1}. ", t)
    pat = re.compile(r"^27\.\s", re.M)          # 旧 26（Mantel）は 27 になっている
    assert len(pat.findall(t)) == 1
    t = pat.sub(ITO + "\n27. ", t, count=1)
    ref.write_text(t)

    # 本文の引用番号を追随させる。26 以上を +1 する。
    n_cit = 0
    for f in BODY:
        p = MD / f
        s = p.read_text()

        def bump(m):
            nonlocal n_cit
            inner = m.group(1)
            if not re.fullmatch(r"[0-9][0-9,–— -]*", inner):
                return m.group(0)
            out, changed = [], False
            for tok in re.split(r"([,–—-])", inner):
                if tok.isdigit():
                    v = int(tok)
                    if v >= 26:
                        v += 1
                        changed = True
                    out.append(str(v))
                else:
                    out.append(tok)
            if changed:
                n_cit += 1
            return "[" + "".join(out) + "]"

        s = re.sub(r"\[([^\]]+)\]", bump, s)
        p.write_text(s)
    print(f'Part 0: 文献を 34 件体系に（本文の引用 {n_cit} 箇所をシフト）')


shift_references()

# §4.6 と Data Availability に Pod-TRECK プロテオームの出典を書く
rp(M, "Group A is defined by detectability in the proteome alone, never by disease direction. "
      "Detection means non-missing in at least 50% of the samples of that dataset.",
      "Group A is defined by detectability in the proteome alone, never by disease direction. "
      "Detection means non-missing in at least 50% of the samples of that dataset. The feline "
      "cortical and medullary proteomes are those of the feline cohort [21]; the Pod-TRECK "
      "proteome is data acquired in the accepted study [26], is not part of GSE299326, and is "
      "described in the Data Availability Statement. The per-gene detection calls behind the "
      "definition, and the proteome tables they were made from, are provided in File S1 so that "
      "the membership of Group A can be checked gene by gene.")
rp(BM, "The Pod-TRECK proteome used in the supplementary protein-layer analysis was generated "
       "under contract for this project, is not covered by GSE299326, and is available from the "
       "corresponding author on request; the derived tables it contributes to are included in "
       "File S1.",
       "The Pod-TRECK proteome used to define Group A and in the supplementary protein-layer "
       "analysis is author-held data acquired in the accepted study [26]. It is not covered by "
       "GSE299326; during review it is supplied to the editorial office and to reviewers on "
       "request, and the derived tables and the per-gene detection calls it contributes to are "
       "included in File S1.")

# =====================================================================
# Part 1. 区間に基づく判定を、記述的な効果量に置き換える。
# 被覆率は無補正 9.4% / SB 49.0%、SB では真の差がゼロでも下限がゼロを超える割合が
# 38-76%。この数字のもとで「何組で解決した」を主要な根拠にはできない。
# =====================================================================
rp(R, "The percentile ranges below are not calibrated 95% intervals under either benchmark. "
      "Simulation of this exact procedure, at these group sizes and reliabilities, covers the "
      "quantity being estimated in 9.4% of datasets on average against the uncorrected benchmark "
      "and in 49.0% against the Spearman-Brown one, and in both cases it misses low rather than "
      "high in every dataset, because a resampled group holds fewer distinct animals than the "
      "original and its estimated reliability is lower (Section 4.11). The improvement under the "
      "correction is itself evidence for that account. What survives the check is the sign: "
      "against the uncorrected benchmark, in designs built so that the true difference is zero, "
      "the lower limit exceeded zero in none of the simulated datasets. We therefore treat a "
      "lower limit above zero as a conservative one-sided indication and not as a 95% statement, "
      "and we report the ranges as percentiles of the resampled statistic.",
      "We report the resampled difference as an effect size and not as a decision. Simulation of "
      "this exact procedure, at these group sizes and reliabilities, covers the quantity being "
      "estimated in 9.4% of datasets on average against the uncorrected benchmark and in 49.0% "
      "against the Spearman–Brown one, and against the corrected benchmark the lower limit "
      "exceeds zero in 38% to 76% of datasets built so that the true difference is zero (Section "
      "4.11). A range with those properties cannot carry a per-pair verdict, whatever it is "
      "called, so no count of how many pairs were resolved is used as evidence below. What we "
      "report is the distribution of the difference: its median and interquartile range within "
      "each group, and the spread across pairs.")

rp(R, "Against the Spearman–Brown benchmark the difference stays above zero in 32 of the 48 "
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
      "zero in 1 of the 48.",
      "The size of the difference follows the split-half reliability of the feline state in the "
      "pair. Against the Spearman–Brown benchmark the 24 cross-species pairs whose feline state "
      "is cortical CKD3/4 or medullary CKD1/2 (reliabilities 0.841 and 0.851) have a median "
      "difference of 0.487 with an interquartile range of 0.437 to 0.561, and the 24 built on "
      "cortical CKD1/2 (0.568) or medullary CKD3/4 (0.616, from five animals) a median of 0.338 "
      "(0.323 to 0.434). The 39 same-species pairs sharing no controls sit at 0.196 (0.117 to "
      "0.360). Under the uncorrected benchmark the three are 0.403 (0.350 to 0.469), 0.215 (0.189 "
      "to 0.299) and 0.119 (0.065 to 0.291): the ordering is the same and only the level moves. "
      "As a sensitivity analysis we also record where the resampled range lies relative to zero "
      "under the uncorrected benchmark, which is the one form for which the simulation of Section "
      "4.11 found no false positive: there the lower limit is above zero for all 24 of the first "
      "group, for 2 of the second and for 24 of the 39 same-species pairs, and a "
      "leave-one-animal-out jackknife puts the difference below zero in 1 of the 48. That "
      "guarantee is limited to the four designs and the reliabilities simulated and we do not "
      "extend it beyond them.")

rp(R, "**The positive claim is limited to the pairs with a reliably measured feline state.** Of "
      "the three groups, only group B - the 24 cat-mouse pairs built on feline cortical CKD3/4 or "
      "medullary CKD1/2 - supports a statement that the observed value falls short of the "
      "benchmark. Group C, the 24 pairs built on the two less reliably measured feline states, "
      "has 16 ranges containing zero and a median difference of 0.338 against 0.196 for "
      "same-species pairs; we do not offer it as evidence in either direction, and no claim below "
      "rests on it.",
      "**How far the difference can be read depends on how well the feline state is measured.** "
      "The 24 cat–mouse pairs built on feline cortical CKD3/4 or medullary CKD1/2 carry the "
      "largest and least dispersed difference of the three groups, and they are the pairs the "
      "sensitivity analysis above places above zero under the one benchmark for which that "
      "statement was checked. The 24 built on the two less reliably measured feline states have a "
      "smaller median and a wider spread that reaches below zero; we do not offer that group as "
      "evidence in either direction, and no claim below rests on it.")

rp(R, "**What separates the classes is the size of the shortfall, not whether there is one.** "
      "Taken as a yes-or-no question, falling short of the benchmark is not a cross-species "
      "phenomenon at all: the range excludes zero for 34 of the 39 same-species pairs (87%) and "
      "for 32 of the 48 cat-mouse pairs (67%). The distributions differ in magnitude.",
      "**What separates the classes is the size of the shortfall, not whether there is one.** "
      "Falling short of the benchmark is not a cross-species phenomenon: same-species pairs fall "
      "short too, and by an amount whose distribution overlaps that of the cross-species pairs. "
      "The distributions differ in magnitude.")
print('Part 1-1: 判定を効果量の記述に置き換え')

# =====================================================================
# Part 1-2. 「ベンチマーク超過は独立性の破れを示す」という断定を削除する。
# 独立誤差のシミュレーションでも、真の応答が同一なら SB を 78.7% で超える。
# 4 ペアが共有構造を持つという観察は記述として残す。
# =====================================================================
rp(R, "Four pairs sit above their Spearman–Brown benchmark, by up to 0.094 (Figure 2C). The "
      "benchmark is a geometric mean of two reliabilities and assumes the two states’ errors are "
      "independent; all four of these pairs share either a cohort or their control samples, and "
      "no cross-species pair — which shares neither animals, nor controls, nor batch — does. We "
      "take 0.094 as a lower bound on the systematic error of the benchmark for pairs with shared "
      "structure, and note that it is not observed where the independence assumption holds.",
      "Four pairs sit above their Spearman–Brown benchmark, by up to 0.094 (Figure 2B). All four "
      "share either a cohort or their control samples; no cross-species pair — which shares "
      "neither animals, nor controls, nor batch — is among them. We report that as an "
      "observation and do not read it as establishing a cause. The benchmark assumes the two "
      "states’ errors are independent, and dependent errors would inflate an observed value "
      "towards its benchmark; but the simulation of Section 4.10 shows that two states with "
      "identical true responses exceed the corrected benchmark in 78.7% of replicates with "
      "independent errors, so exceeding it is not on its own evidence that independence has "
      "failed.")

rp(M, "We accordingly read an observed value above the corrected benchmark as evidence that the "
      "independence assumption has failed rather than as an unusually similar pair, and rest the "
      "claims of Section 2.2 on the animal-level distribution of the paired difference rather "
      "than on a point comparison with either form.",
      "An observed value above the corrected benchmark is therefore consistent with two states "
      "whose true responses are nearly identical as well as with a failure of independence, and "
      "we do not use an exceedance to decide between them. The claims of Section 2.2 rest on the "
      "animal-level distribution of the paired difference rather than on a point comparison with "
      "either form.")

rp(D, "Benchmark estimates assume the two states' errors are independent. Simulation places the "
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
      "comparison with either form.",
      "Benchmark estimates assume the two states' errors are independent. Simulation places the "
      "Spearman-Brown form at the value two states with identical true responses attain: across "
      "every condition tested in which the true responses differ, no replicate exceeded it, while "
      "two states with identical true responses exceeded it in 78.7% of replicates even with "
      "independent errors (Section 4.10, Figure S1). Exceeding the benchmark is therefore not on "
      "its own evidence that independence has failed. We record that the four pairs above their "
      "corrected benchmark, by up to 0.094, all share either a cohort or their control samples, "
      "and that no cat-mouse pair, which shares neither animals, nor controls, nor batch, is "
      "among them, without drawing a cause from it. The uncorrected form is a different matter: "
      "it is estimated from half samples and compared with a full-sample observation, which "
      "places it below the calibrated value, so shortfalls measured against it are understated. "
      "Claims here are made on the animal-level distribution of the paired difference rather than "
      "on a point comparison with either form.")
print('Part 1-2: 独立性の断定を削除')

# Part 1-1 の続き。Abstract・Conclusions・Figure 2 キャプション・§3.2 を、判定ではなく
# 効果量の言い方にそろえる。
rp(F, "the classes differ in the size of the shortfall, which can be resolved only for the 24 "
      "pairs with a reliably measured feline state (median 0.487 against 0.196 within species).",
      "the classes differ in the size of the shortfall, which is largest for the 24 pairs with a "
      "reliably measured feline state (median 0.487 against 0.196 within species).")
rp(C, "And how far an observed similarity falls below a reliability-derived benchmark can be "
      "judged only where the reference state is itself estimated well: under animal-level "
      "resampling the shortfall was resolved for the 24 cat–mouse pairs built on the two reliably "
      "measured feline states and not for 16 of the others.",
      "And how far an observed similarity falls below a reliability-derived benchmark is "
      "interpretable only where the reference state is itself estimated well: under animal-level "
      "resampling the shortfall is largest and least dispersed for the 24 cat–mouse pairs built "
      "on the two reliably measured feline states, and the resampled range is too poorly "
      "calibrated to support a per-pair verdict.")
rp(R, "and the paired difference is plotted with its 95% interval. Resampling three animals from "
      "three leaves fewer than two distinct animals about one replicate in nine, and no benchmark "
      "can be formed there; every interval shown is conditional on the replicates in which one "
      "can (a median of 88.2% of draws per cat–mouse pair).",
      "and the paired difference against the Spearman–Brown benchmark is plotted with the 2.5th "
      "to 97.5th percentile range of the resampled values. That range is not a calibrated "
      "interval — simulation of the same procedure covers the quantity it estimates in 49.0% of "
      "datasets and its lower limit exceeds zero in 38% to 76% of datasets built so that the true "
      "difference is zero (Section 4.11) — so it is shown to convey spread and no count of pairs "
      "is read from it. Resampling three animals from three leaves fewer than two distinct "
      "animals about one replicate in nine, and no benchmark can be formed there; every range "
      "shown is conditional on the replicates in which one can (a median of 88.2% of draws per "
      "cat–mouse pair).")
rp(R, "An interval containing zero means the difference could not be resolved at this sample "
      "size, not that the pair reached its benchmark.",
      "A range containing zero does not mean the pair reached its benchmark.")
rp(D, "it happens about as often within species as across them - the range excludes zero for 87% "
      "of same-species pairs and 67% of cat-mouse pairs - and what separates the two classes is "
      "how large the shortfall is.",
      "it happens within species as well, and the two distributions overlap - and what separates "
      "the two classes is how large the shortfall is.")
print('Part 1-1: Abstract / Conclusions / Figure 2 / §3.2 をそろえた')

# =====================================================================
# Part 2. 交換可能性を満たさない p 値を本文と Table 4 から外す。
# §2.2 は「本文では p を付けない」と宣言しているのに §2.6・Table 4・§3.2 で使っていた。
# AUC はこのパネルの記述量として残す。
# =====================================================================
rp(R, "Separation of the feline from the mouse states over the whole panel persists (AUC 0.7909, "
      "0.7616, 0.7431, 0.8094; p = 0.0049, 0.0082, 0.0093, 0.0038).",
      "Separation of the feline from the mouse states over the whole panel persists (AUC 0.7909, "
      "0.7616, 0.7431, 0.8094), reported as a description of each panel and without a p value, "
      "for the reason given in Section 2.2.")
rp(R, "| AUC over the 75 cross-dataset pairs (state-level permutation; distinct from the 87-pair "
      "value of Section 2.3) | 0.7909 (p = 0.0049) | 0.7616 (p = 0.0082) | 0.7431 (p = 0.0093) | "
      "0.8094 (p = 0.0038) |",
      "| AUC over the 75 cross-dataset pairs (distinct from the 87-pair value of Section 2.3; no "
      "permutation test, Section 2.2) | 0.7909 | 0.7616 | 0.7431 | 0.8094 |")
rp(D, "Separating those two classes by permuting species labels over states gives AUC 0.7909 "
      "(exact p = 0.0049); the value of 0.813 quoted in Section 3.3 is a different quantity, "
      "taken over the wider set of 87 pairs that share no controls rather than the 75 "
      "cross-dataset pairs, and it carries no permutation test.",
      "Those two classes separate at AUC 0.7909 over the 75 cross-dataset pairs; the value of "
      "0.813 quoted in Section 3.3 is a different quantity, taken over the wider set of 87 pairs "
      "that share no controls. Neither carries a permutation test, for the reason given in "
      "Section 2.2.")
rp(D, "Species and dataset are therefore completely confounded, and although permuting species "
      "labels over states gives p = 0.0049, that null includes assignments in which a single "
      "dataset is split between species — configurations that cannot occur. The strictest test, "
      "permuting which dataset is designated feline, has a floor of p = 1/3 with three datasets "
      "and returns exactly that.",
      "Species and dataset are therefore completely confounded. Permuting species labels over "
      "states enumerates assignments in which a single dataset is split between species — "
      "configurations that cannot occur — so the resulting null is not the null of any hypothesis "
      "the design can test, and we attach no p value to the separation anywhere in the article. "
      "The strictest permutation, over which dataset is designated feline, has a floor of p = 1/3 "
      "with three datasets and returns exactly that, which is the honest statement of what this "
      "design can support.")

# 補足 Table S4 の扱いを、妥当な検定として読ませない形にする。
rp(S, "**Table S4.** State-label permutation results, reported as exploratory reference only. The "
      "null enumerates all C(16,4) = 1,820 assignments of the feline label to four of the sixteen "
      "states and therefore includes assignments the design cannot produce, because the feline "
      "data come from one cohort and species is perfectly confounded with dataset.",
      "**Table S4.** State-label permutation results. These are not a valid test of any "
      "hypothesis and no p value from them is used anywhere in the article; they are tabulated "
      "only so that a reader who wants the numbers behind that statement can see them. The null "
      "enumerates all C(16,4) = 1,820 assignments of the feline label to four of the sixteen "
      "states, which includes assignments the design cannot produce, because the feline data come "
      "from one cohort and species is perfectly confounded with dataset; a null containing "
      "impossible configurations does not become usable by being called exploratory.")
print('Part 2: 交換可能性を満たさない p 値を本文と Table 4 から削除')

rp(D, "The median difference is 0.196 within species against 0.487 and 0.338 in the two "
      "cross-species groups, and permuting which four states are feline places the observed "
      "difference in medians, +0.239, against a null median of -0.012 (exact p = 0.0005) - a p "
      "value that, like the one attached to the area under the curve, is computed over "
      "relabellings of states drawn from datasets that are themselves confounded with species, "
      "and so describes the observed split rather than testing a species effect.",
      "The median difference is 0.196 within species against 0.487 and 0.338 in the two "
      "cross-species groups, a gap between the class medians of +0.239. We attach no p value to "
      "it: the only permutation available relabels states drawn from datasets that are themselves "
      "confounded with species, so its null contains configurations the design cannot produce "
      "(Section 3.5, Table S4).")
rp(D, "a median difference of 0.215 there, against 0.119 among same-species pairs, is not "
      "something we read as a shortfall.",
      "a median difference of 0.338 there, against 0.196 among same-species pairs, is not "
      "something we read as a shortfall.")
print('Part 2: §3.2 の群間差からも p を外した')

# =====================================================================
# Part 3-1. Δ の名称を「平均 log2 発現量の差」に統一し、尺度の記述を正す。
# 主要解析はコサイン類似度であって順位相関ではない。
# =====================================================================
rp(M, "All comparisons here are between states measured on the same scale within a dataset or are "
      "rank-based across datasets, and the sensitivity analysis of Table S2 gives the effect of "
      "removing the expression filter that excludes the genes nearest the offset (main effect "
      "−0.067).",
      "We therefore call Δ the difference of mean log2 expression throughout, and not a log2 fold "
      "change. The principal comparisons are cosine similarities between Δ vectors, which are not "
      "invariant to a monotone change of scale, so the offset and the expression filter can move "
      "them; rank-based statistics appear only in the Mantel analysis of Section 4.9 and in the "
      "Spearman ρ of Section 2.4. Table S5 reports how the principal cosine quantities and the "
      "pathway-aggregated AUC move when the offset and the low-expression filter are varied, and "
      "Table S2 the effect on the Spearman ρ.")
rp(M, "This is a standardised effect size rather than a log2 fold change and it reorders genes.",
      "This is a standardised effect size rather than a difference of mean log2 expression, and "
      "it reorders genes.")
rp(R, "Δlog2 fold changes against its own control group (Section 4.4).",
      "differences of mean log2 expression against its own control group (Section 4.4).")
rp(R, "we fixed it in advance — one-to- one Ensembl orthologues, a group-mean FPKM ≥ 1 filter, and "
      "a plain log2 fold change — and report all eight combinations of the three choices as a "
      "sensitivity analysis (Table S2).",
      "we fixed it in advance — one-to-one Ensembl orthologues, a group-mean FPKM ≥ 1 filter, and "
      "the plain difference of mean log2 expression — and report all eight combinations of the "
      "three choices as a sensitivity analysis (Table S2), with the effect on the principal "
      "cosine quantities in Table S5.")
rp(R, "The main effects on the feline-versus-Pod-TRECK Spearman ρ were −0.196 for replacing the "
      "log2 fold change with a per-gene standardised effect, −0.067 for removing the expression "
      "filter, and +0.012 for using upper-case symbol matching in place of orthologues.",
      "The main effects on the feline-versus-Pod-TRECK Spearman ρ were −0.196 for replacing the "
      "difference of mean log2 expression with a per-gene standardised effect, −0.067 for "
      "removing the expression filter, and +0.012 for using upper-case symbol matching in place "
      "of orthologues.")
print('Part 3-1: Δ の名称と尺度の記述を統一')

# =====================================================================
# Part 5-1. AUC はペアの分離であってペア以外の分類ではない、と読めるようにする。
# =====================================================================
rp(F, "The clearest result is that one preprocessing step decides the answer: averaging genes "
      "within pathways distinguishes the feline states from the mouse states at AUC 0.868 when "
      "each state's mean change across genes is kept, and at 0.511, chance here, when it is "
      "removed first.",
      "The clearest result is that one preprocessing step decides the answer. Scoring every pair "
      "of states and asking how well the score tells a feline-mouse pair from a same-species "
      "pair, pathway-averaged data separate the two kinds of pair at AUC 0.868 when each state's "
      "mean change across genes is kept, and at 0.511 when it is removed first.")
rp(C, "Averaging genes within pathways distinguishes the states of this feline cohort from those "
      "of the two mouse models at AUC 0.868 or at 0.511 according to whether each state’s mean "
      "change across genes is removed before the averaging.",
      "Averaging genes within pathways separates cat–mouse pairs of states from same-species "
      "pairs at AUC 0.868 or at 0.511 according to whether each state’s mean change across genes "
      "is removed before the averaging. The quantity classified is the pair, not the state.")
rp(R, "To compare the two on equal terms we used the same quantity and the same pairs throughout: "
      "within each pathway, the cos θ of all 120 state pairs with no reference axis fixed, "
      "restricted as in Section 2.2 to the 87 pairs that share no control samples, and scored as "
      "the AUC separating the 39 within-species from the 48 cross-species pairs.",
      "To compare the two on equal terms we used the same quantity and the same pairs throughout: "
      "within each pathway, the cos θ of all 120 state pairs with no reference axis fixed, "
      "restricted as in Section 2.2 to the 87 pairs that share no control samples, and scored as "
      "the AUC separating the 39 within-species from the 48 cross-species pairs. The unit being "
      "classified is the pair of states and not the state, so every AUC in this section is the "
      "separability of two kinds of pair and not the accuracy of assigning a state to a species.")

# =====================================================================
# Part 5-2. 0.511 を「識別できない」と断定しない。
# =====================================================================
rp(R, "the separation falls to 0.511, indistinguishable from chance in this gene set (Figure 4B).",
      "the separation falls to 0.511 (Figure 4B). That is close to the 0.5 a coin would give, and "
      "with 39 against 48 pairs the estimate is not precise enough to establish that no "
      "discriminative information survives; what it does establish is that the separation carried "
      "by the same pathways before the centring is largely gone.")
rp(D, "It decides the answer completely. With centring, within- and cross-species pairs separate "
      "at AUC 0.511, which is chance, both classes having risen to a median of about 0.83.",
      "It decides the answer. With centring, within- and cross-species pairs separate at AUC "
      "0.511, close to the 0.5 a coin would give, both classes having risen to a median of about "
      "0.83; we read that as the loss of most of the separation and not as proof that none "
      "remains.")
print('Part 5-1 / 5-2: AUC の対象と 0.511 の読みを直した')

# Abstract を 200 語以内に戻す（Part 5-1 の書き換えで 213 語になった）。
rp(F, "The clearest result is that one preprocessing step decides the answer. Scoring every pair "
      "of states and asking how well the score tells a feline-mouse pair from a same-species "
      "pair, pathway-averaged data separate the two kinds of pair at AUC 0.868 when each state's "
      "mean change across genes is kept, and at 0.511 when it is removed first. That step is "
      "rarely reported. Two further properties limit what such a score settles.",
      "One preprocessing step decides the answer. Asked how well a score tells a cat-mouse pair "
      "of states from a same-species pair, pathway-averaged data separate the two kinds of pair "
      "at AUC 0.868 when each state's mean change across genes is kept, and at 0.511 when it is "
      "removed first. That step is rarely reported. Two further properties limit such a score.")
rp(F, "And falling short of a reliability-derived benchmark is not peculiar to cat-mouse pairs, "
      "occurring about as often within species; the classes differ in the size of the shortfall, "
      "which is largest for the 24 pairs with a reliably measured feline state (median 0.487 "
      "against 0.196 within species).",
      "And falling short of a reliability-derived benchmark is not peculiar to cat-mouse pairs, "
      "occurring within species too; the classes differ in the size of the shortfall, largest for "
      "the 24 pairs with a reliably measured feline state (median 0.487 against 0.196).")

# Table S5 を新設する（Part 3-2 の尺度感度分析）。
rp(S, "**Table S3.** Split-half reliability of each of the sixteen states",
      "**Table S5.** Sensitivity of the principal cosine quantities and the pathway-aggregated "
      "AUC to the input scale. The Δ matrix was rebuilt with the offset of the log2 transform set "
      "to 1, 0.5 and 0.1 and the group-mean FPKM filter set to 1, 0 and 5, giving nine "
      "combinations, and four quantities recomputed in each: the α, cos θ and ν of the two feline "
      "CKD1/2 states, the gene-level AUC separating within-species from cross-species pairs over "
      "the 87 pairs sharing no controls, and the same AUC after aggregation to pathway means with "
      "and without per-state centring "
      "(results/round7/scale_sensitivity.tsv).\n\n"
      "**Table S3.** Split-half reliability of each of the sixteen states")
print('Abstract を 200 語以内に戻し、Table S5 を新設')

rp(BM, "Table S3, the split-half reliability of every state "
       "(results/reliability/reliability.tsv).",
       "Table S3, the split-half reliability of every state "
       "(results/reliability/reliability.tsv). Table S5, the sensitivity of the principal cosine "
       "quantities and the pathway-aggregated AUC to the offset of the log2 transform and the "
       "low-expression filter (results/round7/scale_sensitivity.tsv).")
print('背表紙に Table S5 を追加')

# =====================================================================
# Part 5-3. 中心化による AUC 変化を個体単位で再標本化した結果を本文に入れる。
# 向きは 98.6% で再現するが、水準はばらつく（中心化後 AUC の 95% 範囲 0.432-0.737）。
# =====================================================================
rp(R, "The same averaging therefore improves the distinction or removes it according to whether "
      "each state’s average change across genes was taken out first.",
      "The same averaging therefore improves the distinction or removes it according to whether "
      "each state’s average change across genes was taken out first.\n\n"
      "Because this is the central result we asked how far it depends on which animals were "
      "sampled. Resampling animals 1,000 times with the design intact — states sharing a control "
      "group receive the same resampled controls, and the feline cortical and medullary samples, "
      "being the same cats, are resampled together (Section 4.11) — and recomputing both AUCs "
      "from each resample, the uncentred value exceeds the centred one in 98.6% of replicates. "
      "The direction of the effect is therefore robust to the exchange of animals. Its size is "
      "less so: the drop has a median of 0.302 and a 2.5th-to-97.5th percentile range of 0.037 to "
      "0.438 against the 0.356 observed, the uncentred AUC ranges from 0.724 to 0.899 and the "
      "centred one from 0.432 to 0.737. In about a fifth of resamples the centred value is not "
      "near 0.5. We therefore state the effect as a direction that reproduces and a magnitude "
      "that is estimated with this much uncertainty, and we do not treat 0.868 and 0.511 as "
      "stable to the exchange of animals "
      "(`results/round7/centering_auc_animal_bootstrap.tsv`).")
rp(D, "It decides the answer. With centring, within- and cross-species pairs separate at AUC "
      "0.511, close to the 0.5 a coin would give, both classes having risen to a median of about "
      "0.83; we read that as the loss of most of the separation and not as proof that none "
      "remains.",
      "It decides the answer. With centring, within- and cross-species pairs separate at AUC "
      "0.511, close to the 0.5 a coin would give, both classes having risen to a median of about "
      "0.83; we read that as the loss of most of the separation and not as proof that none "
      "remains. Under an animal-level resampling that preserves the design, the uncentred value "
      "exceeds the centred one in 98.6% of replicates, so the direction is not an artefact of "
      "which cats and mice were sampled; the size of the gap is less certain, the drop running "
      "from 0.037 to 0.438 across resamples (Section 2.3).")
print('Part 5-3: 中心化効果の個体再標本化を反映')

# =====================================================================
# Part 6. 再現性情報を確定する。
# =====================================================================
rp(BM, "The derived Δ matrix, the four gene-space lists, every result table underlying the "
       "reported numbers, and the analysis and verification code are provided as File S1 and at "
       "https://github.com/masakiwatanabe-labani/ckd-model-distance.",
       "The derived Δ matrix, the gene-space lists, the per-gene proteome detection calls behind "
       "the definition of Group A, every result table underlying the reported numbers, and the "
       "analysis, simulation and verification code are provided as File S1 and at "
       "https://github.com/masakiwatanabe-labani/ckd-model-distance. The commit used for the "
       "results reported here is tagged in that repository, and the tagged tree is deposited in "
       "Zenodo; the concept DOI issued on deposition is given with the repository link there.")
rp(BM, "File S1, organised as four directories: `data/` with the 16-state Δ matrix, the "
       "control-group expression table, the four gene-space lists and the orthologue and gene-set "
       "reference files; `code/` with the analysis scripts, the figure scripts and the "
       "verification script; `results/` with every result table underlying a reported number, one "
       "directory per analysis; and `manuscript/` with the figure files. A README in the root "
       "lists which script produces which table and which table supplies which reported number.",
       "File S1, organised as four directories. `data/` holds the 16-state Δ matrix, the "
       "control-group expression table, the gene-space lists, the per-gene proteome detection "
       "calls and the proteome tables they were made from, and the orthologue and gene-set "
       "reference files. `code/` holds the analysis scripts, the simulation scripts, the figure "
       "scripts and the verification script. `results/` holds every result table underlying a "
       "reported number, one directory per analysis. `manuscript/` holds the figure files. A "
       "README in the root lists which script produces which table and which table supplies which "
       "reported number, and gives the software versions and the one command that re-runs the "
       "verification.")
print('Part 6: Data Availability と File S1 の記述を確定')

# =====================================================================
# Part 3-2. 尺度を変えたときに主要なコサイン量と集約後 AUC がどう動くかを本文に書く。
# =====================================================================
rp(R, "The main effects on the feline-versus-Pod-TRECK Spearman ρ were −0.196 for replacing the "
      "difference of mean log2 expression with a per-gene standardised effect, −0.067 for "
      "removing the expression filter, and +0.012 for using upper-case symbol matching in place "
      "of orthologues.",
      "The main effects on the feline-versus-Pod-TRECK Spearman ρ were −0.196 for replacing the "
      "difference of mean log2 expression with a per-gene standardised effect, −0.067 for "
      "removing the expression filter, and +0.012 for using upper-case symbol matching in place "
      "of orthologues.\n\n"
      "Because the principal comparisons are cosine similarities rather than rank statistics, we "
      "also rebuilt the Δ matrix with the offset of the log2 transform set to 1, 0.5 and 0.1 and "
      "the group-mean FPKM filter set to 1, 0 and 5, and recomputed the quantities the "
      "conclusions rest on in each of the nine combinations (Table S5). The gene-level AUC "
      "separating within-species from cross-species pairs moves from 0.813 to between 0.813 and "
      "0.821; the pathway-aggregated AUC without centring from 0.868 to between 0.863 and 0.877, "
      "and with centring from 0.511 to between 0.510 and 0.545; the difference the centring makes "
      "stays between 0.326 and 0.357 against 0.356 at the chosen settings. The ordering inversion "
      "of Section 2.1 holds in all nine, as does the crossing of α = 1. The offset therefore "
      "shifts these quantities by a few hundredths and changes none of the conclusions, which is "
      "what the argument of Section 4.4 requires but does not by itself guarantee "
      "(`results/round7/scale_sensitivity.tsv`).")
print('Part 3-2: 尺度感度の結果を §2.4 に')


# =====================================================================
# Part 4. マッチ空間を、保存された 1 つの対応表から定義し直した値に差し替える。
# 旧 matched Group B（1,618）は対応表が残っておらず、方向を変えて作り直すと件数が
# 変わっていた（A→B 1,632、B→A 1,578、対称制限 1,509）。決定的な 1 回の
# マッチングで両空間を 1,632 遺伝子にそろえ、対応表を results/round7 に保存した。
# =====================================================================
def _table4():
    p = _P(R)
    t = p.read_text()
    old = """| Quantity | Group A (2,016) | All 1:1 orthologues (7,897) | Matched Group B (1,618) | Matched Group A (1,578) |
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
| Lowest Spearman–Brown benchmark, cat–mouse | 0.765 | 0.780 | 0.760 | 0.770 |
| Cat–mouse pairs above their benchmark | 0 / 48 | 0 / 48 | 0 / 48 | 0 / 48 |
| AUC over the 75 cross-dataset pairs (distinct from the 87-pair value of Section 2.3; no permutation test, Section 2.2) | 0.7909 | 0.7616 | 0.7431 | 0.8094 |
| Pathway-level AUC, median | 0.760 | 0.782 | 0.771 | 0.779 |
| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.626 | 0.542 |"""
    new = """| Quantity | Group A (2,016) | All 1:1 orthologues (7,897) | Matched Group B (1,632) | Matched Group A (1,632) |
|---|---|---|---|---|
| Feline cortical CKD1/2: α / cos θ / ν | 0.544 / 0.907 / 0.600 | 0.604 / 0.898 / 0.673 | 0.546 / 0.876 / 0.623 | 0.566 / 0.901 / 0.629 |
| Feline medullary CKD1/2: α / cos θ / ν | 1.020 / 0.790 / 1.291 | 0.796 / 0.730 / 1.091 | 0.987 / 0.734 / 1.345 | 0.947 / 0.776 / 1.221 |
| α and cos θ ordered oppositely (the two states above) | yes | yes | yes | yes |
| ν needed for α > 1, i.e. 1/cos θ (medullary) | 1.266 | 1.371 | 1.362 | 1.289 |
| margin, ν − 1/cos θ (medullary) | +0.025 | −0.280 | −0.017 | −0.069 |
| α above the reference’s own 1.000 | yes | no | no | no |
| Mantel, elapsed time, direction | +0.599 (p = 0.0040) | +0.637 (p = 0.0011) | +0.646 (p = 0.0009) | +0.613 (p = 0.0031) |
| cos θ against log time, per state | +0.074 [−0.72, +0.72] | +0.627 [−0.00, +0.96] | +0.147 [−0.61, +0.77] | +0.210 [−0.55, +0.78] |
| Cat–mouse cos θ, maximum | 0.548 | 0.572 | 0.438 | 0.534 |
| Lowest Spearman–Brown benchmark, cat–mouse | 0.765 | 0.780 | 0.765 | 0.769 |
| Cat–mouse pairs above their benchmark | 0 / 48 | 0 / 48 | 0 / 48 | 0 / 48 |
| AUC over the 75 cross-dataset pairs (distinct from the 87-pair value of Section 2.3; no permutation test, Section 2.2) | 0.7909 | 0.7616 | 0.7708 | 0.8102 |
| Pathway-level AUC, median | 0.760 | 0.782 | 0.780 | 0.781 |
| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.624 | 0.529 |"""
    assert t.count(old) == 1, f"Table 4 が 1 つ見つからない: {t.count(old)}"
    p.write_text(t.replace(old, new))
    print('Part 4: Table 4 のマッチ空間を 1,632 / 1,632 に差し替え')


_table4()

rp(R, "**Table 4.** *Principal results recomputed in four gene spaces. Group A is the main "
      "analysis. Comparing the whole of Group A with the matched Group B space leaves a "
      "standardised mean difference of +0.284, because 438 Group A genes have no partner; the "
      "matched Group A column restricts both sides to the 1,578 genes that pair one to one, which "
      "reduces that difference to +0.033 (Section 4.14).",
      "**Table 4.** *Principal results recomputed in four gene spaces. Group A is the main "
      "analysis. The two matched columns are the two sides of one saved 1:1 matching on "
      "control-group expression rank, so they hold the same 1,632 genes each and are balanced to "
      "a standardised mean difference of −0.017, against +1.236 between the unmatched sets; the "
      "384 Group A genes with no partner are the high-expression ones (Section 4.14).")

# ---- §2.6 と §4.14 を、保存された対応表の数値にそろえる ----
rp(R, "Because that restriction is not neutral we repeated the four principal analyses in three "
      "further spaces: all one-to-one orthologues with no proteome restriction (7,897 genes), and "
      "a Group B subset matched to Group A on control-group expression (1,618 genes). The "
      "matching is incomplete by construction: Group A holds 554 genes above the 90th expression "
      "percentile against 156 in Group B, so a standardised mean difference of +0.284 remains "
      "between the whole of Group A and the matched Group B space (Section 4.14). Comparing those "
      "two therefore still carries the Group A genes that have no partner. We added a third for "
      "that reason: matched Group A, the 1,578 Group A genes that pair one-to-one with a matched "
      "Group B gene on control-group expression rank within a caliper of 0.02 (Section 4.14). "
      "Restricting both sides removes almost all of the residual imbalance — the standardised "
      "mean difference falls from +0.300 to +0.033 — and the 438 Group A genes left out are the "
      "high-expression ones, at a median expression rank of 0.937. The comparison is given in "
      "Table 4.",
      "Because that restriction is not neutral we repeated the four principal analyses in three "
      "further spaces. The first is all one-to-one orthologues with no proteome restriction "
      "(7,897 genes). The other two are the two sides of a single one-to-one matching of Group A "
      "genes to Group B genes on control-group expression rank, caliper 0.02, computed without "
      "randomisation and saved as a pair list (Section 4.14). The matching is incomplete by "
      "construction: Group A holds 476 genes above the 90th expression percentile against 176 in "
      "Group B, so 384 of the 2,016 Group A genes find no partner, and those are the "
      "high-expression ones, at a median expression rank of 0.949. The 1,632 pairs that do form "
      "give a matched Group A space and a matched Group B space of 1,632 genes each, balanced on "
      "expression to a standardised mean difference of −0.017 against +1.236 between the "
      "unmatched sets. Because both spaces come from the same pair list, the two columns of Table "
      "4 are directly comparable. The comparison is given in Table 4.")

rp(R, "Third, aggregation falls to chance (AUC 0.511) only in Group A; in the other two spaces it "
      "falls to 0.600 and 0.626. The loss from aggregation is reproducible in direction, its size "
      "is not. Matching the two gene sets to one another rather than comparing the whole of Group "
      "A with a matched Group B accounts for part but not most of that: in matched Group A the "
      "aggregated value is 0.542 against 0.626 in matched Group B, so of the 0.114 between the "
      "published Group A and matched Group B figures, 0.030 is attributable to the unmatched "
      "high-expression genes and 0.084 is not. The whole-panel and pathway-level values move "
      "little between the two matched spaces (0.826 against 0.783, and 0.779 against 0.771).",
      "Third, aggregation falls closest to chance in Group A (AUC 0.511); in the other three "
      "spaces it falls to 0.600, 0.624 and 0.529. The loss from aggregation is reproducible in "
      "direction, its size is not. The two matched spaces isolate what the gene selection "
      "contributes, because they hold the same number of genes matched pair by pair on "
      "expression: matched Group A gives 0.529 and matched Group B 0.624. Of the 0.113 between "
      "the Group A and matched Group B figures, 0.018 is attributable to the 384 unmatched "
      "high-expression genes and 0.095 is not. The whole-panel and pathway-level values move "
      "little between the two matched spaces (0.830 against 0.802, and 0.781 against 0.780).")

rp(R, "The pairwise Mantel association between direction and elapsed time holds (+0.599, +0.637, "
      "+0.594, +0.614; p = 0.0040, 0.0011, 0.0027, 0.0031). No cat–mouse pair reaches its "
      "benchmark in any space — the highest observed value stays below the lowest Spearman–Brown "
      "benchmark (0.548 vs 0.765; 0.572 vs 0.780; 0.453 vs 0.760; 0.530 vs 0.770) and none of the "
      "48 pairs exceeds it, under either form of the benchmark.",
      "The pairwise Mantel association between direction and elapsed time holds (+0.599, +0.637, "
      "+0.646, +0.613; p = 0.0040, 0.0011, 0.0009, 0.0031). No cat–mouse pair reaches its "
      "benchmark in any space — the highest observed value stays below the lowest Spearman–Brown "
      "benchmark (0.548 vs 0.765; 0.572 vs 0.780; 0.438 vs 0.765; 0.534 vs 0.769) and none of the "
      "48 pairs exceeds it, under either form of the benchmark.")

rp(R, "over all orthologues the coefficient is +0.627 (permutation p = 0.031), against +0.074 in "
      "Group A, +0.175 in matched Group B and +0.224 in matched Group A,",
      "over all orthologues the coefficient is +0.627 (permutation p = 0.031), against +0.074 in "
      "Group A, +0.147 in matched Group B and +0.210 in matched Group A,")

rp(R, "Separation of the feline from the mouse states over the whole panel persists (AUC 0.7909, "
      "0.7616, 0.7431, 0.8094), reported as a description of each panel and without a p value, "
      "for the reason given in Section 2.2.",
      "Separation of the feline from the mouse states over the whole panel persists (AUC 0.7909, "
      "0.7616, 0.7708, 0.8102), reported as a description of each panel and without a p value, "
      "for the reason given in Section 2.2.")

rp(R, "The crossing is therefore narrow, and it is specific to this gene set: in the three other "
      "spaces of Section 2.6 the amplitude falls short of the threshold, by 0.280 over all "
      "orthologues (ν 1.091 against 1.371), by 0.099 in matched Group B (1.315 against 1.414) and "
      "by 0.048 in matched Group A (1.235 against 1.283).",
      "The crossing is therefore narrow, and it is specific to this gene set: in the three other "
      "spaces of Section 2.6 the amplitude falls short of the threshold, by 0.280 over all "
      "orthologues (ν 1.091 against 1.371), by 0.017 in matched Group B (1.345 against 1.362) and "
      "by 0.069 in matched Group A (1.221 against 1.289).")

rp(M, "The same matching builds the gene space used in Section 2.6: one-to-one orthologues "
      "outside Group A were matched to Group A on the mean control-group expression percentile, "
      "caliper 0.02, without replacement, taking the scarcest Group A genes first. Group A cannot "
      "be matched completely, because it holds 554 genes above the 90th expression percentile "
      "against 156 in Group B; 1,618 of 2,016 Group A genes find a partner and a standardised "
      "mean difference of +0.284 remains. Results in that space are read with that residual "
      "imbalance in mind.",
      "The same estimator builds the two matched gene spaces used in Section 2.6, with one "
      "difference that matters for reproducibility: the matching is run once, without "
      "randomisation, and the resulting pair list is saved "
      "(`results/round7/matched_pairs.tsv`). Group A genes are processed in ascending order of "
      "their expression percentile and each takes the nearest unused Group B gene within a "
      "caliper of 0.02; ties are broken by gene name, so the pairing is the same on every run. "
      "Group A cannot be matched completely, because it holds 476 genes above the 90th expression "
      "percentile against 176 in Group B: 1,632 of the 2,016 find a partner and the 384 that do "
      "not sit at a median expression percentile of 0.949. The matched Group A space is the 1,632 "
      "Group A genes of that pair list and the matched Group B space its 1,632 partners, so the "
      "two are the same size and balanced to a standardised mean difference of −0.017 against "
      "+1.236 between the unmatched sets. An earlier matched Group B space of 1,618 genes was "
      "built with a randomised processing order whose pair list was not retained; it is "
      "superseded by this one and is included in File S1 for comparison. The matching of Section "
      "2.4, which compares Δρ between the two gene sets, is a separate application of the same "
      "estimator to the genes with complete values in the state pair being compared, and its 1,754 "
      "pairs are not the pairs used here.")
print('Part 4: §2.6 と §4.14 を対応表の数値にそろえた')

rp(R, "The direction of the loss reproduces in the two other gene spaces of Section 2.6; its size "
      "does not, the aggregated values there being 0.600 and 0.626.",
      "The direction of the loss reproduces in the three other gene spaces of Section 2.6; its "
      "size does not, the aggregated values there being 0.600, 0.624 and 0.529.")
print('§2.3 の空間数を 4 にそろえた')
