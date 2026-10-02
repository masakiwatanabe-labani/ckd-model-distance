# Mean-centering alters pathway-based transcriptomic comparisons of feline and mouse
kidney disease

Analysis code and derived results for:

> **Mean-Centering Alters Pathway-Based Transcriptomic Comparisons of Feline and Mouse
> Kidney Disease**
> Masaki Watanabe, Takeru Sasaki, Ryuya Nakagawa and Nobuya Sasaki
> Laboratory of Laboratory Animal Science and Medicine, School of Veterinary Medicine,
> Kitasato University, Towada, Aomori, Japan
>
> Under revision. Citation details will be added on acceptance.

## Two stages, two tags

This repository was first published at an earlier, exploratory stage of the analysis (tag
`v0.1-exploratory`). The analysis reported in the submitted manuscript is a substantially
extended version, tagged `v1.0-submitted`. The analysis added while answering the reviewers
is tagged `revision-1`. All three tags are kept: the earlier history is part of the record,
and nothing in it has been rewritten.

What `revision-1` adds, in `projection/code/revision1/` and `projection/results/roundR/`:

what | where
--- | ---
the 75 cross-dataset pair subset | `results/roundR/` (pair subsets), `centering_auc_animal_bootstrap.py`
four definitions of the response vector Δ | `delta_definitions.py` → `results/roundR/delta_definitions*.tsv`
five pathway scoring methods, run in R | `pathway_scoring_methods.py`, `export_for_r.py`, `reference_impl_all.R`
the same five checked against GSEApy | `reference_impl_cells.py`, `reference_impl_compare.py`
two random pathway-set designs | `random_set_designs.py` → `results/roundR/random_set_designs*.tsv`
composition and redundancy of the pathway set | `pathway_set_design.py`
the proteomic detection calls behind Group A | `af_detection_calls.py`, `apply_AV1.py` → `results/roundR/detection_rule.tsv`
duplicate-symbol sensitivity | `duplicate_symbol_sensitivity.py`
the benchmark simulation, summarised by design | `be_simulation_summary.py` → `results/roundR/simulation_exceedance.tsv`
the coordinates actually plotted in Figures S1–S5 | `bd_figures.py` → `manuscript_C/figure_data/*.tsv`
checks of the figures against the text | `bd_check_figure_text.py`
the bibliography added in revision, from Crossref | `bi_add_references.py` → `results/roundR/crossref_records.json`

The two stages live in different directories.

directory | stage | what it is
--- | --- | ---
`projection/` | submitted | every number, figure and table in the manuscript
`src/`, `results/` | exploratory | the earlier stage, superseded

**Paths cited in the manuscript are relative to `projection/`.** Where the article says
`results/round7/scale_sensitivity.tsv`, the file is
`projection/results/round7/scale_sensitivity.tsv`.

The files under `results/` are the output of the earlier stage. Two of them name quantities
that also appear, with different values, in the manuscript, because the method was not yet
fixed when they were produced: `summary.md` reports a module-level cross-species Spearman
correlation of 0.92, and `fig3_mantel_values.csv` reports a Mantel correlation of 0.553 between
pairwise distance and elapsed time over 14 states. The manuscript reports the corresponding
recomputation over 12 states in a fixed gene space, and treats agreement measured after
aggregation as a consequence of the preprocessing rather than as evidence that the species
agree. Neither file is cited in the manuscript.

## What the analysis does

It takes one feline spontaneous chronic kidney disease state as a fixed reference axis and
compares sixteen feline, podocyte-injury and ischaemia-reperfusion states against it, asking
what a similarity score of that kind can settle. The three results are that a projection onto a
reference axis ranks direction and amplitude together and can invert the ranking; that how far
an observed similarity falls below a reliability-derived benchmark is interpretable only where
the reference state is itself well measured; and that averaging genes within pathways separates
cat-mouse pairs of states from same-species pairs at AUC 0.868 or 0.511 according to whether
each state's mean change across genes is removed first.

## Data

**No primary data are in this repository.** See [`data/README.md`](data/README.md) for where to
get each input, what to name it, where to put it and its MD5 checksum.

source | accession | how
--- | --- | ---
feline spontaneous CKD | GSE303653 | supplementary workbook, manual download
feline proteome | PXD066590 | ProteomeXchange
Pod-TRECK podocyte injury, transcriptome | GSE299326 | automatic
murine ischaemia-reperfusion | GSE98622 | automatic
human tubulointerstitial | GSE104954 | automatic
orthologues | Ensembl BioMart | `src/00_fetch_refs.py`
gene sets | Enrichr: MSigDB_Hallmark_2020, GO_Biological_Process_2023, KEGG_2021_Human,
Reactome_Pathways_2024 | `src/00_fetch_refs.py`

Two inputs are deliberately absent and cannot be regenerated from this repository.

- **The Pod-TRECK proteome.** It is author-held data from a parallel study that is accepted but
  not yet published, so its abundance tables are supplied on request rather than deposited here.
  What the manuscript needs in order for the reader to check the definition of Group A is the
  per-gene detection call, and that is published, as
  `projection/results/round7/groupA_detection_calls.tsv`: for every gene, how many samples of
  each proteome it was detected in and whether it therefore enters Group A.
- **Gene-set collections and orthologue tables.** Retrieved from Enrichr and Ensembl BioMart,
  whose redistribution terms are not ours to grant. `src/00_fetch_refs.py` downloads them and
  warns if a collection has changed size since the analysis was run.

## Requirements

Python 3.9 or newer.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` lists direct dependencies. `requirements-lock.txt` is the exact environment
(`pip freeze`) used to produce the published numbers. The numbers here were produced with
Python 3.9.6, numpy 2.0.2, pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6, matplotlib 3.9.4 and
gseapy 1.3.1.

## Running the submitted analysis

```bash
cd projection
V=../.venv/bin/python

$V build_delta_matrix.py       # the 16-state Δ matrix and the gene-space lists   ~2 min
$V build_precision.py          # per-gene standard errors                         ~3 min
$V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
     --genes groupA_intersection.txt --time time_map.tsv --out results/alpha_groupA
$V cos_matrix.py               # all 120 pairs, shared-control split              ~2 min
$V reliability.py              # split-half reliability, 500 splits per state     ~8 min
$V ceiling_validation.py       # the benchmark in both forms                      ~3 min
$V mantel_s1.py                # Mantel and partial Mantel                        ~3 min
$V pathway_separation.py       # pathway-level and aggregated separation          ~6 min
$V code/revision1/pathway_centering_steps.py   # the four centring/aggregation conditions
$V code/revision1/pathway_leads.py             # which state leads how many pathways
$V code/revision1/build_provenance.py          # results/provenance.tsv
$V pair_uncertainty.py         # animal-level resampling                          ~10 min
$V verify_numbers.py           # the check described below                        ~1 min
$V make_figures.py             # Figures 1 to 6                                   ~1 min
```

About 40 minutes end to end on a laptop, once the reference tables are in place. The
simulations behind Section 4.10 and Section 4.11 are separate and slower; they are in
`projection/code/revision1/` and each writes into `projection/results/round5/` or `round7/`.

The earlier stage still runs as it did: `make refs` then `make analysis` from the repository
root, about 25 minutes in total.

## Checking the numbers

```bash
cd projection && ../.venv/bin/python verify_numbers.py
```

It prints how many checks ran, how many agreed, how many disagreed and which groups of checks
were skipped for want of material this repository does not publish, and exits non-zero if
anything disagrees. With the manuscript sources and the reference gene sets in place it makes
**529 checks**; a bare clone carries neither, and there it makes **206 checks** and names the
three groups it skipped; after `python src/00_fetch_refs.py` it makes **219**. The checks the
revision added all compare the manuscript or the submitted files against the result files, so
they are among the ones a clone skips. The checks cover seven things:

1. **every number reported in the article** against the result table it came from, to a stated
   tolerance;
2. **every number printed inside a figure**, recorded by `make_figures.py` at the moment it is
   drawn into `projection/results/revision1/pathway_reactome/figure_printed_numbers.tsv`, against
   the caption and the body text;
3. **every cell of Tables 1 to 4**, extracted from the manuscript source, against the result
   files and against the body text, including a rule that no cell may concatenate two numbers;
4. **the structure of the manuscript**: that section, figure and table cross-references resolve,
   that the reference list is continuous and fully cited, that no supplementary item is cited
   without a caption, that no sentence is duplicated and that no working note remains;
5. **that the result files are regenerated by the code**, for the tables behind Section 2.3. The
   four aggregation conditions and the per-state leading shares are recomputed inside the same
   run, from the Δ matrix and the gene-set definitions, and compared cell by cell with the stored
   tables. Altering a single cell of
`projection/results/revision1/pathway_reactome/cos_before_after.tsv`
   fails the script;
6. **that every file under `projection/results/` has a generating script**, read from
   `projection/results/provenance.tsv` (see below);
7. **every summary value the article states about a supplementary figure**, recomputed from the
   coordinates that figure actually plots. `bd_figures.py` records those coordinates straight
   off the matplotlib artists into `projection/manuscript_C/figure_data/`, and
   `bd_check_figure_text.py` recomputes the medians, interquartile ranges, AUCs, exceedance
   rates and counts from them. A figure redrawn from different numbers fails the script.

Quantities that appear in more than one section are registered so that correcting one and
leaving the other stale fails the script.

## Where each result file comes from

`projection/results/provenance.tsv` lists every file under `projection/results/` with the script
that writes it and how that was established (an exact filename in a write call, or an f-string
pattern for names built at run time). It is generated:

```bash
cd projection && ../.venv/bin/python code/revision1/build_provenance.py
```

Two files are listed as `intermediate, not used by the article`: they were produced during the
work, no script here writes them, and no reported number depends on them. They are kept under
`projection/results/intermediate/`, with a README saying what they hold and which regenerable file
supersedes each. Every other file has a generating script, and `verify_numbers.py` fails if that
stops being true.

## Which file each reported quantity comes from

The article does not cite result files by path. This table is the mapping instead: for every
quantity `verify_numbers.py` checks, the file it is read from and the script that writes that
file. It is generated, not maintained by hand —

```bash
cd projection && ../.venv/bin/python verify_numbers.py --map \
  && ../.venv/bin/python code/revision1/build_readme_map.py
```

writes `projection/results/round9/quantity_to_file.tsv` (one row per quantity) and
`readme_file_map.tsv` / `readme_file_map.md` (one row per file), of which the table below is a
copy. The generating-script column is read from `provenance.tsv`, so the two tables cannot drift
apart. Quantity labels are the internal ones the verification script uses; the section numbers in
them refer to `NUMBERS.md`, the working record, and not to the sections of the article.

| Result file | Produced by | Quantities it supplies |
|---|---|---|
| `projection/results/aa_genespace/genespace_comparison.tsv` | `projection/aa_summary.py` | §15 遺伝子数 orthologue; §15 遺伝子数 matched B; §15 Group A intersection（主解析） 猫髄質 alpha; §15 Group A intersection（主解析） cos x log_time; §15 全 1:1 orthologue 猫髄質 alpha; §15 全 1:1 orthologue cos x log_time; §15 発現量マッチ Group B 猫髄質 alpha; §15 発現量マッチ Group B cos x log_time |
| `projection/results/alpha_groupA/projection.tsv` | `projection/alpha_decomp.py` | §6 cat_med_CKD12 alpha; §6 cat_med_CKD12 cos; §6 cat_med_CKD12 nrm; §6 mouse_2W alpha; §9 complete case 遺伝子数; §2.7 Group A intersection（主解析） 髄質 nrm - 1/cos |
| `projection/results/alpha_groupA/time_association.tsv` | `projection/alpha_decomp.py` | §12 cos x log_time（状態ごと）; §12 nrm x log_time（状態ごと） |
| `projection/results/alpha_groupB_matched/projection.tsv` | `projection/alpha_decomp.py` | §2.7 発現量マッチ Group B 髄質 nrm - 1/cos |
| `projection/results/alpha_ortholog_all/projection.tsv` | `projection/alpha_decomp.py` | §2.7 全 1:1 orthologue 髄質 nrm - 1/cos |
| `projection/results/alpha_union/projection.tsv` | `projection/alpha_decomp.py` | §9 union complete case; §9 GroupA intersection n; §9 GroupA union n |
| `projection/results/check1/check1_nn_matched.tsv` | `projection/check1_geneset.py / projection/verify_numbers.py` | §1 Δρ added−core（対照） |
| `projection/results/check1/supplementary_sensitivity.tsv` | `projection/check1b_recipe.py / projection/verify_numbers.py` | §7 S1 rho; §7 S7 rho |
| `projection/results/control_analysis/auc_permutation.tsv` | `projection/auc_permutation.py` | §4 AUC 観測; §4 AUC p |
| `projection/results/control_analysis/cos_distributions.tsv` | `projection/cos_matrix.py / projection/verify_numbers.py` | §3 共有対照あり 中央; §3 共有対照なし 中央 |
| `projection/results/groupAB_control_log2cpm/check1_by_decile.tsv` | `projection/check1_geneset.py` | §10 第1十分位 GroupA; §10 第1十分位 GroupB |
| `projection/results/groupAB_control_log2cpm/check1_nn_matched.tsv` | `projection/check1_geneset.py / projection/verify_numbers.py` | §1 Δρ 発現量マッチ後; §1   CI lo; §1   CI hi |
| `projection/results/groupAB_control_log2cpm/check1_summary.tsv` | `projection/check1_geneset.py / projection/verify_numbers.py` | §1 ρ GroupA 素; §1 ρ GroupB 素; §1 Δρ 素 |
| `projection/results/groupAB_control_log2cpm/check1_unmatched_delta.tsv` | `projection/check1_geneset.py` | §1 Δρ 素 点推定; §1 Δρ 素 CI lo; §1 Δρ 素 CI hi |
| `projection/results/groupAB_delta_se_moderated/check1_nn_matched.tsv` | `projection/check1_geneset.py / projection/verify_numbers.py` | §1 Δρ delta_se_moderated マッチ後 |
| `projection/results/groupAB_delta_se_splithalf/check1_nn_matched.tsv` | `projection/check1_geneset.py / projection/verify_numbers.py` | §1 Δρ delta_se_splithalf マッチ後 |
| `projection/results/human_exploratory/cat_vs_mouse_to_human.tsv` | `projection/human_extrapolation.py` | §11 AUC（マッチなし）; §11 p（マッチなし） |
| `projection/results/human_exploratory/human_cos_summary.tsv` | `projection/human_extrapolation.py` | §11 cat x human 中央; §11 mouse x human 中央 |
| `projection/results/human_exploratory/mapping_matched_test.tsv` | `projection/human_mapping_matched.py` | §11 caliper2 cat x human; §11 caliper2 mouse x human; §11 caliper2 差; §11 caliper2 AUC; §11 caliper2 pid_cat; §11 caliper2 pid_mouse |
| `projection/results/mantel_s1/distance_distributions.tsv` | `projection/mantel_s1.py` | §12 D_cos 種間 中央; §12 D_cos 種間 最小 |
| `projection/results/mantel_s1/mantel_s1.tsv` | `projection/mantel_s1.py` | §12 Mantel 時間 D_rho; §12 Mantel 時間 D_cos; §12 Mantel 時間 D_nrm; §12 Mantel 入口 D_rho(12) |
| `projection/results/pair_uncertainty/pair_uncertainty.tsv` | `projection/pair_uncertainty.py` | §16 異種ペア数; §16 区間が0を除外する異種ペア; §16 信頼性の高い猫状態を含むペア; §16 同 区間が0を除外; §16 同 gap 中央; §16 低信頼の猫状態を含むペアで0を除外; §16 ジャックナイフで gap<=0; §16 同種(対照非共有)で0を除外 |
| `projection/results/pathway/aggregation_steps.tsv` | `projection/code/revision1/pathway_centering_steps.py` | §2.5 1. cos of Δ as analysed 種内; §2.5 1. cos of Δ as analysed 種間; §2.5 2. Δ centred per state (= Pearson r of Δ) 種内; §2.5 2. Δ centred per state (= Pearson r of Δ) 種間; §2.5 3. pathway means of centred Δ 種内; §2.5 3. pathway means of centred Δ 種間; §2.5 4. pathway means of uncentred Δ 種内; §2.5 4. pathway means of uncentred Δ 種間 |
| `projection/results/pathway/cos_before_after.tsv` | `projection/code/revision1/pathway_centering_steps.py` | §2.5 AUC Δ そのまま; §2.5 AUC 中心化のみ; §2.5 AUC 中心化+経路平均; §2.5 AUC 中心化なし経路平均 |
| `projection/results/pathway/pathway_cos.tsv` | `projection/code/revision1/pathway_cos.py / projection/pathway_cos.py` | §13 経路数; §13 セル数; §13 経路サイズ中央; §2.5 マウス×経路セル数; §2.5 n>=50 のセル数; §2.5 cos>=0.9 のセル数（条件なし）; §2.5 条件を満たし cos>=0.9 のセル数 |
| `projection/results/pathway/pathway_separation_paired.tsv` | `projection/code/revision1/pathway_separation.py / projection/pathway_separation.py / projection/verify_numbers.py` | §13 経路単位 AUC 中央; §13 経路単位 AUC Q1; §13 経路単位 AUC Q3; §13 AUC>=0.90 の経路数; §13 AUC<=0.60 の経路数 |
| `projection/results/pathway/separation_summary.tsv` | `projection/code/revision1/pathway_separation.py / projection/pathway_separation.py` | §13 集約 AUC; §13 集約 AUC（コホート除外）; §13 制限後のペア数 |
| `projection/results/protein/projection_protein.tsv` | `projection/protein_decomp.py` | §14 蛋白数; §14 cos 猫皮質CKD1/2; §14 cos 猫髄質CKD1/2; §14 alpha D14; §14 alpha D21; §14 alpha 猫髄質CKD1/2 |
| `projection/results/protein/protein_vs_rna.tsv` | `projection/protein_decomp.py` | §14 protein vs RNA cos 順位一致; §14 cos 中央 protein; §14 cos 中央 RNA |
| `projection/results/reliability/auc_permutation_cos_corrected.tsv` | `projection/auc_permutation.py` | §4 AUC 補正後; §4 AUC 補正後 p |
| `projection/results/reliability/ceiling_comparison.tsv` | `projection/ceiling_validation.py` | §10 天井 同一DS raw; §10 天井 同一DS SB; §10 天井 同種異DS raw; §10 マッチ成立数; §10 マッチ後 発現量順位 A; §10 マッチ後 発現量順位 B |
| `projection/results/reliability/pairs_ceilings.tsv` | `projection/ceiling_validation.py` | §3 異種 cos 中央; §3 異種 cos 最大; §3 同種異DS cos 中央; §3 天井raw 異種 最小; §3 天井SB 異種 最小; §3 天井までの最小余裕; §3 系統誤差の下限; Table 1 ペア数 Within dataset; Table 1 cos 中央値 Within dataset; Table 1 cos 最大 Within dataset; Table 1 ペア数 Same species, different dataset; Table 1 cos 中央値 Same species, different dataset; Table 1 cos 最大 Same species, different dataset; Table 1 ペア数 Cross-species (cat vs mouse); Table 1 cos 中央値 Cross-species (cat vs mouse); Table 1 cos 最大 Cross-species (cat vs mouse) |
| `projection/results/reliability/r_curve.tsv` | `projection/ceiling_validation.py` | §10 AUC 帰無中央値 |
| `projection/results/reliability/reliability.tsv` | `projection/ceiling_validation.py / projection/code/revision1/ceiling_decomposition_sim.py / projection/reliability.py / projection/verify_numbers.py` | §10 信頼性 cat_CKD34; §10 信頼性 cat_med_CKD12; §10 信頼性 cat_CKD12; §10 信頼性 cat_med_CKD34; SB 基準 (v) 真の応答が異なるときの超過（中心化）; SB 基準 (v) 真の応答が異なるときの超過（非中心化）; 無補正基準 (iii) の平均超過; 標本数の寄与 (ii)-(i); 非中心化の寄与 (iii)-(ii); SB 基準での非中心化の寄与 (v)-(iv); SB 基準 真の応答が同一のときの超過 |
| `projection/results/reliability/sb_validation.tsv` | `projection/ceiling_validation.py / projection/verify_numbers.py` | §10 SB検証 最小; §10 SB検証 最大; §10 SB検証 中央 |
| `projection/results/revision1/cosine_attenuation_sim.tsv` | `projection/code/revision1/cosine_attenuation_sim.py` | §4.10 真の cos 0.7 での最大超過; §4.10 実データの平均シフト中央 |
| `projection/results/revision1/feline_filter_sensitivity.tsv` | `projection/code/revision1/feline_filter_check.py` | §4.3 フィルタで上がる r_half 猫皮質CKD1/2 |
| `projection/results/revision1/gap_by_group.tsv` | `projection/code/revision1/gap_by_group.py` | §2.2 同種 gap 中央; §2.2 高信頼 異種 gap 中央; §2.2 低信頼 異種 gap 中央; §2.2 同種で区間が0を除外; §2.2 SB 同種 gap 中央; §2.2 SB 高信頼 gap 中央; §2.2 SB 低信頼 gap 中央 |
| `projection/results/revision1/gap_group_permutation.tsv` | `projection/code/revision1/gap_by_group.py` | §2.2 置換 観測差（無補正）; §2.2 置換 p（無補正）; §2.2 置換 観測差（SB）; §2.2 置換 p（SB） |
| `projection/results/revision1/genespace_with_reactome.tsv` | `projection/code/revision1/genespace_with_reactome.py` | §2.6 経路 AUC Group A; §2.6 経路 AUC 全ortholog; §2.6 経路 AUC マッチB; §2.6 集約 AUC 全ortholog; §2.6 集約 AUC マッチB; Table 4 集約後 AUC Group A; Table 4 経路単位 AUC 中央値 Group A; Table 4 集約後 AUC all 1:1 orthologues; Table 4 経路単位 AUC 中央値 all 1:1 orthologues; Table 4 集約後 AUC matched Group B; Table 4 経路単位 AUC 中央値 matched Group B; Table 4 集約後 AUC matched Group A; Table 4 経路単位 AUC 中央値 matched Group A |
| `projection/results/revision1/groupA_expression_by_species.tsv` | `projection/code/revision1/feline_filter_check.py` | §4.3 Group A 発現分位 猫; §4.3 Group A 発現分位 Pod-TRECK; §4.3 Group A 発現分位 IRI |
| `projection/results/revision1/leads_by_state_breakdown.tsv` | `projection/code/revision1/pathway_leads.py` | 再生成 IRI 2 h が首位の経路数; 再生成 うち Reactome |
| `projection/results/revision1/pathway_reactome/cos_before_after.tsv` | `projection/code/revision1/pathway_centering_steps.py` | §2.3(422) cos_raw 種内; §2.3(422) cos_raw 種間; §2.3(422) cos_centred 種内; §2.3(422) cos_centred 種間; §2.3(422) cos_agg_centred 種内; §2.3(422) cos_agg_centred 種間; §2.3(422) cos_agg_raw 種内; §2.3(422) cos_agg_raw 種間; §2.3(422) AUC 中心化+経路平均; §2.3(422) AUC 中心化なし経路平均; 再生成 セット数 3 コレクション; 再生成 セット数 Reactome 込み; 再生成 §2.3 AUC 遺伝子単位; 再生成 §2.3 AUC 中心化のみ; 再生成 §2.3 AUC 中心化+経路平均; 再生成 §2.3 AUC 中心化なし経路平均; 再生成 §2.5 AUC 中心化なし経路平均（193） |
| `projection/results/revision1/pathway_reactome/pathway_separation_paired.tsv` | `projection/code/revision1/pathway_separation.py / projection/pathway_separation.py / projection/verify_numbers.py` | §2.3 経路数（4コレクション）; §2.3 経路単位 AUC 中央; §2.3 同 四分位下; §2.3 同 四分位上; §2.3 AUC>=0.90 の経路; §2.3 AUC<=0.60 の経路; §2.3 AUC 0.90 以上の経路数; §2.3 AUC 0.60 以下の経路数; §2.3 全パネル以上の経路の割合 |
| `projection/results/round5/interval_coverage_sim.tsv` | `projection/code/revision1/interval_coverage_sim.py` | 被覆率 無補正 平均; 被覆率 SB 平均 |
| `projection/results/round5/interval_null_rate_SB.tsv` | `projection/code/revision1/interval_null_rate.py` | 無補正基準での偽陽性率 最大 |
## What this repository does not contain, and why

thing | why | how to obtain it
--- | --- | ---
the manuscript text, in any form | it is the article; the version of record is the journal's | the published article
the reviewers' comments and the response letter | peer review is confidential | not available
the two third-party primary files (Li et al. supplementary workbook; the GSE98622 supplementary file) | not ours to redistribute | `data/README.md` names each file, its origin, its MD5 and its SHA-256; `src/00b_fetch_external.py` fetches the GEO one
the Pod-TRECK proteome abundance tables | author-held data from a parallel study, and a contract analysis whose redistribution terms are still being checked | the per-gene detection calls that Group A rests on **are** published, as `projection/results/roundR/groupA_detection_calls.tsv`; the abundance tables are supplied on request
the gene-set collections (`*.gmt`) | MSigDB's KEGG collection restricts redistribution | `python src/00_fetch_refs.py` regenerates all four in a few minutes
per-sample expression and abundance matrices | effectively the primary data in another container | as above, per dataset

Scripts whose content is the manuscript text itself (the editing scripts that carry the old and
the new wording) are also kept out, for the same reason as the manuscript.

## Licence

Code is under the **MIT** licence (`LICENSE`). Derived data — everything under
`projection/results/`, `results/` and the gene-space lists — are under **CC BY 4.0**
(`LICENSE-DATA`). Primary data carry the terms of their own repositories.
