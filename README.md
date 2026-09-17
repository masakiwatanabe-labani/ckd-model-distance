# Measurement reliability and preprocessing limit what cross-species model comparisons can establish in chronic kidney disease

Analysis code and derived results for:

> **Measurement Reliability and Preprocessing Limit What Cross-Species Model Comparisons Can
> Establish in Chronic Kidney Disease**
> Masaki Watanabe, Takeru Sasaki, Ryuya Nakagawa and Nobuya Sasaki
> Laboratory of Laboratory Animal Science and Medicine, School of Veterinary Medicine,
> Kitasato University, Towada, Aomori, Japan
>
> Submitted to *International Journal of Molecular Sciences*. The DOI will be added here on
> acceptance.

## Two stages, two tags

This repository was first published at an earlier, exploratory stage of the analysis (tag
`v0.1-exploratory`). The analysis reported in the manuscript is a substantially extended
version; the code and results corresponding to the submitted manuscript are tagged
`v1.0-submitted`. Both tags are kept: the exploratory history is part of the record.

The two stages live in different directories.

directory | stage | what it is
--- | --- | ---
`projection/` | submitted | every number, figure and table in the manuscript
`src/`, `results/` | exploratory | the earlier stage, superseded

**Paths cited in the manuscript are relative to `projection/`.** Where the article says
`results/round7/scale_sensitivity.tsv`, the file is `projection/results/round7/scale_sensitivity.tsv`.

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
gene sets | Enrichr: MSigDB_Hallmark_2020, GO_Biological_Process_2023, KEGG_2021_Human, Reactome_Pathways_2024 | `src/00_fetch_refs.py`

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

It prints a count of agreements and disagreements and exits non-zero if anything disagrees. It
makes **440 checks** and they cover four things:

1. **every number reported in the article** against the result table it came from, to a stated
   tolerance;
2. **every number printed inside a figure**, recorded by `make_figures.py` at the moment it is
   drawn into `projection/results/revision1/pathway_reactome/figure_printed_numbers.tsv`, against
   the caption and the body text;
3. **every cell of Tables 1 to 4**, extracted from the manuscript source, against the result
   files and against the body text, including a rule that no cell may concatenate two numbers;
4. **the structure of the manuscript**: that section, figure and table cross-references resolve,
   that the reference list is continuous and fully cited, that no supplementary item is cited
   without a caption, that no sentence is duplicated and that no working note remains.

Quantities that appear in more than one section are registered so that correcting one and
leaving the other stale fails the script.

## Licence

Code is under the **MIT** licence (`LICENSE`). Derived data — everything under
`projection/results/`, `results/` and the gene-space lists — are under **CC BY 4.0**
(`LICENSE-DATA`). Primary data carry the terms of their own repositories.
