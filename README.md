# Cross-species distances in chronic kidney disease

Analysis code for:

> **[Manuscript title]**
> [Author list]
>
> *Citation details will be added on acceptance.*

## What this code does

It measures how far apart chronic kidney disease states are from one another —
across species (cat, mouse, human) and across mouse models — using rank-based
distances between differential-expression profiles projected into a common human
gene space. It then asks what actually structures those distances: the
compartment where disease begins, the model used, or the time elapsed since
injury.

## Data

**No primary data are included in this repository.** Some are third-party
copyright, some belong in the repository that hosts them. Four of the six
inputs download automatically; two need a browser.

See **[`data/README.md`](data/README.md)** for where to get each file, what to
name it, where to put it, and its MD5 checksum.

## Requirements

Python 3.9 or newer. Install with:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` lists direct dependencies. `requirements-lock.txt` is the
exact environment (`pip freeze`) used to produce the published numbers; install
from that if you need bit-for-bit agreement.

## Running it

```bash
# 1. reference tables (Ensembl BioMart, UniProt, MSigDB) and public GEO series
make refs            # ~10 min, network required, writes data/ref/ and data/external/

# 2. obtain the two manual downloads - see data/README.md
#    (feline supplementary workbook, KPMP DataLake_DEPs.txt)

# 3. everything else
make analysis        # ~15 min on a laptop, no network required
```

`make all` runs both stages. Individual steps are plain scripts and can be run
on their own, in numerical order, e.g. `python src/02_de.py`.

`data/ref/` is created on first run of `src/00_fetch_refs.py` and is not stored
in this repository. It includes MSigDB gene set collections; retrieving them
means accepting the MSigDB terms of use, and the KEGG collection in particular
carries redistribution conditions, which is why we do not ship a copy.

Random seeds, thresholds and the number of permutations are all set in
`config/config.yaml`; gene set definitions are in `config/modules.yaml`.

## Output

Everything is written to `results/`. Figures for the paper are in
`results/figures/`, exploratory figures that were not used are in
`results/supplementary/`, and `results/summary.md` accumulates the principal
numbers from each step in the order they were produced (with a machine-readable
`results/summary.json` alongside).

### Where each figure and table comes from

| Paper item | File | Produced by |
|---|---|---|
| Figure 1 | `results/figures/Fig1_design.png` / `.pdf` | `src/22_fig1.py` |
| Figure 2 | `results/figures/Fig2_concordance.png` / `.pdf` | `src/23_fig2.py` |
| Figure 3 | `results/figures/Fig3_distance_structure.png` / `.pdf` | `src/17_fig3.py` |
| Figure 4 | `results/figures/Fig4_distance_distributions.png` / `.pdf` | `src/19_distance_distributions.py` |
| Figure 5 | `results/figures/Fig5_iri_trajectory.png` / `.pdf` | `src/27_fig5.py` |
| Figure 6 | `results/figures/Fig6_module_level.png` / `.pdf` | `src/26_fig6.py` |

Every figure is written twice, as PNG (300 dpi) and as PDF. The PDFs embed
fonts as TrueType (`pdf.fonttype = 42`), so all labels remain live text and can
be selected and edited in Illustrator or Inkscape rather than being converted to
outlines. Style and output format are set in one place, `src/lib_figure.py`.
| Table 1 | `results/table1_datasets.csv` + `table1_footnotes.txt` | `src/21_table1.py` |
| Table 2 | `results/table2_candidates.csv` + `table2_footnotes.txt` | `src/28_table2.py` |

### Principal result files

| Result | File |
|---|---|
| All differential expression contrasts | `results/de_all.csv.gz` |
| Cross-species concordance and within-species benchmarks | `results/crossspecies_concordance.csv`, `within_species_benchmarks.csv` |
| Permutation nulls and their diagnostics | `results/crossspecies_permutation.csv`, `null_diagnostics.csv` |
| Species x onset-compartment decomposition | `results/design_decomposition_rna.csv`, `design_decomposition_prot.csv` |
| Distance distributions by pair class | `results/distance_distributions_16state.csv`, `distance_distribution_summary.csv` |
| Elapsed time versus onset compartment (Mantel) | `results/time_axis_full_pairs_no6mo.csv`, `fig3_mantel_values.csv` |
| Module-level distances and power | `results/module_distance_pairs_*.csv`, `module_power_curve.csv` |
| Candidate validation in human cohorts | `results/candidates_human_validation.csv`, `ercb_candidate_validation.csv` |
| References used in the manuscript | `results/references.bib` |

## Pipeline

| Step | Script | Purpose |
|---|---|---|
| 00 | `00_fetch_refs.py` | Orthologues, synonyms, secretome, gene sets |
| 00b | `00b_fetch_external.py` | GSE98622, GSE79443, GSE104954 |
| 01 | `01_load.py` | Read primary data, resolve gene identifiers |
| 02 | `02_de.py` | Moderated *t* contrasts, leave-one-out stability |
| 03 | `03_crossspecies.py` | Cross-species concordance, benchmarks, permutations |
| 04 | `04_modules.py` | Module scores and coherence tests |
| 05-07 | `05_stage.py`, `06_compartment.py`, `07_candidates.py` | Disease stage, compartment, candidate selection |
| 09-12 | `09_dataset_matrix.py` … `12_model_distance.py` | Distance matrices, human validation, ERCB, model distances |
| 13-16 | `13_time_axis.py` … `16_entry_model_checks.py` | Elapsed time versus onset compartment |
| 18-21 | `18_monotonic_genes.py` … `21_table1.py` | Trajectory classes, distributions, null diagnostics, Table 1 |
| 24-28 | `24_module_level_distance.py` … `28_table2.py` | Module-level distances, power, Table 2 |
| 17, 22-27 | `17_fig3.py`, `22_fig1.py`, `23_fig2.py`, `26_fig6.py`, `27_fig5.py` | Paper figures |
| 08 | `08_figures.py` | Exploratory figures (superseded by Figures 1-6; some labels are not in English) |

## Licence

Code is released under the MIT Licence; see [`LICENSE`](LICENSE). The licence
covers this code only. The primary data are not redistributed here and remain
subject to the terms of their respective sources.
