# Intermediate files kept for the record

These files were produced during the analysis but no script in this repository writes them, and
no number reported in the article comes from them. They are kept because they are part of the
working record, and moved here so that `results/provenance.tsv` does not have to carry an
unexplained gap.

| File | What it holds | Why it is not needed |
|---|---|---|
| `aggregation_distributions.tsv` | Minimum, quartiles, median and maximum of cos θ before and after aggregation to 193 pathway means, by pair class | The quartiles and medians are in `results/pathway/aggregation_steps.tsv`, which `code/revision1/pathway_centering_steps.py` regenerates. Only the minimum and maximum are unique to this file, and the article reports neither |
| `mapping_strata_cos.tsv` | cat-to-human, mouse-to-human and cat-to-mouse cos θ under three ways of splitting genes by orthologue percent identity | The values the article reports for Section 2.7 are in `results/human_exploratory/human_cos_summary.tsv`, which `human_extrapolation.py` writes and `verify_numbers.py` checks. This file's own strata are not reported |

If either becomes needed again, the computation behind it has to be written; the script that
produced them was never committed.
