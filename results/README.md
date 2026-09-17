# This directory is the exploratory stage, not the article

These files are the output of the first stage of this project, the state the repository was in at
tag `v0.1-exploratory`. **They are not what the article reports.**

## If you came here from the manuscript

The article cites result files as `results/...`. Those paths are relative to `projection/`, not to
the repository root. So:

article says | the file is
--- | ---
`results/round7/scale_sensitivity.tsv` | `projection/results/round7/scale_sensitivity.tsv`
`results/reliability/reliability.tsv` | `projection/results/reliability/reliability.tsv`
`results/check1/supplementary_sensitivity.tsv` | `projection/results/check1/supplementary_sensitivity.tsv`

and so on for every path the article names. **Everything the article reports is under
`projection/results/`.** Nothing in this directory is cited by the article.

## What is here, and why it is kept

The exploratory stage fixed a different method, and two of its files carry numbers that look like
quantities the article also reports but are not the same quantity:

- `summary.md` gives a module-level cross-species Spearman correlation of 0.92. The article treats
  agreement measured after aggregation as a consequence of the preprocessing rather than as
  evidence that the species agree, and reports the corresponding analysis differently.
- `fig3_mantel_values.csv` gives a Mantel correlation of 0.553 between pairwise distance and
  elapsed time over 14 states. The article reports the recomputation over 12 states in a fixed gene
  space.

They are kept because the exploratory history is part of the record, not because they supersede or
are superseded by anything in the article. See the root `README.md` for the relation between the
two stages.
