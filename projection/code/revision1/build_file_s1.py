"""Part 6. File S1 を組み立てる。

第三者が数値を再現できることを目的とする。含めるもの:
  data/      Δ 行列、対照群発現量、遺伝子空間のリスト、Group A の検出判定と出典、
             オルソログ表と遺伝子セット GMT
  code/      解析・シミュレーション・作図・検証のスクリプト
  results/   本文の数値の裏づけになる結果表
  manuscript/ 図（PDF と 150 dpi PNG）
  README.md  どのスクリプトがどの表を作り、どの表がどの数値を支えるか
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
DEST = Path.home() / "Desktop" / "CKD_xspecies_round7" / "File_S1"

DATA = [
    (HERE / "delta_matrix.tsv", "data/delta_matrix.tsv"),
    (HERE / "delta_matrix_aux.tsv", "data/delta_matrix_aux.tsv"),
    (HERE / "control_log2cpm.tsv", "data/control_log2cpm.tsv"),
    (HERE / "delta_se_moderated.tsv", "data/delta_se_moderated.tsv"),
    (HERE / "delta_se_splithalf.tsv", "data/delta_se_splithalf.tsv"),
    (HERE / "gene_provenance.tsv", "data/gene_provenance.tsv"),
    (HERE / "time_map.tsv", "data/time_map.tsv"),
    (HERE / "groupA_intersection.txt", "data/genespace_groupA_intersection.txt"),
    (HERE / "groupA_union.txt", "data/genespace_groupA_union.txt"),
    (HERE / "aa_ortholog_all.txt", "data/genespace_all_orthologues.txt"),
    (HERE / "aa_groupA_matched.txt", "data/genespace_matched_groupA.txt"),
    (HERE / "aa_groupB_matched2.txt", "data/genespace_matched_groupB.txt"),
    (HERE / "aa_groupB_matched.txt", "data/genespace_matched_groupB_superseded.txt"),
    (HERE / "ensembl_ortholog_genes.txt", "data/ensembl_ortholog_genes.txt"),
    (HERE / "results" / "round7" / "groupA_detection_calls.tsv",
     "data/groupA_detection_calls.tsv"),
    (HERE / "results" / "round7" / "groupA_detection_sources.tsv",
     "data/groupA_detection_sources.tsv"),
    (HERE / "results" / "round7" / "matched_pairs.tsv", "data/matched_gene_pairs.tsv"),
]
REF_FILES = ["orthologs_cat2mouse.tsv", "orthologs_cat2human.tsv", "orthologs_mouse2human.tsv",
             "geneset_gobp.gmt", "geneset_kegg.gmt", "geneset_hallmark.gmt",
             "geneset_reactome.gmt"]
CODE_TOP = ["build_delta_matrix.py", "alpha_decomp.py", "cos_matrix.py", "reliability.py",
            "ceiling_validation.py", "mantel_s1.py", "auc_permutation.py",
            "pathway_separation.py", "pathway_analysis.py", "pathway_cos.py",
            "check1_geneset.py", "pair_uncertainty.py", "human_checks.py",
            "human_extrapolation.py", "human_mapping_matched.py", "build_precision.py",
            "protein_decomp.py", "protein_ceiling_check.py", "aa_summary.py",
            "make_figures.py", "verify_numbers.py"]
RESULT_DIRS = ["alpha_groupA", "control_analysis", "reliability", "mantel_s1", "pathway",
               "check1", "groupAB_control_log2cpm", "pair_uncertainty", "revision1",
               "round5", "round7", "aa_genespace"]


def copy(src: Path, rel: str) -> bool:
    if not src.exists():
        return False
    dst = DEST / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def main() -> int:
    if DEST.exists():
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)
    missing = []
    for src, rel in DATA:
        if not copy(src, rel):
            missing.append(str(src.name))
    for fn in REF_FILES:
        if not copy(ROOT / "data" / "ref" / fn, f"data/ref/{fn}"):
            missing.append(fn)
    for fn in CODE_TOP:
        if not copy(HERE / fn, f"code/{fn}"):
            missing.append(fn)
    for p in sorted((HERE / "code" / "revision1").glob("*.py")):
        copy(p, f"code/revision1/{p.name}")
    for p in sorted((HERE / "code" / "revision1").glob("*.sh")):
        copy(p, f"code/revision1/{p.name}")
    copy(ROOT / "src" / "lib_figure.py", "code/lib_figure.py")
    copy(ROOT / "src" / "lib_stats.py", "code/lib_stats.py")
    copy(ROOT / "src" / "lib_ebayes.py", "code/lib_ebayes.py")
    n_res = 0
    for d in RESULT_DIRS:
        src = HERE / "results" / d
        if not src.exists():
            missing.append(f"results/{d}")
            continue
        for p in sorted(src.rglob("*")):
            if p.is_file() and p.suffix in (".tsv", ".md", ".txt", ".log"):
                copy(p, f"results/{d}/{p.relative_to(src)}")
                n_res += 1
    n_fig = 0
    for p in sorted((HERE / "manuscript" / "figures").glob("*")):
        if p.suffix in (".pdf", ".png"):
            copy(p, f"manuscript/figures/{p.name}")
            n_fig += 1

    files = sorted(p for p in DEST.rglob("*") if p.is_file())
    size = sum(p.stat().st_size for p in files)
    readme = README.format(n_files=len(files), mb=size / 1e6, n_res=n_res, n_fig=n_fig)
    (DEST / "README.md").write_text(readme)
    print(f"File S1: {len(files) + 1} files, {size / 1e6:.1f} MB")
    if missing:
        print("見つからなかったもの:", missing)
    return 0


README = """# File S1

Everything needed to re-derive every number reported in the article. {n_files} files,
{mb:.1f} MB.

## Layout

directory	contents
`data/`	the 16-state Δ matrix, the control-group expression table, the per-gene standard errors, the gene-space lists, the per-gene proteome detection calls behind Group A, the matched gene pairs, and the orthologue tables and gene-set GMT files under `data/ref/`
`code/`	the analysis scripts, with the round-by-round revision scripts under `code/revision1/`
`results/`	{n_res} result tables, one directory per analysis; every number in the article comes from one of these
`manuscript/figures/`	{n_fig} figure files, vector PDF and 150 dpi PNG

## How to check a number

Run `python code/verify_numbers.py` from a copy of the repository. It re-derives every reported
value from the tables in `results/`, compares each with the wording in the manuscript, and prints
the count of agreements and disagreements. It is the same script the authors run.

## Which script produces which table

script	writes	supports
`build_delta_matrix.py`	`data/delta_matrix.tsv`, the gene-space lists, `data/groupA_detection_calls.tsv`	Sections 4.1 to 4.6
`alpha_decomp.py`	`results/alpha_groupA/projection.tsv`	Section 2.1, Figure 1, Table 4
`cos_matrix.py`	`results/control_analysis/`	Section 2.2, Table 1, Figure 2A
`reliability.py`, `ceiling_validation.py`	`results/reliability/`	Section 2.2, Table 1, Table S1, Table S3, Figure 2B,C
`pair_uncertainty.py`, `code/revision1/pair_uncertainty_both_ceilings.py`	`results/pair_uncertainty/`, `results/revision1/pair_uncertainty_both_ceilings.tsv`	Section 2.2, Figure 2D
`mantel_s1.py`, `code/revision1/mantel_subsets.py`	`results/mantel_s1/`, `results/round5/mantel_subsets.tsv`	Section 2.5, Table 3, Figure 6
`code/revision1/pathway_cos.py`, `pathway_separation.py`	`results/revision1/pathway_reactome/`	Section 2.3, Figures 3 and 4
`check1_geneset.py`	`results/check1/`, `results/groupAB_*/`	Section 2.4, Table 2, Figure 5
`code/revision1/build_matched_spaces.py`	`results/round7/matched_pairs.tsv`, the matched gene-space lists	Section 2.6, Section 4.14, Table 4
`code/revision1/ceiling_decomposition_sim.py`	`results/round5/ceiling_decomposition_sim.tsv`	Section 4.10, Figure S1
`code/revision1/interval_coverage_sim.py`, `interval_null_rate.py`	`results/round5/`, `results/round6/`	Section 4.11
`code/revision1/reversal_animal_bootstrap.py`	`results/round5/reversal_animal_bootstrap.tsv`	Section 2.1
`code/revision1/centering_auc_animal_bootstrap.py`	`results/round7/centering_auc_animal_bootstrap.tsv`	Section 2.3
`code/revision1/scale_sensitivity.py`	`results/round7/scale_sensitivity.tsv`	Section 4.4, Table S5
`make_figures.py`, `code/revision1/figure_s1_decomposition.py`	`manuscript/figures/`	all figures
`verify_numbers.py`	nothing; it checks	every reported number

## Software

Python 3.9.6 with numpy 2.0.2, pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6, matplotlib 3.9.4
and gseapy 1.3.1. The random seed is fixed at 20260826 except where a script states its own.

## Data not included here

The four public datasets are at GEO under GSE303653, GSE299326, GSE98622 and GSE104954, and the
feline proteome at ProteomeXchange under PXD066590. The Pod-TRECK proteome is author-held data
from an accepted study and is supplied on request; the tables derived from it, including the
detection calls used to define Group A, are in `data/`.
"""


if __name__ == "__main__":
    sys.exit(main())
