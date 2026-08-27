"""Table 1: 本解析で用いた全データセット / 疾患状態の一覧。"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("21_table1")
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]

C = ["Species", "Model / cohort", "Onset compartment", "Layer", "Accession", "Platform",
     "Disease group (n)", "Control group (n)", "Days since injury", "Contrast used here"]

R = []


def add(**kw):
    R.append({c: kw.get(c.replace(" ", "_").replace("/", "").replace("(", "")
                        .replace(")", "").replace("-", "_"), "") for c in C})


def row(sp, mc, oc, ly, ac, pf, dg, cg, dsi, con):
    R.append(dict(zip(C, [sp, mc, oc, ly, ac, pf, dg, cg, dsi, con])))


# ---------------- cat ----------------
# n は data/interim/*_grp.parquet の実測値。
CAT_SRC = "Li et al. 2025, Commun Biol; suppl. S3 (doi:10.1038/s42003-025-09164-8)"
for ly, plat, tis, dn, cn, con in [
        ("transcriptome", "RNA-seq (not deposited in GEO)", "cortex", 7, 6, "cat_rna_ctx_late"),
        ("transcriptome", "RNA-seq (not deposited in GEO)", "medulla", 5, 6, "cat_rna_med_late"),
        ("proteome", "LC-MS/MS (row-centred values)", "cortex", 7, 6, "cat_prot_ctx_late"),
        ("proteome", "LC-MS/MS (row-centred values)", "medulla", 6, 6, "cat_prot_med_late")]:
    row("Felis catus", "Natural CKD (IRIS)", "tubular", ly,
        CAT_SRC if con == "cat_rna_ctx_late" else "same as above", plat,
        f"IRIS 3/4 {tis} (n={dn})", f"Control {tis} (n={cn})",
        "not defined (natural onset)", con)

# ---------------- Pod-TRECK ----------------
for d, lab, n in [(5, "5D", 3), (14, "2W", 3), (21, "3W", 3)]:
    row("Mus musculus", "Pod-TRECK (podocyte DT injury)", "glomerular", "transcriptome",
        "GSE299326", "RNA-seq (FPKM table)", f"DT 100 ng, {lab} (n={n})",
        "Non-induced control (n=3)", str(d), f"m_rna_{lab.lower()}")
row("Mus musculus", "Pod-TRECK (podocyte DT injury)", "glomerular", "transcriptome",
    "GSE299326", "RNA-seq (FPKM table)", "DT 250 ng, 2W (n=3)", "Non-induced control (n=3)",
    "14", "m_rna_2w250")
for d, lab in [(14, "Day14"), (21, "Day21")]:
    row("Mus musculus", "Pod-TRECK (podocyte DT injury)", "glomerular", "proteome",
        "PK25055 (in-house; not deposited)", "LC-MS/MS", f"{lab} (n=3)", "Control (n=3)",
        str(d), f"m_prot_d{d}")

# ---------------- IRI (GSE98622) ----------------
IRI = [("IRI2h", 2 / 24, "0.083", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI4h", 4 / 24, "0.167", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI24h", 1, "1", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI48h", 2, "2", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI72h", 3, "3", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI7d", 7, "7", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI14d", 14, "14", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI28d", 28, "28", "SHAM4h + SHAM24h (n=6)", "GPL13112"),
       ("IRI12m", 365, "365", "SHAM12m (n=3)", "GPL13112"),
       ("IRI6mN", 180, "180", "NORM3m + NORM9m + NORM15m (n=9)", "GPL19057")]
for lab, _, days, ctrl, gpl in IRI:
    n = 4 if lab == "IRI6mN" else 3
    plat = ("Illumina HiSeq 2000 (GPL13112)" if gpl == "GPL13112"
            else "Illumina NextSeq 500 (GPL19057)")
    row("Mus musculus", "Unilateral IRI (whole kidney)", "tubular", "transcriptome",
        "GSE98622", plat, f"{lab} (n={n})", ctrl, days, lab)

# ---------------- UUO (GSE79443) ----------------
for lab, days in [("2D", "2"), ("8D", "8")]:
    row("Mus musculus", "UUO (right ureteric ligation)", "tubular", "transcriptome",
        "GSE79443", "Illumina HiSeq 2000 (GPL13112)", f"{lab} obstructed kidney (n=3)",
        "Sham-operated kidneys (n=4)*", days, f"UUO_{lab}")

# ---------------- human ----------------
for comp, nm in [("TI", "tubulointerstitium"), ("G", "glomerulus")]:
    row("Homo sapiens", f"KPMP regional proteomics ({nm})", "glomerular (DKD/HTN)", "proteome",
        "KPMP DataLake_DEPs.txt (doi:10.48698/mg7h-bc51)", "LC-MS/MS on LMD tissue",
        "CKD (pre-computed contrast)", "Healthy reference tissue (HRT)",
        "not defined (natural onset)", f"CKD.vs.HRT.in.{comp}")
for lab, dg, cg, gpl in [
        ("GPL22945: DN vs LD", "Diabetic nephropathy (n=7)", "Living donor (n=18)", "GPL22945"),
        ("GPL22945: RPGN vs LD", "ANCA vasculitis (n=21)", "Living donor (n=18)", "GPL22945"),
        ("GPL24120: HT vs LD+TN", "Hypertensive nephropathy (n=20)",
         "Living donor + tumour nephrectomy (n=5)", "GPL24120")]:
    row("Homo sapiens", "ERCB tubulointerstitium", "varies by diagnosis", "transcriptome",
        "GSE104954", f"Affymetrix, Brainarray ENTREZG CDF ({gpl})", dg, cg,
        "not defined (natural onset)", lab)

t = pd.DataFrame(R)[C]
t.to_csv(RES / "table1_datasets.csv", index=False)
log.info("Table 1: %d 行", len(t))

FOOT = [
    "* UUO (GSE79443): only the obstructed (right) kidney was deposited; contralateral "
    "kidneys of UUO animals are not available. Sham-operated kidneys were therefore used as "
    "controls. Sham comprises left and right kidneys from 2 animals (4 samples); "
    "left-right difference was negligible (Spearman rho = 0.992, median |difference| = "
    "0.19 log2), so all 4 were pooled, but they represent 2 animals (pseudo-replication).",
    "GSE98622 is split across two sequencing platforms. Contrasts were always built within "
    "a platform; controls were never pooled across platforms.",
    "Feline proteome values are centred per protein (row mean ~ 0), so only rank-based "
    "statistics were used for that layer.",
    "Onset compartment refers to where the disease process begins, not to the compartment "
    "sampled: KPMP CKD is predominantly diabetic/hypertensive (glomerular onset) although "
    "TI and G compartments were both sampled.",
]
(RES / "table1_footnotes.txt").write_text("\n\n".join(FOOT), encoding="utf-8")

append_summary("21_table1 / Table 1", {
    "行数": len(t),
    "ファイル": "results/table1_datasets.csv, results/table1_footnotes.txt",
    "脚注": FOOT,
}, cfg)
