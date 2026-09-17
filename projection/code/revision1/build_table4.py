"""A-3. Table 4 を結果ファイルから作り直し、1 セル 1 値で書き出す。

これまで α / cos θ / ν を 1 セルに連結していたため、docx 側で新旧の値が混ざる事故が
起きた。行を分け、区間も推定値と別の行にして、どのセルにも数値が 2 つ以上並ばない形にする。
新しい数値は作らない。既存の結果ファイルの値をそのまま並べ替えるだけ。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
R = HERE / "results"
OUT_MD = HERE / "results" / "round8"
OUT_MD.mkdir(parents=True, exist_ok=True)

SPACES = [("Group A (2,016)", "", "Group A"),
          ("All 1:1 orthologues (7,897)", "_ortholog_all", "all 1:1 orthologues"),
          ("Matched Group B (1,632)", "_groupB_matched2", "matched Group B"),
          ("Matched Group A (1,632)", "_groupA_matched", "matched Group A")]
ALPHA_DIR = {"": "alpha_groupA", "_ortholog_all": "alpha_ortholog_all",
             "_groupB_matched2": "alpha_groupB_matched2",
             "_groupA_matched": "alpha_groupA_matched"}


def main() -> int:
    GS = pd.read_csv(R / "revision1" / "genespace_with_reactome.tsv", sep="\t").set_index("space")
    cols, rows = [], {}

    def put(label, values):
        rows[label] = values

    per = {}
    for col, tag, gsname in SPACES:
        cols.append(col)
        pr = pd.read_csv(R / ALPHA_DIR[tag] / "projection.tsv", sep="\t", index_col=0)
        ta = pd.read_csv(R / ALPHA_DIR[tag] / "time_association.tsv", sep="\t")
        t = ta[(ta.subset == "all_model_states") & (ta.term == "cos x log_time")].iloc[0]
        mn = pd.read_csv(R / f"mantel_s1{tag}" / "mantel_s1.tsv", sep="\t")
        d = mn[(mn.subset == "12 model states (time defined)")
               & (mn.predictor == "elapsed time |Δlog10 h|")
               & (mn.dissimilarity == "D_cos (1-cos)")].iloc[0]
        ce = pd.read_csv(R / f"reliability{tag}" / "pairs_ceilings.tsv", sep="\t")
        cs = ce[ce["class"] == "cross_species"]
        au = pd.read_csv(R / f"control_analysis{tag}" / "auc_permutation.tsv", sep="\t")
        c, m = pr.loc["cat_CKD12"], pr.loc["cat_med_CKD12"]
        thr = 1.0 / float(m["cos"])
        per[col] = {
            "Feline cortical CKD1/2: α": f"{float(c.alpha):.3f}",
            "Feline cortical CKD1/2: cos θ": f"{float(c['cos']):.3f}",
            "Feline cortical CKD1/2: ν": f"{float(c.nrm):.3f}",
            "Feline medullary CKD1/2: α": f"{float(m.alpha):.3f}",
            "Feline medullary CKD1/2: cos θ": f"{float(m['cos']):.3f}",
            "Feline medullary CKD1/2: ν": f"{float(m.nrm):.3f}",
            "α and cos θ ordered oppositely (the two states above)":
                "yes" if (m.alpha > c.alpha and m["cos"] < c["cos"]) else "no",
            "ν needed for α > 1, i.e. 1/cos θ (medullary)": f"{thr:.3f}",
            "margin, ν − 1/cos θ (medullary)": f"{float(m.nrm) - thr:+.3f}".replace("-", "−"),
            "α above the reference’s own 1.000": "yes" if float(m.alpha) > 1.0 else "no",
            "Mantel, elapsed time, direction (ρ)": f"{float(d.mantel_rho):+.3f}",
            "Mantel, elapsed time, direction (permutation p)": f"{float(d.p_perm):.4f}",
            "cos θ against log time, per state": f"{float(t.rho):+.3f}",
            "same, 95% interval lower": f"{float(t.boot_lo):+.2f}".replace("-", "−"),
            "same, 95% interval upper": f"{float(t.boot_hi):+.2f}".replace("-", "−"),
            "Cat–mouse cos θ, maximum": f"{float(cs.cos.max()):.3f}",
            "Lowest Spearman–Brown benchmark, cat–mouse": f"{float(cs.ceiling_SB.min()):.3f}",
            "Cat–mouse pairs above their benchmark, of 48":
                f"{int((cs.cos > cs.ceiling_SB).sum())}",
            "AUC over the 75 cross-dataset pairs (no permutation test, Section 2.2)":
                f"{float(au.auc_obs.iloc[0]):.4f}",
            "Pathway-level AUC, median": f"{float(GS.loc[gsname, 'pathway_auc_median']):.3f}",
            "AUC after aggregation, centred (Section 2.3)":
                f"{float(GS.loc[gsname, 'aggregated_auc_centred']):.3f}",
        }
    order = list(per[cols[0]].keys())
    lines = ["| Quantity | " + " | ".join(cols) + " |", "|" + "---|" * (len(cols) + 1)]
    tsv = ["\t".join(["quantity"] + cols)]
    bad = []
    for q in order:
        vals = [per[c][q] for c in cols]
        lines.append("| " + " | ".join([q] + vals) + " |")
        tsv.append("\t".join([q] + vals))
        for c, v in zip(cols, vals):
            n = len([tok for tok in v.replace("−", "-").split()
                     if tok.lstrip("+-").replace(".", "", 1).isdigit()])
            if n >= 2:
                bad.append((q, c, v))
    (OUT_MD / "table4.md").write_text("\n".join(lines) + "\n")
    (OUT_MD / "table4.tsv").write_text("\n".join(tsv) + "\n")
    print("\n".join(lines))
    print(f"\n行 {len(order)} x 列 {len(cols)}")
    print("連結セル（空白区切りの数値が 2 つ以上）:", bad if bad else "なし")
    assert not bad, bad
    return 0


if __name__ == "__main__":
    sys.exit(main())
