"""Table 2: 候補分子の4データセット横断まとめ（3.7 用）。

"not measured"（そのプラットフォーム/パネルに存在しない）と
"not significant"（測定されたが有意でない）を厳密に区別する。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary, alias_lookup  # noqa: E402

log = get_logger("28_table2")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
EXT = ROOT / "data" / "external"

GENES = ["FAM3D", "PTN", "ANXA3", "NPNT", "ANGPTL2", "THBS1", "SULF1", "ITGB6",
         "TXNIP", "SPP1", "MGP", "FBLN1", "NEDD4L", "HSD11B2"]
ERCB_CONTRAST = "GPL22945: 全CKD vs LD"
NM = "not measured"

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))
CAT2ALIAS = dict(zip(OMAP.cat_symbol, OMAP.mouse_aliases))


def fmt(lfc, p, sig=0.05):
    if lfc is None:
        return NM
    star = " *" if (p is not None and np.isfinite(p) and p < sig) else ""
    return f"{lfc:+.2f} ({p:.3g}){star}"


def kpmp_ti():
    df = pd.read_csv(EXT / "KPMP" / "DataLake_DEPs.txt", sep="\t")
    d = df[df.Comparison == "CKD.vs.HRT.in.TI"].copy()
    gn = d.Gene_name.astype(str).str.strip()
    keep = ~(d.Gene_name.isna() | (gn == "") | (gn.str.lower() == "nan"))
    d, gn = d[keep], gn[keep]
    d = d.assign(_s=gn.str.upper(), _a=d.LogFC.abs())
    d = d.sort_values(["_s", "Adj_pvalue", "_a"], ascending=[True, True, False])
    return d.drop_duplicates("_s").set_index("_s")[["LogFC", "Adj_pvalue"]]


def main():
    cat = DE.loc["cat_rna_ctx_late"]
    mp21 = DE.loc["m_prot_d21"]
    TI = kpmp_ti()
    ercb = pd.read_csv(RES / "ercb_all_contrasts.csv.gz", index_col=0)
    ercb = ercb[ercb.contrast == ERCB_CONTRAST]

    rows = []
    for g in GENES:
        rec = {"Gene": g}
        # --- cat cortex late ---
        if g in cat.index:
            rec["Cat cortex, late (log2FC, q)"] = fmt(float(cat.loc[g, "lfc"]),
                                                      float(cat.loc[g, "q"]))
            c_sign = np.sign(cat.loc[g, "lfc"])
        else:
            rec["Cat cortex, late (log2FC, q)"] = NM
            c_sign = None
        # --- mouse protein D21 (alias 総当たり) ---
        hit = alias_lookup(mp21.index, CAT2ALIAS.get(g), CAT2MOUSE.get(g)) \
            or alias_lookup(mp21.index, None, g)
        if hit is not None:
            rec["Mouse protein, Day 21 (log2FC, q)"] = fmt(float(mp21.loc[hit, "lfc"]),
                                                           float(mp21.loc[hit, "q"]))
            m_sign = np.sign(mp21.loc[hit, "lfc"])
        else:
            rec["Mouse protein, Day 21 (log2FC, q)"] = NM
            m_sign = None
        # --- human KPMP TI (protein) ---
        if g in TI.index:
            rec["Human KPMP TI, protein (log2FC, adj p)"] = fmt(float(TI.loc[g, "LogFC"]),
                                                                float(TI.loc[g, "Adj_pvalue"]))
            k_sign = np.sign(TI.loc[g, "LogFC"])
        else:
            rec["Human KPMP TI, protein (log2FC, adj p)"] = NM
            k_sign = None
        # --- human ERCB TI (transcript) ---
        if g in ercb.index:
            rec["Human ERCB TI, transcript (log2FC, q)"] = fmt(float(ercb.loc[g, "lfc"]),
                                                               float(ercb.loc[g, "q"]))
            e_sign = np.sign(ercb.loc[g, "lfc"])
        else:
            rec["Human ERCB TI, transcript (log2FC, q)"] = NM
            e_sign = None

        # --- 3種方向一致 ---
        # ヒトは KPMP TI を第一とし、未測定なら ERCB TI で代替する。
        h_sign, h_src = (k_sign, "KPMP") if k_sign is not None else (e_sign, "ERCB")
        if c_sign is None or m_sign is None or h_sign is None:
            miss = [n for n, v in [("cat", c_sign), ("mouse", m_sign), ("human", h_sign)]
                    if v is None]
            rec["Direction across three species"] = f"not evaluable ({'/'.join(miss)} {NM})"
        elif c_sign == m_sign == h_sign:
            rec["Direction across three species"] = \
                f"concordant ({'+' if c_sign > 0 else '-'}, human = {h_src})"
        else:
            s = lambda v: "+" if v > 0 else "-"  # noqa: E731
            rec["Direction across three species"] = \
                f"discordant (cat {s(c_sign)} / mouse {s(m_sign)} / human {s(h_sign)}, {h_src})"
        rows.append(rec)

    t = pd.DataFrame(rows)
    t.to_csv(RES / "table2_candidates.csv", index=False)
    log.info("Table 2: %d 行", len(t))

    FOOT = [
        f"Human ERCB values are the '{ERCB_CONTRAST}' contrast (all CKD diagnoses pooled versus "
        "living-donor controls) on platform GPL22945, chosen because that platform carries 18 of the "
        "21 control biopsies; contrasts were never built across platforms. Diagnosis-specific "
        "contrasts are in results/ercb_candidate_validation.csv.",
        "'not measured' means the gene is absent from that platform or protein panel; "
        "'not significant' is never abbreviated this way and can be read from the adjusted p value. "
        "An asterisk marks adjusted p < 0.05.",
        "Direction concordance uses cat cortical transcript, mouse Day-21 protein and human protein "
        "(KPMP TI); where the gene is absent from the KPMP panel the human ERCB transcript is "
        "substituted, and the source used is stated in the cell.",
        "Mouse symbols were resolved through the full alias list, so FAM3D is matched to Oit1 in the "
        "transcript annotation and to Fam3d in the proteome.",
    ]
    (RES / "table2_footnotes.txt").write_text("\n\n".join(FOOT), encoding="utf-8")
    append_summary("28_table2 / Table 2", {
        "行数": len(t), "ERCBコントラスト": ERCB_CONTRAST,
        "ファイル": "results/table2_candidates.csv, results/table2_footnotes.txt",
        "脚注": FOOT}, cfg)


if __name__ == "__main__":
    main()
