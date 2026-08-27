"""種間一致度の定量。

要点:
  - オルソログ表で cat→mouse を対応させる（大文字一致との差分も出す）。
  - 帰無分布は発現量十分位マッチの並べ替え。
  - 「種内ベンチマーク」を必ず併記する。種間rho単体では解釈できないため。
      同一個体 RNA-protein / 皮質-髄質 / 隣接時点 など。
  - 傷害軸遺伝子を除外しても rho が保たれるかを確認（少数の炎症遺伝子による見かけの一致でないこと）。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import (load_config, get_logger, decile_matched_permutation,  # noqa: E402
                       append_summary)

log = get_logger("03_crossspecies")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]

DE = pd.read_parquet(INT / "de_all.parquet")
ABU = pd.read_parquet(INT / "abundance.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))


def lfc(dataset: str) -> pd.Series:
    s = DE.loc[dataset, "lfc"]
    return s[np.isfinite(s)]


def to_mouse_space(s: pd.Series, how: str = "ortholog") -> pd.Series:
    """ネコの系列をマウス遺伝子シンボル空間へ写像する。"""
    if how == "ortholog":
        idx = pd.Series(s.index).map(CAT2MOUSE)
        out = s.copy()
        out.index = idx
        out = out[out.index.notna()]
    else:  # naive: 大文字一致（比較用）
        out = s.copy()
        out.index = [str(i).upper() for i in out.index]
    return out[~out.index.duplicated()]


def norm_index(s: pd.Series) -> pd.Series:
    out = s.copy()
    out.index = [str(i).upper() for i in out.index]
    return out[~out.index.duplicated()]


def concord(cat_ds: str, mouse_ds: str, how: str = "ortholog", label: str = None):
    x = norm_index(to_mouse_space(lfc(cat_ds), how))
    y = norm_index(lfc(mouse_ds))
    idx = x.index.intersection(y.index)
    x, y = x.loc[idx], y.loc[idx]
    rho, p = stats.spearmanr(x, y)
    m = (x.abs() > np.log2(1.2)) & (y.abs() > np.log2(1.2))
    agree = float((np.sign(x[m]) == np.sign(y[m])).mean()) if m.sum() > 20 else np.nan
    return {"label": label or f"{cat_ds} x {mouse_ds}", "mapping": how,
            "n": len(idx), "rho": float(rho), "p": float(p),
            "dir_agree": agree, "n_agree_set": int(m.sum())}


def within_species(a: str, b: str, label: str):
    """同一種内の比較（物差し用）。両方ともそのままの空間で比較する。"""
    x, y = norm_index(lfc(a)), norm_index(lfc(b))
    idx = x.index.intersection(y.index)
    rho, p = stats.spearmanr(x.loc[idx], y.loc[idx])
    return {"label": label, "n": len(idx), "rho": float(rho), "p": float(p)}


INJURY_AXIS = set("""COL1A1 COL1A2 COL3A1 COL4A1 FN1 POSTN SPARC LUM DCN THBS2 TNC LOX TIMP1 MMP2
SERPINE1 ACTA2 TAGLN VIM CCN2 SPP1 BGN C1QA C1QB C1QC C1R C1S C3 CFB CFH CFD C3AR1 C5AR1 SERPING1
CD68 AIF1 ITGAM ITGAX CSF1R LYZ TYROBP FCER1G MRC1 MSR1 CTSS CTSK LAPTM5 PTPRC IL1B CCL2 CXCL10
HAVCR1 LCN2 CLU KRT8 KRT18 KRT19 VCAM1 CDKN1A ANXA2 S100A4 S100A6
SLC34A1 SLC5A2 SLC22A6 SLC22A8 LRP2 CUBN AQP1 SLC4A4 SLC9A3
NDUFA1 NDUFA4 NDUFS1 SDHA SDHB UQCRC1 COX4I1 COX5A ATP5F1A ATP5F1B
CPT1A CPT2 ACADM ACADL ACADVL HADHA HADHB ACOX1 PPARA PPARGC1A""".split())


def main():
    rows = []
    for mouse_ds in ["m_rna_5d", "m_rna_2w", "m_rna_3w", "m_rna_2w250"]:
        for cat_ds in ["cat_rna_ctx_late", "cat_rna_ctx_early", "cat_rna_med_late", "cat_rna_med_early"]:
            for how in ["ortholog", "naive"]:
                rows.append(concord(cat_ds, mouse_ds, how))
    for mouse_ds in ["m_prot_d14", "m_prot_d21"]:
        for cat_ds in ["cat_prot_ctx_late", "cat_prot_ctx_early", "cat_prot_med_late"]:
            for how in ["ortholog", "naive"]:
                rows.append(concord(cat_ds, mouse_ds, how))
    # 層をまたぐ比較
    rows.append(concord("cat_rna_ctx_late", "m_prot_d21", "ortholog", "catRNA x mouseProt"))
    rows.append(concord("cat_prot_ctx_late", "m_rna_2w", "ortholog", "catProt x mouseRNA"))
    xs = pd.DataFrame(rows)
    xs.to_csv(RES / "crossspecies_concordance.csv", index=False)

    # マッピング法の差分
    piv = xs[xs.mapping.notna()].pivot_table(index="label", columns="mapping", values=["rho", "n"])
    piv.to_csv(RES / "crossspecies_mapping_comparison.csv")
    log.info("オルソログ vs 大文字一致の rho 差分:\n%s",
             (piv[("rho", "ortholog")] - piv[("rho", "naive")]).describe())

    # ---------- 種内ベンチマーク ----------
    bench = [
        within_species("m_prot_d14", "m_prot_d21", "マウス protein D14 vs D21"),
        within_species("m_rna_2w", "m_rna_2w250", "マウス RNA 2W DT100 vs DT250"),
        within_species("m_rna_2w", "m_rna_3w", "マウス RNA 2W vs 3W"),
        within_species("m_rna_5d", "m_rna_2w", "マウス RNA 5D vs 2W"),
        within_species("cat_rna_ctx_late", "cat_rna_ctx_early", "ネコ皮質 RNA 晩期 vs 早期"),
        within_species("cat_rna_ctx_late", "cat_rna_med_late", "ネコ 皮質 vs 髄質 RNA"),
        within_species("cat_prot_ctx_late", "cat_prot_med_late", "ネコ 皮質 vs 髄質 protein"),
        within_species("cat_rna_ctx_late", "cat_prot_ctx_late", "ネコ皮質 RNA vs protein"),
        within_species("m_rna_2w", "m_prot_d14", "マウス RNA2W vs protein D14"),
        within_species("m_rna_3w", "m_prot_d21", "マウス RNA3W vs protein D21"),
    ]
    bdf = pd.DataFrame(bench)
    bdf.to_csv(RES / "within_species_benchmarks.csv", index=False)

    # ---------- 十分位マッチ順列検定 ----------
    perm = []
    for cat_ds, mouse_ds, abu in [
        ("cat_rna_ctx_late", "m_rna_2w", "mouse_rna"),
        ("cat_prot_ctx_late", "m_prot_d21", "mouse_prot"),
        ("cat_rna_med_late", "m_rna_2w", "mouse_rna"),
    ]:
        x = norm_index(to_mouse_space(lfc(cat_ds)))
        y = norm_index(lfc(mouse_ds))
        a = norm_index(ABU[abu].dropna())
        r = decile_matched_permutation(x, y, a, n_perm=cfg["permutation"]["n_perm_rho"],
                                       n_decile=cfg["permutation"]["n_decile"], seed=SEED)
        r["label"] = f"{cat_ds} x {mouse_ds}"
        perm.append(r)
        log.info("%s: rho=%.3f null=%.4f±%.4f z=%.1f p=%.4g",
                 r["label"], r["rho"], r["null_mean"], r["null_sd"], r["z"], r["p_perm"])
    pd.DataFrame(perm).to_csv(RES / "crossspecies_permutation.csv", index=False)

    # ---------- 傷害軸を除いても一致は保たれるか ----------
    drop_rows = []
    for cat_ds, mouse_ds in [("cat_rna_ctx_late", "m_rna_2w"), ("cat_prot_ctx_late", "m_prot_d21")]:
        x = norm_index(to_mouse_space(lfc(cat_ds)))
        y = norm_index(lfc(mouse_ds))
        idx = x.index.intersection(y.index)
        keep = [g for g in idx if g not in INJURY_AXIS]
        drop_rows.append({
            "label": f"{cat_ds} x {mouse_ds}",
            "rho_all": float(stats.spearmanr(x.loc[idx], y.loc[idx])[0]),
            "rho_wo_injury": float(stats.spearmanr(x.loc[keep], y.loc[keep])[0]),
            "n_all": len(idx), "n_kept": len(keep),
        })
    pd.DataFrame(drop_rows).to_csv(RES / "crossspecies_injury_dropout.csv", index=False)

    append_summary("03_crossspecies", {
        "種間 RNA（皮質晩期 x マウス2W, ortholog）":
            xs.query("label=='cat_rna_ctx_late x m_rna_2w' and mapping=='ortholog'")["rho"].round(3).tolist(),
        "種間 protein（皮質晩期 x D21, ortholog）":
            xs.query("label=='cat_prot_ctx_late x m_prot_d21' and mapping=='ortholog'")["rho"].round(3).tolist(),
        "種内ベンチマーク": {r["label"]: round(r["rho"], 3) for r in bench},
        "順列検定": {r["label"]: f"rho={r['rho']:.3f}, z={r['z']:.1f}, p={r['p_perm']:.4g}" for r in perm},
        "傷害軸除外": drop_rows,
    }, cfg)


if __name__ == "__main__":
    main()
