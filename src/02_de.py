"""全コントラストの差次的発現／存在量解析。

方針:
  - n=3 の素の t 検定は使わず、moderated t（経験ベイズ分散縮小）を用いる。
  - マウスRNAは群平均FPKM>=1 のフィルタ後に log2(FPKM+1)。
  - すべてのコントラストで leave-one-out の符号安定性を併記する。
  - ネコの「進行コントラスト」（CKD3/4 vs CKD1/2）を明示的に作る。ここが病期依存性の本体。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import (load_config, get_logger, moderated_ttest, hedges_g,  # noqa: E402
                       leave_one_out_contrast, append_summary)

log = get_logger("02_de")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
RES.mkdir(parents=True, exist_ok=True)
TH = cfg["thresholds"]


def read(name):
    df = pd.read_parquet(INT / f"{name}.parquet")
    return df.iloc[:, 0] if df.shape[1] == 1 and name.endswith("grp") else df


def contrast(mat: pd.DataFrame, groups: pd.Series, a: str, b: str, label: str) -> pd.DataFrame:
    A = mat[groups[groups == a].index]
    B = mat[groups[groups == b].index]
    res = moderated_ttest(A, B)
    res["hedges_g"] = hedges_g(A, B)
    loo = leave_one_out_contrast(A, B)
    res["loo_min"], res["loo_max"] = loo["loo_min"], loo["loo_max"]
    res["sign_stable"] = loo["sign_stable"]
    res["contrast"] = label
    log.info("%s: n=%d, prior_df=%.2f, q<%.2f: %d, |lfc|>%.2f&q<%.2f: %d",
             label, len(res), res.attrs.get("prior_df", np.nan), TH["q_sig"],
             int((res.q < TH["q_sig"]).sum()), TH["lfc_strong"],
             TH["q_sig"], int(((res.q < TH["q_sig"]) & (res.lfc.abs() > TH["lfc_strong"])).sum()))
    return res


def main():
    out: dict[str, pd.DataFrame] = {}

    # ---------- cat ----------
    for tis in ["ctx", "med"]:
        rna = read(f"cat_rna_{tis}")
        grp = read(f"cat_rna_{tis}_grp").squeeze()
        rna = rna[rna.std(axis=1) > 0]
        out[f"cat_rna_{tis}_early"] = contrast(rna, grp, "CKD1/2", "Control", f"cat_rna_{tis}_early")
        out[f"cat_rna_{tis}_late"] = contrast(rna, grp, "CKD3/4", "Control", f"cat_rna_{tis}_late")
        out[f"cat_rna_{tis}_prog"] = contrast(rna, grp, "CKD3/4", "CKD1/2", f"cat_rna_{tis}_prog")

        prot = read(f"cat_prot_{tis}")
        pgrp = read(f"cat_prot_{tis}_grp").squeeze()
        prot = prot.dropna(thresh=5)
        out[f"cat_prot_{tis}_early"] = contrast(prot, pgrp, "CKD1/2", "Control", f"cat_prot_{tis}_early")
        out[f"cat_prot_{tis}_late"] = contrast(prot, pgrp, "CKD3/4", "Control", f"cat_prot_{tis}_late")
        out[f"cat_prot_{tis}_prog"] = contrast(prot, pgrp, "CKD3/4", "CKD1/2", f"cat_prot_{tis}_prog")

    # ---------- mouse RNA ----------
    mr = read("mouse_rna")
    mgrp = read("mouse_rna_grp").squeeze()
    gm = mr.T.groupby(mgrp).mean().T
    keep = gm.max(axis=1) >= TH["fpkm_min"]
    log.info("マウスRNA発現フィルタ: %d / %d 遺伝子を保持", int(keep.sum()), len(mr))
    mrl = np.log2(mr[keep] + 1)
    for cond, label in [("DT100ng-5D", "m_rna_5d"), ("DT100ng-2W", "m_rna_2w"),
                        ("DT100ng-3W", "m_rna_3w"), ("DT250ng-2W", "m_rna_2w250")]:
        if (mgrp == cond).sum() == 0:
            log.warning("群 %s が見つかりません", cond)
            continue
        out[label] = contrast(mrl, mgrp, cond, "Ctrl", label)
    # マウスの進行コントラスト（2W→3W が本当に「進行」かを検証するため）
    out["m_rna_prog_5d2w"] = contrast(mrl, mgrp, "DT100ng-2W", "DT100ng-5D", "m_rna_prog_5d2w")
    out["m_rna_prog_2w3w"] = contrast(mrl, mgrp, "DT100ng-3W", "DT100ng-2W", "m_rna_prog_2w3w")
    out["m_rna_dose"] = contrast(mrl, mgrp, "DT250ng-2W", "DT100ng-2W", "m_rna_dose")

    # ---------- mouse protein ----------
    mp = read("mouse_prot").dropna(thresh=6)
    pgrp = read("mouse_prot_grp").squeeze()
    out["m_prot_d14"] = contrast(mp, pgrp, "Day14", "Ctrl", "m_prot_d14")
    out["m_prot_d21"] = contrast(mp, pgrp, "Day21", "Ctrl", "m_prot_d21")
    out["m_prot_prog"] = contrast(mp, pgrp, "Day21", "Day14", "m_prot_prog")

    # 存在量（順列帰無分布の十分位マッチ用）
    abundance = {
        "mouse_rna": np.log2(mr[keep] + 1).mean(axis=1),
        "mouse_prot": mp.mean(axis=1),
        "cat_rna_ctx": read("cat_rna_ctx").mean(axis=1),
        "cat_rna_med": read("cat_rna_med").mean(axis=1),
        "cat_prot_ctx": read("cat_prot_ctx").mean(axis=1),
        "cat_prot_med": read("cat_prot_med").mean(axis=1),
    }
    pd.DataFrame(abundance).to_parquet(INT / "abundance.parquet")

    allde = pd.concat({k: v for k, v in out.items()}, names=["dataset", "feature"])
    allde.to_parquet(INT / "de_all.parquet")
    allde.reset_index().to_csv(RES / "de_all.csv.gz", index=False, compression="gzip")

    summary = {}
    for k, v in out.items():
        sig = int(((v.q < TH["q_sig"]) & (v.lfc.abs() > TH["lfc_strong"])).sum())
        summary[k] = f"n={len(v)}, q<0.05&|FC|>1.5: {sig}, prior_df={v.attrs.get('prior_df', float('nan')):.1f}"
    append_summary("02_de / コントラスト別サマリ", summary, cfg)
    log.info("saved de_all")


if __name__ == "__main__":
    main()
