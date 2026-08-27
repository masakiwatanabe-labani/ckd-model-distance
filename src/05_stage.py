"""病期トラジェクトリ解析。

検証対象:
  A. ネコ皮質は IRIS1/2 で分子病態がほぼ完成しているか（早期起動 vs 後期起動）
  B/C. 髄質は非単調か（早期ピーク→後退、符号反転）
  D. マウスの時間軸はネコの進行軸と一致するか（5D→2W、2W→3W、用量）
  E. マウスはネコのどちらの病期に対応するか（偏相関）

注意:
  - 早期起動／後期起動は「晩期に有意な遺伝子の、早期での達成率」で定義する。
    「早期に個々の遺伝子が有意」という意味ではない。原著は皮質早期DEGを6個としており、
    両者は矛盾しない（検出力の差）。この点は必ず結果に明記すること。
  - 低発現バイアスの確認（後期起動集合が発現下限に偏っていないか）を必ず出す。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("05_stage")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
TH = cfg["thresholds"]

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))


def tab(ds: str) -> pd.DataFrame:
    return DE.loc[ds].copy()


def classify(tissue: str) -> pd.DataFrame:
    E, L = tab(f"cat_rna_{tissue}_early"), tab(f"cat_rna_{tissue}_late")
    P = tab(f"cat_rna_{tissue}_prog")
    idx = E.index.intersection(L.index)
    t = pd.DataFrame({
        "lfc_early": E.loc[idx, "lfc"], "p_early": E.loc[idx, "p"], "q_early": E.loc[idx, "q"],
        "lfc_late": L.loc[idx, "lfc"], "p_late": L.loc[idx, "p"], "q_late": L.loc[idx, "q"],
        "lfc_prog": P.loc[idx, "lfc"], "p_prog": P.loc[idx, "p"],
        "prog_sign_stable": P.loc[idx, "sign_stable"],
    }).dropna(subset=["lfc_early", "lfc_late"])

    sig_late = (t.lfc_late.abs() > TH["lfc_strong"]) & (t.q_late < TH["q_sig"])
    ratio = t.lfc_early / t.lfc_late
    same = np.sign(t.lfc_early) == np.sign(t.lfc_late)
    t["class"] = "other"
    t.loc[sig_late & same & (ratio >= TH["early_onset_ratio"]), "class"] = "early_onset"
    t.loc[sig_late & (ratio < TH["late_onset_ratio"]), "class"] = "late_onset"
    t.loc[sig_late & same & (t.lfc_early.abs() > 1.5 * t.lfc_late.abs()), "class"] = "transient_peak"
    both_big = (t.lfc_early.abs() > TH["lfc_strong"]) & (t.lfc_late.abs() > TH["lfc_strong"])
    t.loc[both_big & ~same, "class"] = "sign_reversal"
    t["ratio_early_late"] = ratio
    return t


def low_expression_check(t: pd.DataFrame, tissue: str) -> dict:
    mat = pd.read_parquet(INT / f"cat_rna_{tissue}.parquet")
    pct = mat.mean(axis=1).rank(pct=True) * 100
    out = {}
    for cls in ["early_onset", "late_onset", "transient_peak", "sign_reversal"]:
        sub = t.index[t["class"] == cls]
        v = pct.reindex(sub).dropna()
        out[cls] = {"n": len(sub), "median_pct": round(float(v.median()), 1) if len(v) else None,
                    "frac_below_25pct": round(float((v < 25).mean()), 3) if len(v) else None}
    out["all_genes_median_pct"] = 50.0
    return out


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series):
    idx = x.index.intersection(y.index).intersection(z.index)
    a, b, c = x.loc[idx], y.loc[idx], z.loc[idx]
    ok = a.notna() & b.notna() & c.notna() & np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    a, b, c = (stats.rankdata(v[ok]) for v in (a, b, c))
    rab = np.corrcoef(a, b)[0, 1]
    rac = np.corrcoef(a, c)[0, 1]
    rbc = np.corrcoef(b, c)[0, 1]
    return float((rab - rac * rbc) / np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))), int(ok.sum())


def to_mouse(s: pd.Series) -> pd.Series:
    idx = pd.Series(s.index).map(CAT2MOUSE)
    out = s.set_axis(idx)
    out = out[out.index.notna()]
    out.index = [str(i).upper() for i in out.index]
    return out[~out.index.duplicated()]


def norm(s: pd.Series) -> pd.Series:
    out = s.copy()
    out.index = [str(i).upper() for i in out.index]
    return out[~out.index.duplicated()]


def main():
    payload = {}
    for tissue in ["ctx", "med"]:
        t = classify(tissue)
        t.to_csv(RES / f"stage_classes_{tissue}.csv")
        counts = t["class"].value_counts().to_dict()
        log.info("%s: %s", tissue, counts)
        payload[f"{tissue} クラス内訳"] = counts
        payload[f"{tissue} 発現量バイアス確認"] = low_expression_check(t, tissue)
        # 符号反転・早期ピークは進行コントラストで支持されるかを確認
        for cls in ["sign_reversal", "transient_peak"]:
            sub = t[(t["class"] == cls)]
            supported = sub[(sub.p_prog < 0.05) & sub.prog_sign_stable.fillna(False)]
            payload[f"{tissue} {cls} のうち進行コントラストで支持"] = f"{len(supported)}/{len(sub)}"
            supported.to_csv(RES / f"stage_{cls}_{tissue}_supported.csv")

    # 早期/晩期の全体相関（皮質は高く、髄質は低いはず）
    for tissue in ["ctx", "med"]:
        E, L = tab(f"cat_rna_{tissue}_early")["lfc"], tab(f"cat_rna_{tissue}_late")["lfc"]
        idx = E.index.intersection(L.index)
        payload[f"{tissue} rho(early,late)"] = round(float(stats.spearmanr(E[idx], L[idx])[0]), 3)

    # ---------- 進行軸の種間比較 ----------
    prog_rows = []
    cat_prog = to_mouse(tab("cat_rna_ctx_prog")["lfc"])
    cat_prog_med = to_mouse(tab("cat_rna_med_prog")["lfc"])
    for mouse_ds, label in [("m_rna_prog_5d2w", "マウス 5D→2W"),
                            ("m_rna_prog_2w3w", "マウス 2W→3W"),
                            ("m_rna_dose", "マウス DT100→DT250"),
                            ("m_rna_2w", "マウス 2W vs Ctrl（参考）")]:
        y = norm(tab(mouse_ds)["lfc"])
        for cat_axis, cat_label in [(cat_prog, "ネコ皮質進行軸"), (cat_prog_med, "ネコ髄質進行軸")]:
            idx = cat_axis.index.intersection(y.index)
            rho, p = stats.spearmanr(cat_axis[idx], y[idx])
            prog_rows.append({"cat_axis": cat_label, "mouse_axis": label,
                              "rho": float(rho), "p": float(p), "n": len(idx)})
    pd.DataFrame(prog_rows).to_csv(RES / "progression_axis_crossspecies.csv", index=False)
    payload["進行軸の種間一致"] = prog_rows

    # ---------- 偏相関：マウスはどちらの病期に対応するか ----------
    part = []
    e = to_mouse(tab("cat_rna_ctx_early")["lfc"])
    l = to_mouse(tab("cat_rna_ctx_late")["lfc"])
    for ds in ["m_rna_5d", "m_rna_2w", "m_rna_3w", "m_rna_2w250"]:
        y = norm(tab(ds)["lfc"])
        r_late, n = partial_spearman(y, l, e)
        r_early, _ = partial_spearman(y, e, l)
        part.append({"mouse": ds, "partial_rho_late_given_early": round(r_late, 3),
                     "partial_rho_early_given_late": round(r_early, 3), "n": n})
    pd.DataFrame(part).to_csv(RES / "stage_partial_correlation.csv", index=False)
    payload["偏相関"] = part
    payload["注意"] = ("早期側の推定は分散が大きく減衰するため、偏相関の非対称性は"
                       "そのまま「晩期に対応」と読んではいけない")

    append_summary("05_stage", payload, cfg)


if __name__ == "__main__":
    main()
