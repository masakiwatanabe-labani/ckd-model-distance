"""CHECK 1b: rho 不一致の原因を「遺伝子空間」以外の要因に分解する。

CHECK 1 で、大文字一致で追加された遺伝子（added）は発現量をマッチさせると
core と区別がつかないことが分かった。つまり遺伝子空間の違いでは
旧稿 0.490 と新稿の値の差は説明できない。

残る候補は Δ の作り方そのもの。新旧で3点が違う:
    (1) 写像       : Ensembl ortholog（旧） vs 大文字一致（新）
    (2) 発現フィルタ: 群平均 FPKM>=1（旧） vs フィルタなし（新）
    (3) Δ の定義   : 群平均の log2 差（旧） vs 遺伝子ごとに z 化してから差（新）

(3) の z 化は new pipeline の compute_delta_vectors(normalize=True) にある
    delta_base = x - mean(ctrl);  z = (delta_base - mean) / sd;  Δ = mean_grp(z) - mean_ctrl(z)
で、これは log2FC ではなく標準化効果量。遺伝子間の順位が変わるので Spearman が動く。

2x2x2 の全セルを回して、どの要因がどれだけ rho を動かすかを出す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
from lib_stats import load_config  # noqa: E402

cfg = load_config(ROOT / "config" / "config.yaml")
INT = ROOT / cfg["paths"]["interim"]
REF = ROOT / cfg["paths"]["ref"]
FPKM_MIN = cfg["thresholds"]["fpkm_min"]


def o2o(fname):
    o = pd.read_csv(REF / fname, sep="\t")
    o = o[o.orthology_type.isin(cfg["orthology"]["keep_types"])]
    o = o.dropna(subset=["src_symbol", "tgt_symbol"]).copy()
    o["rank"] = (-o.confidence.fillna(0)) * 1000 - o.perc_id.fillna(0)
    return dict(zip(*zip(*o.sort_values("rank").drop_duplicates("src_symbol")
                         [["src_symbol", "tgt_symbol"]].to_numpy())))


C2H, M2H = o2o("orthologs_cat2human.tsv"), o2o("orthologs_mouse2human.tsv")


def read(n):
    d = pd.read_parquet(INT / f"{n}.parquet")
    return d.iloc[:, 0] if n.endswith("grp") else d


def delta(mat, grp, case_label, ctrl_label, zscore):
    """Δ ベクトル。zscore=False なら群平均の差（log2FC）、True なら遺伝子ごとz化後の差。"""
    cm = grp == ctrl_label
    gm = grp == case_label
    X = mat.loc[:, mat.columns.isin(grp.index)]
    if not zscore:
        return X.loc[:, gm[X.columns].to_numpy()].mean(axis=1) - \
               X.loc[:, cm[X.columns].to_numpy()].mean(axis=1)
    base = X.sub(X.loc[:, cm[X.columns].to_numpy()].mean(axis=1), axis=0).fillna(0)
    z = base.sub(base.mean(axis=1), axis=0).div(base.std(axis=1, ddof=1) + 1e-8, axis=0)
    return z.loc[:, gm[X.columns].to_numpy()].mean(axis=1) - \
           z.loc[:, cm[X.columns].to_numpy()].mean(axis=1)


def to_space(s, species, how):
    m = {"cat": C2H, "mouse": M2H}[species]
    idx = [str(m[g]).upper() if (how == "ortholog" and g in m) else str(g).upper()
           for g in s.index]
    if how == "ortholog":
        keep = [g in m for g in s.index]
        s, idx = s[keep], [i for i, k in zip(idx, keep) if k]
    out = pd.Series(s.to_numpy(), index=idx)
    return out[~out.index.duplicated(keep="first")]


def main():
    cat = read("cat_rna_ctx"); cgrp = read("cat_rna_ctx_grp").squeeze()
    cat = cat[cat.std(axis=1) > 0]
    mr = read("mouse_rna"); mgrp = read("mouse_rna_grp").squeeze()

    gm = mr.T.groupby(mgrp).mean().T
    filt = gm.max(axis=1) >= FPKM_MIN

    rows = []
    for how in ["ortholog", "naive"]:
        for use_filter in [True, False]:
            m = np.log2(mr[filt] + 1) if use_filter else np.log2(mr + 1)
            for z in [False, True]:
                x = to_space(delta(cat, cgrp, "CKD3/4", "Control", z), "cat", how)
                y = to_space(delta(m, mgrp, "DT100ng-2W", "Ctrl", z), "mouse", how)
                idx = x.index.intersection(y.index)
                xv, yv = x.loc[idx], y.loc[idx]
                ok = np.isfinite(xv) & np.isfinite(yv)
                rows.append({"mapping": how,
                             "fpkm_filter": "on" if use_filter else "off",
                             "delta": "z-scored" if z else "log2FC",
                             "n": int(ok.sum()),
                             "rho": float(stats.spearmanr(xv[ok], yv[ok])[0])})
                print(f"{how:9s} filter={'on ' if use_filter else 'off'} "
                      f"{'z-scored' if z else 'log2FC  '}  n={int(ok.sum()):6d}  "
                      f"rho={rows[-1]['rho']:+.4f}")

    df = pd.DataFrame(rows)
    df.insert(0, "variant", [f"S{i+1}" for i in range(len(df))])
    df["primary"] = ((df.mapping == "ortholog") & (df.fpkm_filter == "on")
                     & (df["delta"] == "log2FC"))
    out = HERE / "results" / "check1"
    out.mkdir(parents=True, exist_ok=True)
    pass

    base_rho = df[df.primary].rho.iloc[0]
    base = base_rho
    new = df[(df.mapping == "naive") & (df.fpkm_filter == "off") & (df["delta"] == "z-scored")].rho.iloc[0]
    newms = df[(df.mapping == "naive") & (df.fpkm_filter == "off") & (df["delta"] == "log2FC")].rho.iloc[0]

    # 各要因の主効果（他2要因を平均した上での差）
    eff = {}
    for col, a, b in [("mapping", "naive", "ortholog"),
                      ("fpkm_filter", "off", "on"),
                      ("delta", "z-scored", "log2FC")]:
        eff[col] = float(df[df[col] == a].rho.mean() - df[df[col] == b].rho.mean())

    df["delta_vs_primary"] = df.rho - base_rho
    lines = ["# Supplementary: Δ の作り方に対する感度分析（cat_CKD34 x mouse_2W）", "",
             "主解析は S1（ortholog / FPKM>=1 / log2FC）に固定する。", "",
             "| 変種 | 写像 | 発現フィルタ | Δ の定義 | n | Spearman rho | 主解析との差 | 主解析 |",
             "|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(f"| {r['variant']} | {r['mapping']} | {r['fpkm_filter']} | {r['delta']} "
                     f"| {r['n']} | {r['rho']:+.4f} | {r['delta_vs_primary']:+.4f} "
                     f"| {'**はい**' if r['primary'] else ''} |")
    lines += ["", "## 参照点", "",
              f"- 主解析 S1（ortholog / filter on / log2FC）: rho = {base:+.4f}（旧稿 0.490）",
              f"- naive / filter off / log2FC:  rho = {newms:+.4f}（新稿 0.398）",
              f"- naive / filter off / z-scored: rho = {new:+.4f}",
              f"- S1 との差: {newms - base:+.4f} / {new - base:+.4f}", "",
              "## 要因ごとの主効果（他2要因で平均した差）", ""]
    for k, v in sorted(eff.items(), key=lambda kv: -abs(kv[1])):
        lines.append(f"- {k}: {v:+.4f}")
    df.round(6).to_csv(out / "supplementary_sensitivity.tsv", sep="\t", index=False)
    (out / "supplementary_sensitivity.md").write_text("\n".join(lines) + "\n")
    print()
    print("\n".join(lines[-8:]))


if __name__ == "__main__":
    main()
