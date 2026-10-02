# -*- coding: utf-8 -*-
"""BE1. ベンチマークのシミュレーションの超過率を、設計別・条件別に出す。

§4.10 の段落が述べる数値はすべてこの表の行に対応させる。手で集計しない。

  panel        図 S4 のどのパネルの量か（i / ii / iii / iv）
  metric       結果ファイルの列名
  design       pooled / feline-like (7 vs 6) / mouse-like (3 vs 6)
  true_cos     応答どうしの真のコサイン
  shift        共通平均のずれ（pooled なら all）
  reliability  集計した信頼性の範囲
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
SRC = HERE / "results" / "round5" / "ceiling_decomposition_sim.tsv"
OUT = HERE / "results" / "roundR" / "simulation_exceedance.tsv"

# 帯は観測された split-half 信頼性の範囲（Methods 4.10）
LO, HI = 0.568, 0.986

PANELS = [("i", "frac_i_pearson_matched_n", "Centered, sizes matched"),
          ("ii", "frac_ii_pearson_half_vs_full", "Centered, half-sample benchmark"),
          ("iii", "frac_iii_cosine_half_vs_full", "Uncentered, half-sample benchmark"),
          ("iv", "frac_v_cosine_SB_vs_full", "Uncentered, corrected benchmark")]


def main() -> int:
    T = pd.read_csv(SRC, sep="\t")
    band = (T.mean_r_half_cosine >= LO) & (T.mean_r_half_cosine <= HI)
    rows = []
    for tc in (0.9, 1.0):
        for tag, col, what in PANELS:
            d = T[band & (T.true_cos == tc)]
            groups = [("pooled", "all", d)]
            for g, e in d.groupby("design"):
                groups.append((g, "all", e))
            # 図が描く shift ごとにも出す（図と本文の対応を取るため）
            for g, e in d.groupby(["design", "shift"]):
                groups.append((g[0], f"{g[1]:.2f}", e))
            for design, shift, e in groups:
                rows.append({"panel": tag, "metric": col, "what": what,
                             "design": design, "true_cos": tc, "shift": shift,
                             "reliability": f"{LO}-{HI}", "n": len(e),
                             "mean": round(float(e[col].mean()), 4),
                             "min": round(float(e[col].min()), 4),
                             "max": round(float(e[col].max()), 4),
                             "spread_pp": round(float(e[col].max() - e[col].min()) * 100, 1)})
    # 非同一応答（true_cos < 1）が補正後ベンチマークを超えないこと
    h = T[band & (T.true_cos < 1.0)]
    rows.append({"panel": "iv", "metric": "frac_v_cosine_SB_vs_full",
                 "what": "Uncentered, corrected benchmark", "design": "pooled",
                 "true_cos": "<1.0", "shift": "all", "reliability": f"{LO}-{HI}",
                 "n": len(h), "mean": round(float(h.frac_v_cosine_SB_vs_full.mean()), 4),
                 "min": round(float(h.frac_v_cosine_SB_vs_full.min()), 4),
                 "max": round(float(h.frac_v_cosine_SB_vs_full.max()), 4),
                 "spread_pp": 0.0})
    # 真のコサイン 0.9 では半標本基準の超過が信頼性の階段関数になる。
    # 帯内平均（39.7 / 39.9%）は格子点の数え上げなので、階段の位置も残す。
    for tag, col in (("ii", "frac_ii_pearson_half_vs_full"),
                     ("iii", "frac_iii_cosine_half_vs_full")):
        d = T[band & (T.true_cos == 0.9)]
        hi, lo = d[d[col] > 0.5], d[d[col] == 0.0]
        assert len(hi) + len(lo) == len(d), "中間の点がある。階段としては書けない"
        rows.append({"panel": tag, "metric": col + " (exceeding grid points)",
                     "what": "reliability where the benchmark is exceeded",
                     "design": "pooled", "true_cos": 0.9, "shift": "all",
                     "reliability": f"{hi.mean_r_half_cosine.min():.3f}-"
                                    f"{hi.mean_r_half_cosine.max():.3f}",
                     "n": len(hi), "mean": round(float(hi[col].mean()), 4),
                     "min": round(float(hi[col].min()), 4),
                     "max": round(float(hi[col].max()), 4), "spread_pp": 0.0})
        rows.append({"panel": tag, "metric": col + " (non-exceeding grid points)",
                     "what": "reliability where the benchmark is never exceeded",
                     "design": "pooled", "true_cos": 0.9, "shift": "all",
                     "reliability": f"{lo.mean_r_half_cosine.min():.3f}-"
                                    f"{lo.mean_r_half_cosine.max():.3f}",
                     "n": len(lo), "mean": 0.0, "min": 0.0, "max": 0.0, "spread_pp": 0.0})

    # 経験的な shift の中央値に最も近い格子点（図が描く shift）
    shifts = sorted(T["shift"].unique())
    drawn = min(shifts, key=lambda s: abs(s - 0.257))

    S = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    S.to_csv(OUT, sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("### §4.10 の段落が引く行（信頼性 0.568–0.986 の帯の中）\n")
    key = S[(S["shift"] == "all") & (S.true_cos.isin([0.9, 1.0]))]
    print(key[["panel", "what", "design", "true_cos", "reliability", "n",
               "mean", "min", "max", "spread_pp"]].to_string(index=False))
    print(f"\n図が描く shift: {drawn}（経験的中央値 0.257 に最も近い格子点。"
          f"格子は {shifts}）")
    print(f"書き出し: {OUT}（{len(S)} 行）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
