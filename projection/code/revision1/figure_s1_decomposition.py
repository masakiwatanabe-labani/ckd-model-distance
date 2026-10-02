"""ベンチマークのシミュレーションの図（C案では Figure S4）。

4 条件を並べる。(i) 標本数をそろえた中心化ピアソン、(ii) 半標本基準 vs 全標本観測の
中心化ピアソン、(iii) 同じ比較を非中心化コサインで、(iv) 同じ比較を Spearman-Brown
補正済み基準で。SB 補正は半標本の信頼性を全標本に外挿する式なので、(iv) が (i) に
戻るなら、標本数の不一致は補正で解消できることになる。

BE2: 本文が設計別の超過率（ネコ型 59.4% / マウス型 98.0%）を報告するので、両設計を
描く。実線がマウス型（3 対 6）、破線がネコ型（7 対 6）。shift は経験的中央値 0.257 に
最も近い格子点 0.26 に固定する。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from make_figures import S1, S2, S3, INK, INK2, FIGW, clean, savefig  # noqa: E402

OUT = HERE / "manuscript" / "figures"
T = pd.read_csv(HERE / "results" / "round5" / "ceiling_decomposition_sim.tsv", sep="\t")

PANELS = [("frac_i_pearson_matched_n", "(i)  Centered, sizes matched"),
          ("frac_ii_pearson_half_vs_full", "(ii)  Centered, half-sample benchmark"),
          ("frac_iii_cosine_half_vs_full", "(iii)  Uncentered, half-sample benchmark"),
          ("frac_v_cosine_SB_vs_full", "(iv)  Uncentered, corrected benchmark")]
PALETTE = [S1, S2, S3, INK2, "#7a5cc4"]

SHIFT = 0.26          # 経験的中央値 0.257 に最も近い格子点
BAND = (0.568, 0.986)  # 観測された split-half 信頼性の範囲
# 実線 = マウス型、破線 = ネコ型。順序は描き順（マウス型を上に出す）
DESIGNS = [("mouse-like (3 vs 6)", "-", "o", 4.0, "full"),
           ("feline-like (7 vs 6)", (0, (3, 2)), "o", 4.2, "none")]

PRINTED: list[dict] = []


def main() -> int:
    d = T[T["shift"] == SHIFT]
    assert set(d.design) == {g for g, *_ in DESIGNS}, sorted(set(d.design))
    true_cos = sorted(d.true_cos.unique(), reverse=True)

    fig, axes = plt.subplots(2, 2, figsize=(FIGW, 5.0), sharey=True, sharex=True)
    for ax, (col, title) in zip(axes.ravel(), PANELS):
        for design, ls, mk, ms, fill in DESIGNS:
            e = d[d.design == design]
            for c, tc in zip(PALETTE, true_cos):
                g = e[e.true_cos == tc].sort_values("mean_r_half_cosine")
                ax.plot(g.mean_r_half_cosine, g[col], marker=mk, markersize=ms,
                        markerfacecolor=c if fill == "full" else "white",
                        markeredgecolor=c, markeredgewidth=1.0,
                        linewidth=1.5, linestyle=ls, color=c,
                        label=f"{design} | true cos {tc:.1f}")
                inside = g[(g.mean_r_half_cosine >= BAND[0]) & (g.mean_r_half_cosine <= BAND[1])]
                if len(inside):
                    PRINTED.append({"panel": title.split()[0].strip("()"), "metric": col,
                                    "design": design, "true_cos": tc,
                                    "n_in_band": len(inside),
                                    "mean_in_band": round(float(inside[col].mean()), 4),
                                    "min_in_band": round(float(inside[col].min()), 4),
                                    "max_in_band": round(float(inside[col].max()), 4)})
        ax.axhline(0.01, color=INK, linestyle=(0, (4, 3)), linewidth=1.0)
        ax.axvspan(*BAND, color="#eceae4", zorder=0)
        ax.set_title(title, loc="left", pad=4)
        ax.set_xlabel("reliability")
        ax.set_xlim(0.1, 1.0)
        clean(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of replicates with\nobserved above benchmark")
    for ax in axes[0, :]:
        ax.set_xlabel("")

    # 凡例は 2 本立て。色が真のコサイン、線種が設計。
    h_col = [Line2D([], [], color=c, linewidth=1.5, marker="o", markersize=4.0,
                    label=f"{tc:.1f}") for c, tc in zip(PALETTE, true_cos)]
    h_des = [Line2D([], [], color=INK2, linewidth=1.5, linestyle=ls, marker=mk,
                    markersize=ms, markerfacecolor=INK2 if fill == "full" else "white",
                    markeredgecolor=INK2, label=design)
             for design, ls, mk, ms, fill in DESIGNS]
    l1 = fig.legend(handles=h_col, title="true cos θ between the responses",
                    loc="lower center", ncol=5, frameon=False, handlelength=1.6,
                    columnspacing=1.4, bbox_to_anchor=(0.5, 0.055))
    fig.add_artist(l1)
    fig.legend(handles=h_des, title="sample-size design", loc="lower center", ncol=2,
               frameon=False, handlelength=2.4, columnspacing=2.0,
               bbox_to_anchor=(0.5, -0.025))

    fig.tight_layout(rect=(0, 0.16, 1, 1), h_pad=1.0, w_pad=0.8)
    print("FigS1:", savefig(fig, OUT / "FigS1_cosine_attenuation")[0].name)
    P = pd.DataFrame(PRINTED)
    pth = HERE / "results" / "roundR" / "figure_s4_drawn_in_band.tsv"
    pth.parent.mkdir(parents=True, exist_ok=True)
    P.to_csv(pth, sep="\t", index=False)
    print(f"  帯の中に描いた値: {pth.name}（{len(P)} 行）")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
