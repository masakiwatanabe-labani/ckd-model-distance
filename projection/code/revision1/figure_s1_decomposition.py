"""Figure S1 を 3 条件の分解で描き直す。

旧 Figure S1 は非中心化コサインの影響しか見ておらず、標本数の差と混ざっていた。
新しい図は 4 条件を並べる。(i) 標本数をそろえた中心化ピアソン、(ii) 半標本基準 vs
全標本観測の中心化ピアソン、(iii) 同じ比較を非中心化コサインで、(iv) 同じ比較を
Spearman-Brown 補正済み基準で。SB 補正は半標本の信頼性を全標本に外挿する式なので、
(iv) が (i) に戻るなら、標本数の不一致は補正で解消できることになる。
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
from make_figures import S1, S2, S3, INK, INK2, FIGW, clean, savefig  # noqa: E402

OUT = HERE / "manuscript" / "figures"
T = pd.read_csv(HERE / "results" / "round5" / "ceiling_decomposition_sim.tsv", sep="\t")

PANELS = [("frac_i_pearson_matched_n", "(i)  Centred, sizes matched"),
          ("frac_ii_pearson_half_vs_full", "(ii)  Centred, half-sample benchmark"),
          ("frac_iii_cosine_half_vs_full", "(iii)  Uncentred, half-sample benchmark"),
          ("frac_v_cosine_SB_vs_full", "(iv)  Uncentred, corrected benchmark")]
PALETTE = [S1, S2, S3, INK2, "#7a5cc4"]


def main() -> int:
    d = T[(T.design == "mouse-like (3 vs 6)") & (T["shift"] == 0.26)]
    true_cos = sorted(d.true_cos.unique(), reverse=True)
    fig, axes = plt.subplots(2, 2, figsize=(FIGW, 4.6), sharey=True, sharex=True)
    for ax, (col, title) in zip(axes.ravel(), PANELS):
        for c, tc in zip(PALETTE, true_cos):
            g = d[d.true_cos == tc].sort_values("mean_r_half_cosine")
            ax.plot(g.mean_r_half_cosine, g[col], marker="o", markersize=4.0,
                    linewidth=1.5, color=c, label=f"{tc:.1f}")
        ax.axhline(0.01, color=INK, linestyle=(0, (4, 3)), linewidth=1.0)
        ax.axvspan(0.568, 0.986, color="#eceae4", zorder=0)
        ax.set_title(title, loc="left", pad=4)
        ax.set_xlabel("reliability")
        ax.set_xlim(0.1, 1.0)
        clean(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of replicates with\nobserved above benchmark")
    for ax in axes[0, :]:
        ax.set_xlabel("")
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, title="true cos θ between the responses", loc="lower center", ncol=5,
               frameon=False, handlelength=1.6, columnspacing=1.4,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.09, 1, 1), h_pad=1.0, w_pad=0.8)
    print("FigS1:", savefig(fig, OUT / "FigS1_cosine_attenuation")[0].name)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
