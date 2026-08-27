"""Fig 3: 距離行列は入口ではなく時間で構造化される。

(A) 14状態の距離行列ヒートマップ（階層クラスタリング順）＋入口・経過日数サイドバー
(B) 入口一致ペア vs 不一致ペアの距離分布
(C) Mantel / 偏Mantel の rho（入口・モデル・時間）
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).parent))
from lib_figure import apply_style, savefig as _savefig  # noqa: E402
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("17_fig3")
apply_style()
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m16", Path(__file__).parent / "16_entry_model_checks.py")
m16 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m16)
m14 = m16.m14

S, DAYS, MODEL, ENTRY = m16.S, m16.DAYS, m16.MODEL, m16.ENTRY
STATES = m16.STATES

LABEL = {"PodTRECK_5D": "PodTRECK 5d", "PodTRECK_2W": "PodTRECK 14d", "PodTRECK_3W": "PodTRECK 21d",
         "UUO_2D": "UUO 2d", "UUO_8D": "UUO 8d", "IRI2h": "IRI 2h", "IRI4h": "IRI 4h",
         "IRI24h": "IRI 24h", "IRI48h": "IRI 48h", "IRI72h": "IRI 72h", "IRI7d": "IRI 7d",
         "IRI14d": "IRI 14d", "IRI28d": "IRI 28d", "IRI12m": "IRI 12mo"}
MCOL = {"PodTRECK": "#4C72B0", "IRI": "#DD8452", "UUO": "#55A868"}
ECOL = {"glomerular": "#8172B3", "tubular": "#C44E52"}


def main():
    D, E, M, T = m16.build_mats(STATES, S, MODEL, ENTRY, DAYS)
    iu = np.triu_indices(len(STATES), 1)

    # ---- 階層クラスタリング ----
    Dz = D.copy(); np.fill_diagonal(Dz, 0.0)
    Z = linkage(squareform(Dz, checks=False), method="average")
    order = dendrogram(Z, no_plot=True)["leaves"]
    ordered = [STATES[i] for i in order]
    log.info("クラスタリング順: %s", " -> ".join(LABEL[s] for s in ordered))

    fig = plt.figure(figsize=(11.4, 5.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, .82], wspace=.34)

    # ================= (A) ヒートマップ =================
    gsA = gs[0].subgridspec(2, 2, height_ratios=[.13, 1], width_ratios=[1, .045],
                            hspace=.10, wspace=.05)
    axbar = fig.add_subplot(gsA[0, 0])
    axh = fig.add_subplot(gsA[1, 0])
    axcb = fig.add_subplot(gsA[1, 1])

    Do = D[np.ix_(order, order)]
    cmap = LinearSegmentedColormap.from_list("d", ["#08306B", "#4292C6", "#DEEBF7", "#FFF5EB"])
    im = axh.imshow(Do, cmap=cmap, vmin=np.nanmin(D[iu]), vmax=np.nanmax(D[iu]))
    axh.set_xticks(range(len(ordered)))
    axh.set_xticklabels([LABEL[s] for s in ordered], rotation=90, fontsize=7.5)
    axh.set_yticks(range(len(ordered)))
    axh.set_yticklabels([LABEL[s] for s in ordered], fontsize=7.5)
    cb = fig.colorbar(im, cax=axcb); cb.set_label("distance (1 - Spearman rho)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # サイドバー: 入口 と 経過日数
    dayv = np.log10([DAYS[s] for s in ordered])
    norm = Normalize(dayv.min(), dayv.max())
    for i, s in enumerate(ordered):
        axbar.add_patch(plt.Rectangle((i, 1.05), 1, .85, color=ECOL[ENTRY[s]]))
        axbar.add_patch(plt.Rectangle((i, .05), 1, .85,
                                      color=plt.cm.viridis(norm(np.log10(DAYS[s])))))
    axbar.set_xlim(0, len(ordered)); axbar.set_ylim(0, 2)
    axbar.set_xticks([]); axbar.set_yticks([1.475, .475])
    axbar.set_yticklabels(["onset compartment", "days (log)"], fontsize=7.5)
    for sp in axbar.spines.values():
        sp.set_visible(False)
    hs = [plt.Rectangle((0, 0), 1, 1, color=ECOL[k]) for k in ECOL]
    hs += [plt.Rectangle((0, 0), 1, 1, color=plt.cm.viridis(x)) for x in (0.0, 1.0)]
    axbar.legend(hs, list(ECOL) + ["2 h", "12 mo"], fontsize=6.8, ncol=4,
                 loc="lower left", bbox_to_anchor=(-0.02, 1.35), frameon=False,
                 handlelength=1.1, columnspacing=1.0)
    axbar.set_title("A  Distance matrix (hierarchically clustered)",
                    loc="left", fontsize=10.5, pad=30, x=-.20)

    # ================= (B) 入口一致 vs 不一致 =================
    axB = fig.add_subplot(gs[1])
    same = D[iu][E[iu] == 0]; diff = D[iu][E[iu] == 1]
    parts = axB.violinplot([same, diff], showextrema=False, widths=.8)
    for pc, c in zip(parts["bodies"], ["#4C72B0", "#C44E52"]):
        pc.set_facecolor(c); pc.set_alpha(.32)
    bp = axB.boxplot([same, diff], widths=.22, patch_artist=True, showfliers=False)
    for p, c in zip(bp["boxes"], ["#4C72B0", "#C44E52"]):
        p.set_facecolor(c); p.set_alpha(.85)
    for i, (v, c) in enumerate(zip([same, diff], ["#4C72B0", "#C44E52"]), start=1):
        axB.scatter(np.random.default_rng(SEED).normal(i, .055, len(v)), v,
                    s=11, color=c, alpha=.65, zorder=3, edgecolor="none")
    axB.set_xticks([1, 2])
    axB.set_xticklabels([f"same onset\ncompartment\n(n={len(same)})", f"different\n\n(n={len(diff)})"],
                        fontsize=8.5)
    axB.set_ylabel("distance (1 - Spearman rho)", fontsize=9)
    u, p = stats.mannwhitneyu(same, diff)
    axB.set_title(f"B  Onset compartment does not separate distances\nMann-Whitney p = {p:.2f}",
                  loc="left", fontsize=10.5)
    axB.grid(axis="y", alpha=.25)

    # パネルC（Mantel 横棒）は Fig 6B に一本化したため削除。
    # Mantel 値自体は fig3_mantel_values.csv に引き続き書き出す。
    res = []
    for nm, X, part in [("onset compartment", E, None), ("model identity", M, None),
                        ("elapsed time", T, None),
                        ("onset compartment | time", E, T), ("model | time", M, T),
                        ("time | onset compartment", T, E)]:
        r, pv = m14.mantel(D, X, seed=SEED, partial=part)
        res.append({"variable": nm, "rho": r, "p": pv})
    rf = pd.DataFrame(res).iloc[::-1]

    _savefig(fig, FIG / "Fig3_distance_structure")
    log.info("saved %s", FIG / "Fig3_distance_structure.png")

    rf.iloc[::-1].to_csv(RES / "fig3_mantel_values.csv", index=False)
    pd.DataFrame({"order": range(1, len(ordered) + 1), "state": ordered,
                  "label": [LABEL[s] for s in ordered],
                  "model": [MODEL[s] for s in ordered],
                  "entry": [ENTRY[s] for s in ordered],
                  "days": [DAYS[s] for s in ordered]}).to_csv(
        RES / "fig3_cluster_order.csv", index=False)

    append_summary("17_fig3 / Fig3", {
        "図": "results/figures/Fig3_distance_structure.png (300 dpi)",
        "クラスタリング順": " -> ".join(LABEL[s] for s in ordered),
        "入口一致ペア": f"n={len(same)}, 中央値 {np.median(same):.3f}",
        "入口不一致ペア": f"n={len(diff)}, 中央値 {np.median(diff):.3f}",
        "Mann-Whitney p": round(float(p), 4),
        "Mantel値": rf.iloc[::-1].round(4).to_dict("records"),
    }, cfg)


if __name__ == "__main__":
    main()
