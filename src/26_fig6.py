"""Fig 6: モジュールレベルへの集約は種間とモデル間の差を消すが、時間軸は残す。

(A) between-model vs cross-species の距離分布を 3 レベルで同一軸に並べる
(B) 時間 / onset compartment の Mantel rho を 3 レベルで横棒
(C) モジュール数 k に対する時間 rho の安定性
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_figure import apply_style, savefig as _savefig  # noqa: E402
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("26_fig6")
apply_style()
cfg = load_config()
RES = Path(cfg["_root"]) / cfg["paths"]["results"]
FIG = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)
SEED = cfg["seed"]

BET = "#4C72B0"; CRS = "#C44E52"
LEVELS = [("gene level", "distance_distributions_16state.csv"),
          ("module\n(curated, 22)", "module_distance_pairs_curated.csv"),
          ("module\n(Hallmark, 50)", "module_distance_pairs_hallmark.csv")]


def main():
    fig = plt.figure(figsize=(16.0, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 1.22, 1.0], wspace=.42)
    axA, axB, axC = (fig.add_subplot(gs[i]) for i in range(3))

    # ---------------- (A) ----------------
    stats_rows, pos, ticks, tlab = [], [], [], []
    rng = np.random.default_rng(SEED)
    for i, (lab, f) in enumerate(LEVELS):
        t = pd.read_csv(RES / f)
        b = t[t.kind == "between_model"].distance.dropna().values
        c = t[t.kind == "cross_species"].distance.dropna().values
        u, p = stats.mannwhitneyu(c, b, alternative="greater")
        stats_rows.append({"level": lab.replace("\n", " "), "between_median": np.median(b),
                           "cross_median": np.median(c), "ratio": np.median(c) / np.median(b),
                           "mw_p": p, "n_between": len(b), "n_cross": len(c)})
        base = i * 3.0
        for k, (vals, col) in enumerate([(b, BET), (c, CRS)]):
            x = base + k
            pos.append(x)
            vp = axA.violinplot([vals], positions=[x], showextrema=False, widths=.78)
            vp["bodies"][0].set_facecolor(col); vp["bodies"][0].set_alpha(.30)
            bp = axA.boxplot([vals], positions=[x], widths=.20, patch_artist=True,
                             showfliers=False, manage_ticks=False)
            bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(.9)
            bp["medians"][0].set_color("black")
            axA.scatter(rng.normal(x, .045, len(vals)), vals, s=8, color=col,
                        alpha=.55, zorder=3, edgecolor="none")
        ticks.append(base + .5); tlab.append(lab)
        axA.text(base + .5, 1.155, f"p = {p:.3f}" if p >= 1e-3 else f"p = {p:.0e}",
                 ha="center", fontsize=8.2, color="#C44E52" if p < .05 else "0.45",
                 fontweight="bold" if p < .05 else "normal")
    axA.set_xticks(ticks); axA.set_xticklabels(tlab, fontsize=8.6)
    axA.set_ylabel("distance (1 - Spearman rho)", fontsize=9.5)
    axA.set_ylim(0, 1.30)
    axA.set_title("A  Aggregation collapses the species / model-choice gap",
                  loc="left", fontsize=10.5)
    axA.grid(axis="y", alpha=.22)
    hs = [plt.Rectangle((0, 0), 1, 1, color=BET, alpha=.75),
          plt.Rectangle((0, 0), 1, 1, color=CRS, alpha=.75)]
    axA.legend(hs, ["between-model (mouse)", "cross-species (cat x mouse)"],
               fontsize=7.6, loc="upper center", bbox_to_anchor=(.5, 1.005),
               frameon=False, ncol=2, handlelength=1.2)

    # ---------------- (B) ----------------
    import json
    blob = json.load(open(RES / "summary.json", encoding="utf-8"))
    mod = blob["24_module_level / モジュールレベルの距離と Mantel"]
    gene = blob["14_time_axis_full / GSE98622全時点での時間軸 vs 入口"][
        "B: IRI6mo除く(14状態, 全て同一プラットフォーム内)"]

    series = [
        ("gene level", gene["時間差 |Δlog10(日数)|"]["rho"], gene["時間差 |Δlog10(日数)|"]["p"],
         gene["入口不一致"]["rho"], gene["入口不一致"]["p"]),
        ("module (curated, 22)", mod["curated (config/modules.yaml)"]["mantel"]["elapsed time (14 mouse states)"]["rho"],
         mod["curated (config/modules.yaml)"]["mantel"]["elapsed time (14 mouse states)"]["p"],
         mod["curated (config/modules.yaml)"]["mantel"]["onset compartment (14 mouse)"]["rho"],
         mod["curated (config/modules.yaml)"]["mantel"]["onset compartment (14 mouse)"]["p"]),
        ("module (Hallmark, 50)", mod["MSigDB Hallmark"]["mantel"]["elapsed time (14 mouse states)"]["rho"],
         mod["MSigDB Hallmark"]["mantel"]["elapsed time (14 mouse states)"]["p"],
         mod["MSigDB Hallmark"]["mantel"]["onset compartment (14 mouse)"]["rho"],
         mod["MSigDB Hallmark"]["mantel"]["onset compartment (14 mouse)"]["p"]),
    ]
    bars = []
    for lab, rt, pt, re_, pe in series:
        bars.append((f"{lab}\n  elapsed time", rt, pt))
        bars.append((f"{lab}\n  onset compartment", re_, pe))
    bars = bars[::-1]
    y = np.arange(len(bars))
    axB.barh(y, [b[1] for b in bars], height=.62, zorder=3,
             color=["#C44E52" if b[2] < .05 else "#C9C9C9" for b in bars])
    axB.set_yticks(y)
    axB.set_yticklabels([b[0] for b in bars], fontsize=7.6, linespacing=1.5)
    for i, b in enumerate(bars):
        axB.text(b[1] + (.018 if b[1] >= 0 else -.018), i, f"{b[1]:+.3f}  p={b[2]:.4f}",
                 va="center", ha="left" if b[1] >= 0 else "right", fontsize=7.4)
    axB.axvline(0, color="0.3", lw=.9)
    axB.set_xlim(-.62, .95); axB.set_xlabel("Mantel rho", fontsize=9.5)
    axB.set_title("B  Time survives aggregation, onset compartment does not\n(red: p < 0.05)",
                  loc="left", fontsize=10.5)
    axB.grid(axis="x", alpha=.22)

    # ---------------- (C) ----------------
    curve = pd.read_csv(RES / "module_power_curve.csv")
    for tag, col, mk in [("curated", "#4C72B0", "o"), ("hallmark", "#DD8452", "s")]:
        d = curve[curve["定義"] == tag].drop_duplicates("k_modules").sort_values("k_modules")
        axC.plot(d.k_modules, d.time_rho_median, mk + "-", color=col, label=f"{tag}", zorder=4)
        axC.fill_between(d.k_modules, d["time_rho_p2.5"], d["time_rho_p97.5"],
                         color=col, alpha=.16, zorder=2)
    axC.axhline(0, color="0.35", lw=1.0, ls="--", zorder=3)
    axC.set_xlabel("number of modules used (k)", fontsize=9.5)
    axC.set_ylabel("Mantel rho (elapsed time)", fontsize=9.5)
    axC.set_title("C  Estimate is unstable below ~15 modules\n(median and 95% interval)",
                  loc="left", fontsize=10.5)
    axC.legend(fontsize=8, frameon=False, loc="lower right")
    axC.grid(alpha=.22)
    axC.annotate("interval crosses zero at k = 5", xy=(5.4, -0.08), xytext=(13, 0.10),
                 fontsize=7.4, color="0.3",
                 arrowprops=dict(arrowstyle="->", color="0.45", lw=.9))
    axC.set_ylim(-0.32, 0.80)

    _savefig(fig, FIG / "Fig6_module_level")
    log.info("saved %s", FIG / "Fig6_module_level.png")

    sdf = pd.DataFrame(stats_rows)
    sdf.to_csv(RES / "fig6_panelA_values.csv", index=False)
    log.info("\n%s", sdf.round(4).to_string(index=False))
    append_summary("26_fig6 / Fig6", {
        "図": "results/figures/Fig6_module_level.png (300 dpi)",
        "パネルA": sdf.round(4).to_dict("records"),
        "パネルB": [{"level": s[0], "time_rho": s[1], "time_p": s[2],
                     "onset_rho": s[3], "onset_p": s[4]} for s in series],
        "パネルC": "results/module_power_curve.csv",
    }, cfg)


if __name__ == "__main__":
    main()
