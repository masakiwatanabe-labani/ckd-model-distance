"""Fig 5: IRI 時系列の軌跡は単調でない。

最早時点（2 h）からの距離を対数時間軸でプロットする。
24-72 h のピーク、14 d の回帰、28 d の再乖離を注記する。

within_model_time.png（モデル内 距離 vs 時間差）は UUO がペア1個、
PodTRECK が3個しかなく関係を評価できないため本図には含めず、補足に回す。
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("27_fig5")
cfg = load_config()
RES = Path(cfg["_root"]) / cfg["paths"]["results"]
FIG = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)

COL = "#DD8452"
LAB = {"IRI2h": "2 h", "IRI4h": "4 h", "IRI24h": "24 h", "IRI48h": "48 h", "IRI72h": "72 h",
       "IRI7d": "7 d", "IRI14d": "14 d", "IRI28d": "28 d", "IRI12m": "12 mo"}


def main():
    t = pd.read_csv(RES / "iri_trajectory.csv").sort_values("days")
    t["label"] = t.state.map(LAB)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(t.days, t.distance_from_earliest, "-", color=COL, lw=1.8, zorder=3)
    ax.scatter(t.days, t.distance_from_earliest, s=62, color=COL,
               edgecolor="white", lw=1.2, zorder=4)

    off = {"2 h": (7, -12), "4 h": (-4, 10), "24 h": (-2, 11), "48 h": (0, 11),
           "72 h": (6, 8), "7 d": (7, 6), "14 d": (4, -15), "28 d": (2, 10), "12 mo": (-14, 11)}
    for _, r in t.iterrows():
        ax.annotate(r.label, (r.days, r.distance_from_earliest), textcoords="offset points",
                    xytext=off.get(r.label, (5, 7)), fontsize=8.6, zorder=5)

    peak = t[t.state.isin(["IRI24h", "IRI48h", "IRI72h"])]
    ax.axvspan(peak.days.min() * .78, peak.days.max() * 1.28, color=COL, alpha=.09, zorder=0)
    ax.text(1.75, 0.945, "acute peak  24 - 72 h", ha="center", fontsize=8.4, color="#A2582C")

    d14 = float(t.loc[t.state == "IRI14d", "distance_from_earliest"].iloc[0])
    d28 = float(t.loc[t.state == "IRI28d", "distance_from_earliest"].iloc[0])
    d12 = float(t.loc[t.state == "IRI12m", "distance_from_earliest"].iloc[0])

    ax.annotate(f"partial return toward the\nearliest state at 14 d ({d14:.2f})",
                xy=(14, d14), xytext=(2.4, 0.20), fontsize=8.4, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0), zorder=6)
    ax.annotate(f"renewed divergence\nat 28 d ({d28:.2f})",
                xy=(28, d28), xytext=(60, 0.44), fontsize=8.4, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0), zorder=6)
    ax.annotate(f"chronic state at 12 mo ({d12:.2f})",
                xy=(365, d12), xytext=(42, 0.26), fontsize=8.4, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0), zorder=6)

    ax.set_xscale("log")
    ax.set_xlabel("days after injury (log scale)", fontsize=10)
    ax.set_ylabel("distance from the earliest state (2 h)\n(1 - Spearman rho)", fontsize=10)
    ax.set_title("The IRI trajectory is not monotonic in time",
                 fontsize=11.5, loc="left")
    ax.set_ylim(-0.05, 1.0)
    ax.grid(alpha=.25)
    fig.tight_layout()
    out = FIG / "Fig5_iri_trajectory.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    log.info("saved %s", out)

    log.info("\n%s", t[["state", "days", "distance_from_earliest"]].round(4).to_string(index=False))
    append_summary("27_fig5 / Fig5", {
        "図": "results/figures/Fig5_iri_trajectory.png (300 dpi, 単一パネル)",
        "軌跡": t[["state", "days", "distance_from_earliest"]].round(4).to_dict("records"),
        "注記した点": {"acute peak 24-72 h": round(float(peak.distance_from_earliest.max()), 4),
                       "14 d": round(d14, 4), "28 d": round(d28, 4), "12 mo": round(d12, 4)},
        "除外": ("within_model_time.png はモデル内ペアが UUO 1個 / PodTRECK 3個しかなく "
                 "関係を評価できないため本図に含めず、補足資料に回す。"),
    }, cfg)


if __name__ == "__main__":
    main()
