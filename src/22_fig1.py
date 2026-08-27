"""Fig 1: (A) 疾患状態の配置図  (B) 解析フロー模式図。"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from lib_figure import apply_style, savefig as _savefig  # noqa: E402
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("22_fig1")
apply_style()
cfg = load_config()
RES = Path(cfg["_root"]) / cfg["paths"]["results"]
FIG = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)

MOUSE = "#4C72B0"; CAT = "#C44E52"; HUMAN = "#9E9E9E"
GLOM_BG = "#EFEAF6"; TUB_BG = "#FBECEC"


def box(ax, x, y, w, h, text, fc, ec, dashed=False, fs=7.4, tc="black", lw=1.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec=ec, lw=lw,
                                linestyle=(0, (3.2, 2.0)) if dashed else "solid", zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, zorder=4, linespacing=1.45)


def panelA(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.text(-0.35, 10.15, "A", fontsize=13, fontweight="bold", va="top")

    # 区画の背景（横軸 = onset compartment）
    ax.add_patch(Rectangle((1.55, 0.35), 3.9, 8.35, fc=GLOM_BG, ec="none", zorder=0))
    ax.add_patch(Rectangle((5.55, 0.35), 4.3, 8.35, fc=TUB_BG, ec="none", zorder=0))
    ax.text(3.50, 8.92, "glomerular onset", ha="center", fontsize=9.6, fontweight="bold")
    ax.text(7.70, 8.92, "tubular onset", ha="center", fontsize=9.6, fontweight="bold")
    ax.text(0.75, 9.55, "onset compartment", ha="left", fontsize=8.4, style="italic", color="0.35")

    # 縦軸 = species
    for y, lab, col in [(6.9, "mouse", MOUSE), (3.55, "cat", CAT), (1.15, "human", HUMAN)]:
        ax.text(0.72, y, lab, ha="center", va="center", fontsize=9.6,
                fontweight="bold", color=col, rotation=90)
    for y in (5.95, 2.35):
        ax.plot([0.25, 9.85], [y, y], color="0.82", lw=.9, zorder=1)

    # --- mouse, glomerular : Pod-TRECK ---
    box(ax, 1.72, 7.55, 3.55, 0.95,
        "Pod-TRECK  (podocyte injury)\nRNA 5 d / 14 d / 21 d   n = 3 each",
        "white", MOUSE)
    box(ax, 1.72, 6.42, 3.55, 0.95,
        "Pod-TRECK proteome\nDay 14 / Day 21   n = 3 each",
        "white", MOUSE)

    # --- mouse, tubular : IRI + UUO ---
    box(ax, 5.72, 7.55, 3.95, 0.95,
        "IRI  (GSE98622)   RNA, n = 3 per state\n2 h  4 h  24 h  48 h  72 h  7 d  14 d  28 d  365 d",
        "white", MOUSE)
    box(ax, 5.72, 6.72, 3.95, 0.62,
        "IRI 180 d  (n = 4, separate platform)", "white", MOUSE, fs=7.0)
    box(ax, 5.72, 6.05, 3.95, 0.55,
        "UUO  (GSE79443)   RNA  2 d / 8 d,  n = 3", "white", MOUSE, fs=7.0)
    ax.text(9.72, 5.62, "days since injury are defined for every mouse state",
            ha="right", fontsize=6.9, style="italic", color=MOUSE)

    # --- cat : tubular, no time origin ---
    box(ax, 5.72, 3.95, 3.95, 1.35,
        "feline spontaneous CKD   (IRIS 3/4 vs control)\n"
        "RNA  cortex n = 7 vs 6 · medulla n = 5 vs 6\n"
        "protein  cortex n = 7 vs 6 · medulla n = 6 vs 6",
        "white", CAT, dashed=True, lw=1.5)
    ax.text(7.70, 3.62, "no defined time origin — outside the injury-time axis",
            ha="center", fontsize=7.4, fontweight="bold", color=CAT)
    ax.text(7.70, 3.34, "(natural onset · cross-sectional · IRIS stage = renal function, not time)",
            ha="center", fontsize=6.6, style="italic", color=CAT)

    # --- human : reference only ---
    box(ax, 1.72, 1.62, 3.55, 0.80,
        "KPMP regional proteome\nCKD vs HRT   in TI  ·  in G", "#F5F5F5", HUMAN, fs=7.0, tc="0.30")
    ax.text(3.50, 1.42, "DKD / hypertensive = glomerular onset",
            ha="center", fontsize=6.4, style="italic", color="0.5")
    box(ax, 3.55, 0.52, 4.35, 0.66,
        "ERCB tubulointerstitium (GSE104954)   DN · RPGN · HT  vs living donor",
        "#F5F5F5", HUMAN, fs=6.9, tc="0.30")
    ax.text(5.72, 0.36, "onset compartment varies by diagnosis — not assigned to either column",
            ha="center", fontsize=6.4, style="italic", color="0.5")
    ax.text(7.75, 1.90, "human cohorts:\nreference only,\noutside the\nspecies x onset\ncomparison",
            ha="center", va="center", fontsize=6.9, style="italic", color="0.45")

    # 凡例
    ax.plot([7.55, 8.15], [9.58, 9.58], color=CAT, lw=1.5, linestyle=(0, (3.2, 2.0)))
    ax.text(8.27, 9.58, "dashed = no time origin", fontsize=6.8, va="center", color="0.3")


def panelB(ax):
    ax.set_xlim(0, 10); ax.set_ylim(1.1, 10); ax.axis("off")
    ax.text(-0.35, 10.15, "B", fontsize=13, fontweight="bold", va="top")

    steps = [
        (8.55, "within-dataset contrast\nmoderated t, disease vs its own control",
         "controls never pooled across platforms or ages"),
        (6.95, "orthologue mapping to human symbol space",
         "one-to-one orthologue, else upper-case fallback; alias list queried exhaustively"),
        (5.35, "rank-based distance\n1 - Spearman rho",
         "ranks only: feline proteome values are centred per protein"),
    ]
    for y, main, sub in steps:
        box(ax, 1.35, y, 7.3, 0.88, main, "white", "0.35", fs=8.2)
        ax.text(5.0, y - 0.20, sub, ha="center", fontsize=6.8, style="italic", color="0.45")
        ax.add_patch(FancyArrowPatch((5.0, y - 0.42), (5.0, y - 0.72),
                                     arrowstyle="-|>", mutation_scale=13, color="0.45", lw=1.1))

    ax.text(5.0, 4.42, "layer separation — the two matrices are never mixed",
            ha="center", fontsize=8.4, fontweight="bold", color="0.2")

    box(ax, 0.75, 2.42, 4.0, 1.62,
        "TRANSCRIPTOME matrix\n\ncat 2  +  mouse 14\n= 16 states\n\nspecies x onset 2x2 complete",
        "#EAF0F8", MOUSE, fs=7.6)
    box(ax, 5.25, 2.42, 4.0, 1.62,
        "PROTEOME matrix\n\ncat 2 + mouse 1 + human 2\n= 5 states\n\nno mouse tubular proteome",
        "#FBECEC", CAT, fs=7.6)

    ax.text(2.75, 2.06, "decomposition of species\nvs onset performed here",
            ha="center", fontsize=6.9, style="italic", color=MOUSE)
    ax.text(7.25, 2.06, "2x2 incomplete -> used only for\ndistance to human cohorts",
            ha="center", fontsize=6.9, style="italic", color=CAT)

    ax.add_patch(Rectangle((0.55, 1.45), 8.9, 0.42, fc="#FFF6E5", ec="#D9A441", lw=1.0, zorder=3))
    ax.text(5.0, 1.66, "cross-layer comparisons are excluded from every distance matrix",
            ha="center", va="center", fontsize=7.4, color="#7A5A12", zorder=4)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(15.2, 7.4))
    panelA(axes[0]); panelB(axes[1])
    fig.tight_layout(w_pad=3.0)
    out = FIG / "Fig1_design"
    _savefig(fig, out)
    log.info("saved %s", out)
    append_summary("22_fig1 / Fig1", {
        "図": "results/figures/Fig1_design.png (300 dpi)",
        "パネルA": ("種 x onset compartment の配置。マウス各状態に傷害後日数、"
                    "ネコは破線枠＋'no defined time origin' で時間軸の外、"
                    "ヒトはグレーで本比較の外であることを明示。"),
        "パネルB": ("within-dataset contrast -> orthologue mapping -> rank-based distance -> "
                    "層分離。蛋白層で2x2が完成しないことを明記。"),
    }, cfg)


if __name__ == "__main__":
    main()
