# -*- coding: utf-8 -*-
"""C案の図をすべて解析パイプラインの出力から描く。

C案の図は本文と表の丸めた値から描き直されていた（キャプション自身が
「Values are rounded as in Table 1」「redrawn from the three-decimal values in Table S1」
「redrawn from the reported summaries」と書いている）。
ここでは結果ファイルだけを入力にし、丸め値は一切使わない。

  Figure 1  中心化と経路平均の四条件（新規。個体再標本化の分布も描く）
  Figure 2  4つの遺伝子空間での α・cos θ・ν（新規）
  Figure 3  対照の割り当てと信頼性ベンチマーク（A系 Figure 2 と同一の描画）
  Figure S1 ベンチマークのシミュレーション（同 Figure S1）
  Figure S2 経路単位の整列（同 Figure 3）
  Figure S3 経路の順位と集約（同 Figure 4）
  Figure S4 オルソログと遺伝子選択（同 Figure 5）
  Figure S5 経過時間と方向・振幅（同 Figure 6）

印字した数値は figure_printed_numbers_C.tsv に記録し、verify_numbers.py が本文と突き合わせる。
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))
import make_figures as MF  # noqa: E402
from lib_figure import EMPH, MAX_WIDTH_IN, apply_style, savefig  # noqa: E402

apply_style()
RES = HERE / "results"
PWD = RES / "revision1" / "pathway_reactome"
OUT = HERE / "manuscript_C" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIGW = MAX_WIDTH_IN
S1, S2, S3 = MF.S1, MF.S2, MF.S3
INK2, MUTED, NEUTRAL = MF.INK2, MF.MUTED, MF.NEUTRAL

COND = [("cos_raw", "Δ as\nanalysed"),
        ("cos_centred", "Δ centered\nper state"),
        ("cos_agg_centred", "pathway means,\ncentered"),
        ("cos_agg_raw", "pathway means,\nuncentered")]

PRINTED: list[tuple[str, str, str, str]] = []


def rec(figure, panel, item, value):
    PRINTED.append((figure, panel, item, str(value)))


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# --------------------------------------------------------------- Figure 1
def figure1():
    """四条件の AUC と中央値、および個体再標本化による差の分布。"""
    BA = pd.read_csv(PWD / "cos_before_after.tsv", sep="\t")
    REP = pd.read_csv(RES / "round7" / "centering_auc_animal_bootstrap_replicates.tsv", sep="\t")

    CC = pd.read_csv(RES / "roundC" / "centering_controls.tsv", sep="\t")
    CC = CC[(CC.collection.str.startswith("422")) & (CC.pairs.str.startswith("87"))]

    fig = plt.figure(figsize=(FIGW, 8.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.74, 0.95, 0.86], width_ratios=[1.0, 0.70])

    # (A) 四条件の AUC
    ax = fig.add_subplot(gs[0, :])
    xs = np.arange(4)
    aucs = [auc(BA[BA.within == 1][c], BA[BA.within == 0][c]) for c, _ in COND]
    ax.bar(xs, aucs, width=0.54, color=[NEUTRAL, NEUTRAL, S2, S1],
           edgecolor="white", linewidth=1.0, zorder=3)
    ax.axhline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
    for x, a in zip(xs, aucs):
        rec("Fig1", "A", f"AUC {COND[x][0]}", f"{a:.3f}")
        ax.text(x, a + 0.018, f"{a:.3f}", ha="center", va="bottom", zorder=5)
    ax.text(-0.44, 0.5, "chance", ha="left", va="bottom")
    ax.set_xticks(xs); ax.set_xticklabels([l for _, l in COND])
    ax.set_ylim(0, 1.06); ax.set_ylabel("AUC, same-species vs cross-species")
    ax.set_title("(A)  Centering before pathway averaging removes the separation;\n"
                 "averaging without it improves on the gene level", loc="left", pad=6)
    clean(ax)

    # (B) 各条件のクラス中央値
    ax = fig.add_subplot(gs[1, 0])
    for w, col, lab, dx in ((1, S1, "within species", -0.11), (0, S3, "cat vs mouse", +0.11)):
        d = BA[BA.within == w]
        for x, (c, _l) in enumerate(COND):
            ax.scatter([x + dx] * len(d), d[c], s=9, color=col, alpha=0.45, zorder=3,
                       label=lab if x == 0 else None)
            med = float(d[c].median())
            rec("Fig1", "B", f"median {c} [{lab}]", f"{med:.3f}")
            ax.plot([x + dx - 0.085, x + dx + 0.085], [med] * 2, color=col, linewidth=2.4,
                    zorder=5, solid_capstyle="butt")
    ax.set_xticks(xs)
    ax.set_xticklabels(["gene\nlevel", "centered", "pathway,\ncentered", "pathway,\nuncentered"])
    ax.set_ylabel(r"$\cos\theta$ between states")
    ax.legend(frameon=False, loc="lower left", handlelength=1.2, borderaxespad=0.3)
    ax.set_title("(B)  Class medians under the same four conditions", loc="left", pad=6)
    clean(ax)

    # (C) 個体再標本化での差の分布（点推定だけを見せない）
    ax = fig.add_subplot(gs[1, 1])
    drop = REP["drop"].to_numpy(float)
    ax.hist(drop, bins=34, color=NEUTRAL, edgecolor="white", linewidth=0.5, zorder=3)
    obs = aucs[3] - aucs[2]
    lo, hi = np.percentile(drop, [2.5, 97.5])
    med = float(np.median(drop))
    frac = float((drop > 0).mean())
    for v, col, lab in ((obs, S2, "observed"), (med, INK2, "resampled median")):
        ax.axvline(v, color=col, linewidth=1.6, zorder=5)
    ax.axvline(0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=4)
    rec("Fig1", "C", "observed uncentered minus centered", f"{obs:.3f}")
    rec("Fig1", "C", "resampled median", f"{med:.3f}")
    rec("Fig1", "C", "resampled 2.5th percentile", f"{lo:.3f}")
    rec("Fig1", "C", "resampled 97.5th percentile", f"{hi:.3f}")
    rec("Fig1", "C", "replicates with uncentered above centered", f"{frac * 100:.1f}%")
    ax.set_xlabel("uncentered − centered AUC\n(1,000 animal resamples)")
    ax.set_ylabel("replicates")
    ax.set_title(f"(C)  Direction holds in {frac * 100:.1f}%;\n"
                 f"size {lo:.3f} to {hi:.3f}", loc="left", pad=6)
    ax.annotate(f"observed {obs:.3f}", xy=(obs, ax.get_ylim()[1] * 0.92),
                xytext=(obs - 0.16, ax.get_ylim()[1] * 0.99), ha="right", va="top",
                arrowprops=dict(arrowstyle="->", color=S2, linewidth=0.9))
    clean(ax)

    # (D) 4つの遺伝子空間での集約 AUC。centered と uncentered を線で結ぶ。
    ax = fig.add_subplot(gs[2, :])
    SPACE_ORDER = ["Group A (2,016)", "All 1:1 orthologues (7,897)",
                   "Matched Group B (1,632)", "Matched Group A (1,632)"]
    SHORT = ["Group A\n2,016 genes", "All 1:1 orthologues\n7,897 genes",
             "Matched Group B\n1,632 genes", "Matched Group A\n1,632 genes"]
    g = {r.space: r for r in CC.itertuples()}
    xs = np.arange(len(SPACE_ORDER))
    ax.axhline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
    ax.text(len(SPACE_ORDER) - 0.34, 0.5, "chance", ha="right", va="bottom")
    for x, sp in zip(xs, SPACE_ORDER):
        r = g[sp]
        lo, hi = float(r.auc_agg_centred), float(r.auc_agg_uncentred)
        ax.plot([x, x], [lo, hi], color=MUTED, linewidth=2.0, zorder=3,
                solid_capstyle="butt")
        ax.scatter([x], [lo], s=44, marker="o", color=S2, zorder=5,
                   label="centered" if x == 0 else None)
        ax.scatter([x], [hi], s=44, marker="o", color=S1, zorder=5,
                   label="uncentered" if x == 0 else None)
        rec("Fig1", "D", f"centered aggregate AUC {sp}", f"{lo:.3f}")
        rec("Fig1", "D", f"uncentered aggregate AUC {sp}", f"{hi:.3f}")
        rec("Fig1", "D", f"difference {sp}", f"{hi - lo:.3f}")
        rec("Fig1", "D", f"pathway sets {sp}", f"{int(r.n_sets):,}")
        ax.text(x - 0.10, lo, f"{lo:.3f}", ha="right", va="center", zorder=6)
        ax.text(x - 0.10, hi, f"{hi:.3f}", ha="right", va="center", zorder=6)
        ax.annotate("", xy=(x + 0.13, hi), xytext=(x + 0.13, lo),
                    arrowprops=dict(arrowstyle="<->", color=INK2, linewidth=0.9,
                                    shrinkA=0, shrinkB=0))
        ax.text(x + 0.18, (lo + hi) / 2, f"+{hi - lo:.3f}", ha="left", va="center",
                zorder=6, **EMPH)
        ax.text(x, 0.055, f"{int(r.n_sets):,} pathway sets", ha="center",
                va="bottom", zorder=6)
    ax.set_xticks(xs); ax.set_xticklabels(SHORT)
    ax.set_xlim(-0.45, len(SPACE_ORDER) - 0.30); ax.set_ylim(0.0, 1.10)
    ax.set_ylabel("aggregate AUC over 87 pairs")
    ax.legend(frameon=False, loc="upper left", handlelength=1.0, borderaxespad=0.3,
              ncol=2, columnspacing=1.2)
    ax.set_title("(D)  The difference holds in every gene space, and its size varies;\n"
                 "the number of pathways reaching 30 genes differs between spaces",
                 loc="left", pad=6)
    clean(ax)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.72, wspace=0.26)
    print("Fig1:", savefig(fig, OUT / "Fig1_centering_and_aggregation")[0].name)


# --------------------------------------------------------------- Figure 2
SPACES = [("Group\nA", "alpha_groupA"),
          ("All\n1:1", "alpha_ortholog_all"),
          ("Match\nB", "alpha_groupB_matched2"),
          ("Match\nA", "alpha_groupA_matched")]
STATES = [("cat_CKD12", "cortex CKD1/2", S1), ("cat_med_CKD12", "medulla CKD1/2", S2)]


def figure2():
    """4つの遺伝子空間での α・cos θ・ν。結果ファイルの値をそのまま描く。

    (A) は角度が意味を持つ模式図なので縦横比を保ち、幅を抑えて中央に置く。
    凡例はパネルの外に1つだけ置く（パネル内だとデータに重なるため）。
    """
    P = {lbl: pd.read_csv(RES / d / "projection.tsv", sep="\t", index_col=0)
         for lbl, d in SPACES}
    fig = plt.figure(figsize=(FIGW, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0])

    # ---- (A) Group A の幾何。α と ν から垂直成分を復元する（丸めた表の値は使わない）
    ax = fig.add_subplot(gs[0, 1])
    pr = P["Group\nA"]
    ax.annotate("", xy=(1.0, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=INK2, linewidth=1.6))
    ax.text(0.02, -0.07, "reference: cat cortex CKD3/4", ha="left", va="top")
    for s, lab, col in STATES:
        a, nrm = float(pr.loc[s, "alpha"]), float(pr.loc[s, "nrm"])
        perp = float(np.sqrt(max(nrm ** 2 - a ** 2, 0.0)))
        ax.annotate("", xy=(a, perp), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=col, linewidth=1.8))
        ax.plot([a, a], [0, perp], color=col, linewidth=0.8, linestyle=(0, (3, 2)))
        rec("Fig2", "A", f"alpha {s}", f"{a:.3f}")
        rec("Fig2", "A", f"nu {s}", f"{nrm:.3f}")
        ax.text(a + 0.05, perp, lab.replace(" CKD1/2", "\nCKD1/2"), ha="left", va="center",
                color=col)
    ax.axhline(0, color=MUTED, linewidth=0.8, zorder=1)
    ax.axvline(1.0, color=MUTED, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xlim(-0.05, 1.55); ax.set_ylim(-0.26, 1.02)
    ax.set_xticks([0.0, 0.5, 1.0, 1.5])
    ax.set_xlabel("component along the reference (α)")
    ax.set_ylabel("component perpendicular\nto the reference")
    ax.set_title("(A)  Group A geometry", loc="left", pad=6)
    ax.set_aspect("equal", adjustable="box")
    clean(ax)

    # ---- (B–D) 4空間での α / cos θ / ν
    handles = None
    for j, (key, lab, ref, ttl) in enumerate(
            [("alpha", "projection α", 1.0, "Projection"),
             ("cos", "cos θ", None, "Direction"),
             ("nrm", "relative amplitude", 1.0, "Amplitude")]):
        ax = fig.add_subplot(gs[1, j])
        xs = np.arange(len(SPACES))
        for s, slab, col in STATES:
            ys = [float(P[lbl].loc[s, key]) for lbl, _d in SPACES]
            ax.plot(xs, ys, color=col, linewidth=1.4, marker="o", markersize=5,
                    markerfacecolor="white", markeredgewidth=1.3, zorder=4, label=slab)
            for x, y in zip(xs, ys):
                rec("Fig2", "BCD"[j], f"{key} {s} {SPACES[x][0]}".replace("\n", " "), f"{y:.3f}")
        if ref is not None:
            ax.axhline(ref, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
        ax.set_xticks(xs)
        ax.set_xticklabels([l for l, _ in SPACES])          # 回転させず水平に
        ax.set_xlim(-0.42, len(SPACES) - 0.58)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=4))  # 刻みの密度をそろえる
        ax.set_ylabel(lab)
        ax.set_title(f"({'BCD'[j]})  {ttl}", loc="left", pad=6)
        clean(ax)
        if handles is None:
            handles, labels = ax.get_legend_handles_labels()

    # ---- 共通の凡例をパネルの外に1つだけ
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               handlelength=1.6, columnspacing=2.2, bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.subplots_adjust(hspace=0.34, wspace=0.44)
    print("Fig2:", savefig(fig, OUT / "Fig2_projection_across_gene_spaces")[0].name)


def dump():
    path = PWD / "figure_printed_numbers_C.tsv"
    lines = ["figure\tpanel\titem\tprinted"] + ["\t".join(r) for r in PRINTED]
    path.write_text("\n".join(lines) + "\n")
    print("printed numbers:", path, f"({len(PRINTED)} rows)")


# A系の図をそのまま C案の番号で使う（同じ描画コード、同じ結果ファイル）
# 補足図の番号は本文の初出順に振り直してある（Part 2-2）。
REUSE = [("Fig2_controls_and_ceiling", "Fig3_controls_and_benchmarks"),
         ("Fig5_confounders", "FigS1_mapping_and_gene_selection"),
         ("Fig3_pathway_alignment", "FigS2_pathway_alignment"),
         ("Fig4_pathway_ranking_and_aggregation", "FigS3_pathway_ranking_and_aggregation"),
         ("FigS1_cosine_attenuation", "FigS4_benchmark_simulation"),
         ("Fig6_time_direction_amplitude", "FigS5_time_direction_amplitude")]


def main() -> int:
    figure1()
    figure2()
    dump()
    # A系の make_figures.py を回し、該当する図を C案の名前で置く
    MF.OUT.mkdir(parents=True, exist_ok=True)
    for fn in (MF.fig1, MF.fig2, MF.fig3, MF.fig4, MF.fig5, MF.fig6):
        fn()
    MF.dump_printed()
    import figure_s1_decomposition as FS1  # noqa: E402
    FS1.main()
    for src, dst in REUSE:
        for ext in (".pdf", ".png"):
            s = MF.OUT / (src + ext)
            if s.exists():
                shutil.copy2(s, OUT / (dst + ext))
    print(f"\nC案の図: {sorted(p.name for p in OUT.glob('*.pdf'))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
