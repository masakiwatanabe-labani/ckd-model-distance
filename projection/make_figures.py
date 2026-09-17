"""本文図の作成（Fig 1-6）。

  Fig 1  §2.1  幾何模式 + alpha / cos / nrm の分解
  Fig 2  §2.3  時間差 vs 方向 / 振幅
  Fig 3  §2.4  共有対照バイアスの降格 + 観測 cos と天井
  Fig 4  §2.5  経路サイズのノイズ床 + 状態 x 経路のヒートマップ
  Fig 5  §2.5  経路間の順位入れ替わり + 経路単位 vs 集約 + サイズ非依存
  Fig 6  §2.6  3つの交絡候補

級数・書体・線幅・文字色はこのファイルでは決めない。すべて src/lib_figure.py の
apply_style() が rcParams で与える。図ごとに fontsize= や文字の color= を書かないこと
（7 点でばらつくのを防ぐため）。太字は EMPH を付けた数値ラベルだけに使う。

凡例の様式は2種類に統一する:
  - 複数パネルで共通の符号化 → 図下部の fig.legend
  - そのパネル固有の凡例       → ax.legend

色は dataviz スキルの reference palette のスロット 1-3 をそのまま使う
（#2a78d6 / #eb6834 / #1baf7a）。全図で形状・位置・直接ラベルが情報を担い、
色だけに意味を持たせないためグレースケールでも読める。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
from lib_figure import EMPH, MAX_WIDTH_IN, apply_style, measure, savefig  # noqa: E402

apply_style()
OUT = HERE / "manuscript" / "figures"
PWD = HERE / "results" / "revision1" / "pathway_reactome"   # Reactome を含む 422 セット
OUT.mkdir(parents=True, exist_ok=True)
RES = HERE / "results"

# 図に印字した数値をそのまま書き出す。本文・キャプションとの照合は verify_numbers.py が行う。
PRINTED: list[tuple[str, str, str, str]] = []


def col_name(c: str) -> str:
    return {"cos_raw": "gene level", "cos_centred": "centred",
            "cos_agg_centred": "pathway means of centred",
            "cos_agg_raw": "pathway means of uncentred"}.get(c, c)


def rec(figure: str, panel: str, item: str, value) -> None:
    """図のテキストとして実際に描いた値を記録する。value は描いた文字列そのもの。"""
    PRINTED.append((figure, panel, item, str(value)))


def dump_printed(path=None) -> Path:
    path = Path(path) if path else PWD / "figure_printed_numbers.tsv"
    lines = ["figure\tpanel\titem\tprinted"]
    lines += ["\t".join(r) for r in PRINTED]
    path.write_text("\n".join(lines) + "\n")
    print("printed numbers:", path, f"({len(PRINTED)} rows)")
    return path


S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, NEUTRAL = "#000000", "#52514e", "#b8b7b2", "#d9d8d4"
# INK は文字色として残っている箇所の保険。図中の文字色は lib_figure の rcParams が黒に
# 統一しており、INK2 と MUTED は線と縁にしか使わない（薄いグレーの文字を作らないため）。
FIGW = MAX_WIDTH_IN   # 全図の幅。MDPI 本文幅 6.3 インチを上限とする

ORDER = [
    ("cat_CKD12", "Cat cortex CKD1/2", "Cat cortex"),
    ("cat_CKD34", "Cat cortex CKD3/4 (reference)", "Cat cortex"),
    ("cat_med_CKD12", "Cat medulla CKD1/2", "Cat medulla"),
    ("cat_med_CKD34", "Cat medulla CKD3/4", "Cat medulla"),
    ("mouse_5D", "Pod-TRECK D5", "Pod-TRECK"),
    ("mouse_2W", "Pod-TRECK D14", "Pod-TRECK"),
    ("mouse_3W", "Pod-TRECK D21", "Pod-TRECK"),
    ("IRI_2h", "IRI 2 h", "IRI"), ("IRI_4h", "IRI 4 h", "IRI"),
    ("IRI_24h", "IRI 24 h", "IRI"), ("IRI_48h", "IRI 48 h", "IRI"),
    ("IRI_72h", "IRI 72 h", "IRI"), ("IRI_7d", "IRI 7 d", "IRI"),
    ("IRI_14d", "IRI 14 d", "IRI"), ("IRI_28d", "IRI 28 d", "IRI"),
    ("IRI_12mo", "IRI 12 mo", "IRI"),
]
LAB = {k: l for k, l, _ in ORDER}
HILITE = "cat_med_CKD12"


def clean(ax, grid="both"):
    """軸まわりの共通処理。線幅・色・級数は rcParams に任せ、ここでは形だけ決める。"""
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if grid:
        ax.grid(axis=grid if grid != "both" else "both", zorder=0)
        ax.set_axisbelow(True)


# ============================================================ Fig 1
def fig1():
    P = pd.read_csv(RES / "alpha_groupA" / "projection.tsv", sep="\t", index_col=0)
    fig = plt.figure(figsize=(FIGW, 7.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.86, 1.6])

    # ---- (A) 順序逆転（主張1）。基準超え（主張2）は同じ図の注記で示す。
    ax = fig.add_subplot(gs[0, :2])
    COMP = "cat_CKD12"
    geom = {}
    for key, col in [(COMP, S3), (HILITE, S2)]:
        c, n = float(P.loc[key, "cos"]), float(P.loc[key, "nrm"])
        geom[key] = (c, n, n * c, n * np.sqrt(1 - c ** 2), col)
    ax.annotate("", xy=(1, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=2.2,
                                shrinkA=0, shrinkB=0), zorder=5)
    for key, (c, n, vx, vy, col) in geom.items():
        ax.annotate("", xy=(vx, vy), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=col, linewidth=2.2,
                                    shrinkA=0, shrinkB=0), zorder=5)
        ax.plot([vx, vx], [0, vy], color=col, linestyle=(0, (3, 3)), linewidth=1.1, zorder=4)
        ax.plot([vx], [0], marker="o", markersize=4.5, color=col, zorder=6)
    # 角度の弧（皮質のほうが基準に近い）
    for key, r in [(COMP, 0.46), (HILITE, 0.70)]:
        c, n, vx, vy, col = geom[key]
        ax.add_patch(Arc((0, 0), r, r, theta1=0, theta2=np.degrees(np.arccos(c)),
                         color=col, linewidth=1.1, zorder=4))
    ax.text(1.06, 0.04, "a  (reference),  " r"$\|a\|$ = 1", ha="left", va="bottom", zorder=7,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    for key, lab, dy, va in [(COMP, "cat cortex CKD1/2", 0.10, "bottom"),
                             (HILITE, "cat medulla CKD1/2", 0.0, "center")]:
        c, _n, vx, vy, col = geom[key]
        ax.text(vx + 0.06, vy + dy, lab + "\n" + r"$\cos\theta$ = " + f"{c:.3f}",
                ha="left", va=va, linespacing=1.3, zorder=7,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    # 射影の足を下の帯で比べる
    for key, y0 in [(COMP, -0.22), (HILITE, -0.40)]:
        c, n, vx, vy, col = geom[key]
        ax.plot([0, vx], [y0, y0], color=col, linewidth=1.6, zorder=4, solid_capstyle="butt")
        ax.plot([vx, vx], [y0 - 0.04, y0 + 0.04], color=col, linewidth=1.6, zorder=4)
        ax.text(vx + 0.04, y0, r"$\alpha$ = " + f"{c * n:.3f}", ha="left", va="center")
    ax.plot([0, 0], [-0.44, -0.18], color=INK2, linewidth=1.0, zorder=4)
    ax.axvline(1.0, color=MUTED, linestyle=(0, (4, 3)), linewidth=1.0, zorder=2)
    ax.set_xlim(-0.08, 1.92); ax.set_ylim(-0.52, 1.18)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(A)  The better-aligned state receives the smaller projection", loc="left", pad=4)
    axt = fig.add_subplot(gs[0, 2]); axt.axis("off")
    axt.text(0.0, 1.00, r"$\alpha = \nu \cos\theta$", ha="left", va="top", **EMPH)
    axt.text(0.0, 0.86,
             "Cat medulla CKD1/2 points\n"
             r"farther off-axis ($\cos\theta$ 0.790" "\n"
             "vs 0.907) and still scores\n"
             "higher: it is 1.291 times the\n"
             "length of a. This ordering\n"
             "holds in all three gene spaces.", ha="left", va="top", linespacing=1.45)
    axt.text(0.0, 0.40,
             "Here it also passes the\n"
             r"reference itself ($\alpha$ = 1.020 >" "\n"
             "1): such a score has no upper\n"
             "reference point. That crossing\n"
             "is specific to this gene set.", ha="left", va="top", linespacing=1.45)

    # ---- (B)(C)(D) 分解
    ypos, y, last = {}, 0.0, None
    for key, _l, grp in ORDER:
        if last is not None and grp != last:
            y += 0.85
        ypos[key] = y; y += 1.0; last = grp
    ymax = y
    axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    panels = [("alpha", r"(B)  $\alpha$", 1.78), ("cos", r"(C)  $\cos\theta$", 1.80),
              ("nrm", r"(D)  $\|v\|\,/\,\|a\|$", 2.15)]
    for ax, (col, title, xmax) in zip(axes, panels):
        for key, _l, _g in ORDER:
            v = P.loc[key, col]
            ax.barh(ymax - ypos[key], v, height=0.72,
                    color=S2 if key == HILITE else (S3 if key == COMP else S1),
                    edgecolor="white", linewidth=0.8, zorder=3)
            ax.plot([P.loc[key, f"{col}_lo95"], P.loc[key, f"{col}_hi95"]],
                    [ymax - ypos[key]] * 2, color=INK2, linewidth=1.1,
                    solid_capstyle="butt", zorder=4)
        ax.axvline(1.0, color=MUTED, linestyle=(0, (4, 3)), linewidth=1.0, zorder=2)
        ax.set_title(title, loc="left", pad=6)
        clean(ax, grid="x")
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(0, xmax)
        ax.set_ylim(ymax - y + 0.4, ymax + 0.7)
    axes[0].set_yticks([ymax - ypos[k] for k, _l, _g in ORDER])
    axes[0].set_yticklabels([l for _k, l, _g in ORDER])
    for ax in axes[1:]:
        ax.set_yticks([])
    hy = ymax - ypos[HILITE]
    for ax, col, txt in [(axes[0], "alpha", r"$\alpha$ = 1.020"),
                         (axes[1], "cos", r"$\cos\theta$ = 0.790"),
                         (axes[2], "nrm", "1.291 $\\times$\nreference")]:
        ax.text(P.loc[HILITE, col] + 0.04, hy, txt,
                ha="left", va="center", zorder=6, linespacing=1.3, **EMPH)
    cy = ymax - ypos[COMP]
    axes[1].text(P.loc[COMP, "cos"] + 0.04, cy, r"$\cos\theta$ = 0.907", ha="left", va="center",
                 zorder=6, **EMPH)
    axes[0].text(P.loc[COMP, "alpha"] + 0.04, cy, r"$\alpha$ = 0.544", ha="left", va="center",
                 zorder=6, **EMPH)
    axes[1].annotate("", xy=(1.12, cy - 0.35), xytext=(1.12, hy + 0.35),
                     arrowprops=dict(arrowstyle="<->", color=INK2, linewidth=1.0), zorder=5)
    axes[1].text(1.17, (cy + hy) / 2, "farther\nfrom the\nreference", ha="left", va="center",
                 zorder=6, linespacing=1.3)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=S1, edgecolor="white"),
                        plt.Rectangle((0, 0), 1, 1, facecolor=S3, edgecolor="white"),
                        plt.Rectangle((0, 0), 1, 1, facecolor=S2, edgecolor="white"),
                        Line2D([0], [0], color=MUTED, linestyle=(0, (4, 3)))],
               labels=["Disease state", "Cat cortex CKD1/2 (better aligned, lower α)",
                       "Cat medulla CKD1/2 (worse aligned, higher α)",
                       "Value of the reference state (= 1)"],
               loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.045), handlelength=1.6, columnspacing=1.4)
    fig.tight_layout(rect=(0, 0.062, 1, 1), h_pad=0.6, w_pad=0.4)
    fig.subplots_adjust(hspace=0.16, wspace=0.10)
    print("Fig1:", savefig(fig, OUT / "Fig1_alpha_decomposition")[0].name)
    plt.close(fig)


# ============================================================ Fig 2
def fig2():
    P = pd.read_csv(RES / "mantel_s1" / "distance_pairs.tsv", sep="\t")
    tm = pd.read_csv(HERE / "time_map.tsv", sep="\t").set_index("state")["hours"]
    M = pd.read_csv(RES / "mantel_s1" / "mantel_s1.tsv", sep="\t")
    M = M[M.subset.str.startswith("12")]

    def grp(a, b):
        ga = "IRI" if a.startswith("IRI") else "PodTRECK"
        gb = "IRI" if b.startswith("IRI") else "PodTRECK"
        return "IRI x IRI" if ga == gb == "IRI" else (
            "Pod-TRECK x Pod-TRECK" if ga == gb else "IRI x Pod-TRECK")

    d = P[P.a.isin(tm.dropna().index) & P.b.isin(tm.dropna().index)].copy()
    d["dt"] = [abs(np.log10(tm[a]) - np.log10(tm[b])) for a, b in zip(d.a, d.b)]
    d["pt"] = [grp(a, b) for a, b in zip(d.a, d.b)]
    style = {"IRI x IRI": (S1, "o"), "IRI x Pod-TRECK": (S2, "s"),
             "Pod-TRECK x Pod-TRECK": (S3, "^")}
    fig, axes = plt.subplots(1, 2, figsize=(FIGW, 3.45))
    for ax, col, title in [
            (axes[0], "D_cos (1-cos)", r"(A)  Direction   $1-\cos\theta$"),
            (axes[1], "D_nrm (|log2 norm ratio|)",
             r"(B)  Amplitude   $|\mathrm{log2}\,(\|\Delta i\|\,/\,\|\Delta j\|)|$")]:
        for pt, (c, mk) in style.items():
            s = d[d.pt == pt]
            ax.scatter(s.dt, s[col], s=36, marker=mk, facecolor=c, edgecolor="white",
                       linewidth=0.8, zorder=4, label=pt)
        r = M[(M.dissimilarity == col) & (M.predictor == "elapsed time |Δlog10 h|")]
        ax.text(0.03, 0.97, f"Mantel $\\rho$ = {float(r.mantel_rho.iloc[0]):+.3f}\n"
                f"$p$ = {float(r.p_perm.iloc[0]):.4f}", transform=ax.transAxes, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=MUTED, linewidth=0.8))
        ax.set_title(title, loc="left", pad=6)
        ax.set_xlabel("|Δ log10 hours| between states")
        clean(ax)
    axes[0].set_ylabel("Dissimilarity between states")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.06), handletextpad=0.4, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    print("Fig6:", savefig(fig, OUT / "Fig6_time_direction_amplitude")[0].name)
    plt.close(fig)


# ============================================================ Fig 3
def fig3():
    SC = pd.read_csv(RES / "control_analysis" / "split_control.tsv", sep="\t")
    DIST = pd.read_csv(RES / "control_analysis" / "cos_distributions.tsv",
                       sep="\t").set_index("set")
    P = pd.read_csv(RES / "reliability" / "pairs_ceilings.tsv", sep="\t")
    test = SC[SC.pair.str.contains("cat_CKD12")].iloc[0]
    iri = SC[~SC.pair.str.contains("cat_CKD12")]

    # 主たる基準は Spearman-Brown 補正済みなので、両方を持つ表から SB 列を使う。
    U = pd.read_csv(RES / "revision1" / "pair_uncertainty_both_ceilings.tsv", sep="\t")
    fig = plt.figure(figsize=(FIGW, 8.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.92])

    # ---- (A) 共有対照バイアスの降格
    ax = fig.add_subplot(gs[0, :])
    for _, r in iri.iterrows():
        ax.plot([0, 1], [r.cos_shared_control, r.cos_split_median], color=NEUTRAL,
                linewidth=0.9, zorder=2)
        ax.scatter([0, 1], [r.cos_shared_control, r.cos_split_median], s=16,
                   marker="o", facecolor="white", edgecolor=INK2, linewidth=0.8, zorder=3)
    ax.plot([0, 1], [test.cos_shared_control, test.cos_split_median], color=S2,
            linewidth=2.4, zorder=5)
    ax.scatter([0, 1], [test.cos_shared_control, test.cos_split_median], s=64,
               marker="D", facecolor=S2, edgecolor="white", linewidth=1.0, zorder=6)
    ax.plot([1, 1], [test.cos_split_lo, test.cos_split_hi], color=S2, linewidth=1.6,
            zorder=5, solid_capstyle="butt")
    pct = float((iri.cos_split_median < test.cos_split_median).mean() * 100)
    ax.annotate(f"cat cortex CKD1/2 × CKD3/4\n{test.cos_shared_control:.3f} to "
                f"{test.cos_split_median:.3f}\n({pct:.0f}th percentile of\nthe IRI pairs)",
                xy=(1.02, test.cos_split_median), xytext=(1.12, 0.60), ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color=S2, linewidth=1.0))
    ax.text(-0.32, 1.13, "IRI time-point pairs\n(same treatment applied)", ha="left", va="top")
    ax.text(0.62, 1.15, f"pairs sharing controls: median {DIST.loc['shared_control_pairs','median']:.3f}"
            f"   ·   not sharing: {DIST.loc['no_shared_control_pairs','median']:.3f}",
            ha="center", va="bottom")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["computed against\nshared controls",
                        "computed against\ndisjoint control halves"])
    ax.set_xlim(-0.34, 1.78); ax.set_ylim(0.05, 1.24)
    ax.set_ylabel(r"$\cos\theta$ between states")
    ax.set_title("(A)  Sharing controls inflates similarity; the same treatment\n"
                 "applied to both removes the outlier", loc="left", pad=6)
    clean(ax, grid="y")

    # ---- (B)(C) 観測 cos vs 天井
    style = {"within_dataset": (S1, "Within dataset", "o"),
             "same_species_diff_dataset": (S2, "Same species, different dataset", "s"),
             "cross_species": (S3, "Cross-species (cat vs mouse)", "^")}
    cs = P[P["class"] == "cross_species"]
    X0 = 0.55
    # 主たる基準は Spearman-Brown 補正済み（較正されている。Section 4.10）なので先に置く。
    for k, (ccol, ctitle) in enumerate([
            ("ceiling_SB", "(B)  Spearman–Brown benchmark"),
            ("ceiling_raw", "(C)  Uncorrected benchmark")]):
        ax = fig.add_subplot(gs[1, k])
        ax.fill_between([X0, 1.02], [X0, 1.02], [1.02, 1.02], color="#f1f1f0", zorder=0)
        ax.plot([X0, 1.02], [X0, 1.02], color=INK2, linewidth=1.2, zorder=3)
        for cls, (col, lab, mk) in style.items():
            d = P[P["class"] == cls]
            ax.scatter(d[ccol], d.cos, s=30, marker=mk, facecolor=col, edgecolor="white",
                       linewidth=0.8, zorder=4, label=lab)
        gap = float((P.loc[P["class"] == "cross_species", ccol] - cs.cos).min())
        xg = float(cs.cos.max() + gap)
        n_above = int((P.cos > P[ccol]).sum())
        rec("Fig2", "BC"[k], f"smallest cross-species gap, {ccol}", f"{gap:.3f}")
        rec("Fig2", "BC"[k], f"pairs above the benchmark, {ccol}", n_above)
        ax.text(0.558, 1.010, "Above the estimated\nbenchmark under the\nindependent-error model",
                ha="left", va="top", linespacing=1.3, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="#f1f1f0", edgecolor="none",
                          alpha=0.9))
        ax.axhline(cs.cos.max(), color=S3, linestyle=(0, (4, 3)), linewidth=1.0, zorder=2)
        ax.annotate("", xy=(xg, xg), xytext=(xg, cs.cos.max()),
                    arrowprops=dict(arrowstyle="<->", color=INK, linewidth=1.0), zorder=6)
        ax.text(xg - 0.012, 0.5 * (xg + cs.cos.max()), f"{gap:.3f}", ha="right", va="center",
                zorder=6, bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                                    edgecolor="none", alpha=0.85))
        if ccol == "ceiling_SB":
            for _, r in P[P.cos > P[ccol]].iterrows():
                ax.scatter([r[ccol]], [r.cos], s=48, marker="o", facecolor="none",
                           edgecolor=INK, linewidth=1.2, zorder=5)
            ax.text(0.558, 0.035, f"{n_above} pairs above the benchmark (circled);\n"
                    "each shares a cohort or its controls.\nNo cross-species pair does.\n"
                    "dashed: highest cross-species " r"$\cos\theta$ = 0.548",
                    ha="left", va="bottom", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="none", alpha=0.88))
        else:
            ax.text(0.558, 0.035, f"{n_above} pairs above the benchmark;\n"
                    "this form sits below the\ncalibrated one (Section 4.10)",
                    ha="left", va="bottom", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="none", alpha=0.88))
        ax.set_title(ctitle, loc="left", pad=6)
        ax.set_xlabel("Reliability-derived benchmark for the pair")
        ax.set_xlim(X0, 1.02); ax.set_ylim(0, 1.02)
        clean(ax)
        if k == 0:
            ax.set_ylabel(r"Observed $\cos\theta$")
        else:
            ax.set_yticklabels([])
        if k == 0:
            h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.028), handletextpad=0.4, columnspacing=1.6)
    # ---- (D) 個体ブートストラップによる「天井 − 観測」の対応のある区間
    ax = fig.add_subplot(gs[2, :])
    sub = U[(U.cls == "cross_species") | (~U.shared_control)].copy()
    sub["cat"] = sub.a.where(sub.a.str.startswith("cat"), sub.b)
    REL = {"cat_CKD34", "cat_med_CKD12"}          # 分割半信頼性 0.841 / 0.851
    sub["grp"] = np.where(sub.cls != "cross_species", 0,
                          np.where(sub.cat.isin(REL), 1, 2))
    sub = sub.sort_values(["grp", "gap_sb_median"]).reset_index(drop=True)
    style3 = {"within_dataset": (S1, "o"), "same_species_diff_dataset": (S2, "s"),
              "cross_species": (S3, "^")}
    for i, r in sub.iterrows():
        col, mk = style3[r.cls]
        ax.plot([i, i], [r.gap_sb_lo, r.gap_sb_hi], color=col, linewidth=1.1, zorder=3,
                solid_capstyle="butt")
        ax.plot([i], [r.gap_sb_median], marker=mk, markersize=3.4, color=col, zorder=4)
    ax.axhline(0, color=INK, linewidth=1.2, zorder=5)
    sizes = [int((sub.grp == g).sum()) for g in (0, 1, 2)]
    edges = np.cumsum(sizes)
    for e in edges[:-1]:
        ax.axvline(e - 0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
    top = ax.get_ylim()[1]
    labs = [(f"same species,\nno shared controls ({sizes[0]})", S1),
            (f"cat vs mouse, reliable\nfeline state ({sizes[1]})", S3),
            (f"cat vs mouse, less reliable\nfeline state ({sizes[2]})", S3)]
    starts = np.concatenate([[0], edges[:-1]])
    for (lab, col), s, n in zip(labs, starts, sizes):
        ax.text(s + n / 2, top * 0.99, lab, ha="center", va="top",
                linespacing=1.3)
    for s, n_ in zip(starts, sizes):                 # 群ごとの gap 中央値を水平線で重ねる
        med = float(sub.iloc[s:s + n_].gap_sb_median.median())
        ax.plot([s - 0.5, s + n_ - 0.5], [med, med], color=INK, linewidth=1.6,
                zorder=6, solid_capstyle="butt")
        rec("Fig2", "D", f"group median gap, group starting at {s}", f"{med:.3f}")
        ax.text(s + n_ - 0.5, med, f" {med:.3f}",
                ha="left", va="center", zorder=7)
    ax.set_ylim(ax.get_ylim()[0], top * 1.22)
    ax.set_ylabel(r"benchmark $-$ observed $\cos\theta$")
    # 判定は撤回したので、ペアごとの可否を含意する文言を図から外す。
    # 群の中央値（黒い水平線）が報告される量で、縦棒は幅を示すだけ。
    ax.set_xlabel("Black lines are the group medians, which are the reported quantity;\n"
                  "the vertical ranges convey spread and are not read pair by pair",
                  linespacing=1.4)
    ax.set_xticks([])
    ax.set_title("(D)  Resampling animals, not genes: the paired difference and its\n"
                 "2.5th to 97.5th percentile range", loc="left", pad=6)
    clean(ax, grid="y")

    fig.tight_layout(rect=(0, 0.028, 1, 1), h_pad=0.8, w_pad=0.6)
    fig.subplots_adjust(hspace=0.46, wspace=0.22)
    print("Fig2:", savefig(fig, OUT / "Fig2_controls_and_ceiling")[0].name)
    plt.close(fig)


# ============================================================ Fig 4
def fig4():
    R = pd.read_csv(PWD / "pathway_cos.tsv", sep="\t")
    N = pd.read_csv(PWD / "random_set_null.tsv", sep="\t")
    fig = plt.figure(figsize=(FIGW, 5.65))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.85, 1.45])

    ax = fig.add_subplot(gs[0])
    b = N.groupby("n_genes").apply(lambda t: (t.hi95 - t.lo95).median())
    ax.plot(b.index, b.values, color=S1, linewidth=2, marker="o", markersize=5,
            markeredgecolor="white", zorder=4)
    ax.axvline(50, color=INK2, linestyle=(0, (4, 3)), linewidth=1.0, zorder=3)
    # 集合の同一性はコレクションと名前の組。Apoptosis / DNA Repair / Fatty Acid Metabolism の
    # 3 名は 2 コレクションに別メンバーで存在するので、名前だけで畳むと 422 が 419 になる。
    R = R.assign(key=R.collection + "|" + R.pathway)
    sizes = R.groupby("key").n_genes.first()
    ax2 = ax.twinx()
    ax2.hist(sizes, bins=np.arange(30, 140, 6), color=NEUTRAL, zorder=1)
    ax2.set_ylabel("pathways")
    ax2.tick_params( length=3)
    ax2.spines["top"].set_visible(False)
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    n_below, n_sets = int((sizes < 50).sum()), len(sizes)
    ax.text(54, 0.37, f"n = 50\n{n_below} of {n_sets} pathways fall below", ha="left", va="top")
    rec("Fig3", "A", "pathway sets in the size histogram", n_sets)
    rec("Fig3", "A", "sets below 50 Group A genes", n_below)
    ax.set_xlabel("Genes per pathway (Group A space)")
    ax.set_ylabel(r"95% width of $\cos\theta$" "\nfor random gene sets")
    ax.set_title("(A)  How much a pathway score can move by chance alone", loc="left", pad=6)
    ax.set_xlim(25, 135)
    clean(ax, grid=None)

    ax = fig.add_subplot(gs[1])
    # Section 4.13: パネルB は従来 3 コレクション（GO BP / KEGG / Hallmark）のうち
    # 50 遺伝子以上の 74 セットに限る。Reactome を足すと 178 列になり判読できない。
    big = R[(R.n_genes >= 50) & (R.collection != "Reactome")]
    piv = big.pivot_table(index="state", columns="key", values="cos")
    order = [k for k, _l, _g in ORDER if k in piv.index]
    piv = piv.loc[order]
    mo = [s for s in order if not s.startswith("cat")]
    piv = piv[piv.loc[mo].mean().sort_values().index]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("div", [S2, "#f6f5f2", S1])
    M_ = piv.to_numpy()
    # imshow はセルを 1 枚のラスタ画像として PDF に埋める。pcolormesh の平坦シェーディングなら
    # 同じ見た目のまま四角形のベクタになるので、図全体をラスタなしで保てる。
    nrow_, ncol_ = M_.shape
    im = ax.pcolormesh(np.arange(ncol_ + 1) - 0.5, np.arange(nrow_ + 1) - 0.5, M_,
                       cmap=cmap, vmin=-1, vmax=1, shading="flat", rasterized=False)
    ax.set_xlim(-0.5, ncol_ - 0.5)
    ax.set_ylim(nrow_ - 0.5, -0.5)          # 行 0 を上に（imshow の origin="upper" と同じ）
    # diverging カラーマップは，グレースケールにすると正負が同じ明度に落ちて区別できない。
    # 負のセル（43 / 1,184）に点を重ね，色に依存せず符号が読めるようにする。
    rec("Fig3", "B", "columns (pathways drawn)", piv.shape[1])
    rec("Fig3", "B", "collections drawn", "+".join(sorted(big.collection.unique())))
    rec("Fig3", "B", "cells drawn", M_.size)
    rec("Fig3", "B", "cells below zero (dotted)", int((M_ < 0).sum()))
    ny_, nx_ = np.where(M_ < 0)
    ax.scatter(nx_, ny_, s=5.5, marker="o", color=INK, zorder=5, linewidths=0)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([LAB[s].replace(" (reference)", " (ref)") for s in piv.index])
    ax.set_xticks([])
    ax.set_xlabel(f"{piv.shape[1]} pathways of the three original collections "
                  r"with $\geq$ 50 Group A genes" "\n"
                  r"(ordered by mean mouse $\cos\theta$)")
    ax.tick_params(axis="y", length=0)
    if "KEGG|Focal adhesion" in piv.columns:
        fa = list(piv.columns).index("KEGG|Focal adhesion")
        # 端の列では中央揃えだと文字がパネルの外（カラーバー側）へはみ出す
        ha_ = "right" if fa > ncol_ * 0.75 else ("left" if fa < ncol_ * 0.25 else "center")
        ax.annotate("Focal adhesion", xy=(fa, -0.55), xytext=(fa, -1.5), ha=ha_, va="bottom", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.9))
    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label(r"$\cos\theta$ to the reference axis")
    cb.ax.tick_params( length=3)
    cb.outline.set_visible(False)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=INK,
                              markersize=3.5)],
              labels=[r"$\cos\theta$ < 0"], frameon=False,
              loc="lower left", bbox_to_anchor=(0.0, 1.005), handlelength=1.4,
              handletextpad=0.4)
    ax.set_title("(B)  Alignment varies by pathway, and the ranking of mouse states with it", loc="left", pad=22)
    fig.tight_layout(rect=(0, 0, 1, 1), h_pad=1.0)
    fig.subplots_adjust(hspace=0.52)
    print("Fig3:", savefig(fig, OUT / "Fig3_pathway_alignment")[0].name)
    plt.close(fig)


# ============================================================ Fig 5
def fig5():
    R = pd.read_csv(PWD / "pathway_cos.tsv", sep="\t")
    SEP = pd.read_csv(PWD / "pathway_separation_paired.tsv", sep="\t")
    SUM = pd.read_csv(PWD / "separation_summary.tsv", sep="\t").set_index("restriction")
    col = "auc[exclude shared controls]"
    a = SEP[col].dropna()
    auc_agg = float(SUM.loc["exclude shared controls", "aggregated_auc"])
    auc_panel = 0.813

    BA = pd.read_csv(PWD / "cos_before_after.tsv", sep="\t")
    fig = plt.figure(figsize=(FIGW, 8.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 0.92], width_ratios=[1.55, 1.0])

    # ---- (A) 順位の入れ替わり
    ax = fig.add_subplot(gs[0, :])
    mouse = [k for k, _l, _g in ORDER if not k.startswith("cat")]
    R = R.assign(key=R.collection + "|" + R.pathway)          # 422 セット。名前だけだと 419 に潰れる
    piv = R[R.state.isin(mouse)].pivot_table(index="key", columns="state", values="cos").dropna()
    ranks = piv.rank(axis=1, ascending=False)          # 1 = 最も基準に近い
    med = ranks.median().sort_values()
    top = piv.idxmax(axis=1).value_counts()
    n_p, k = ranks.shape
    Rj = ranks.sum(axis=0).to_numpy()
    W = 12 * ((Rj - Rj.mean()) ** 2).sum() / (n_p ** 2 * (k ** 3 - k))
    rec("Fig4", "A", "pathways ranked", n_p)
    rec("Fig4", "A", "mouse states ranked", k)
    rec("Fig4", "A", "Kendall W", f"{W:.3f}")
    data = [ranks[s].to_numpy() for s in med.index]
    bp = ax.boxplot(data, vert=False, widths=0.62, patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, linewidth=1.6),
                    whiskerprops=dict(color=INK2, linewidth=0.9),
                    capprops=dict(color=INK2, linewidth=0.9),
                    boxprops=dict(facecolor=NEUTRAL, edgecolor=INK2, linewidth=0.8))
    for i, s in enumerate(med.index, start=1):
        if s in top.index:
            rec("Fig4", "A", f"leads share {s}", f"{top[s]/n_p:.0%}")
            ax.text(12.6, i, f"leads {top[s]/n_p:.0%}", ha="left", va="center")
    ax.set_yticklabels([LAB[s] for s in med.index])
    ax.set_xlabel(r"Rank among the 12 mouse states within a pathway (1 = closest to the axis)")
    ax.set_xlim(0.4, 15.4)
    ax.set_xticks(range(1, 13))
    ax.set_title("(A)  No mouse state is consistently closest: ranks reshuffle\n"
                 f"across pathways   (Kendall's $W$ = {W:.3f})", loc="left", pad=6)
    clean(ax, grid="x")
    ax.tick_params(axis="y", length=0)

    # ---- (B) 経路単位 vs 集約
    ax = fig.add_subplot(gs[1, 0])
    ax.hist(a, bins=np.arange(0.30, 1.03, 0.035), color=S1, edgecolor="white",
            linewidth=0.8, zorder=3)
    ym = ax.get_ylim()[1]
    ax.vlines(auc_agg, 0, ym, color=S2, linewidth=2.2, zorder=4)
    ax.vlines(auc_panel, 0, ym, color=INK, linestyle=(0, (4, 3)), linewidth=1.4, zorder=4)
    ax.set_ylim(0, ym * 1.55)
    handles_ = [Line2D([0], [0], color=S2, linewidth=2.2)]
    ax.legend(handles=[Line2D([0], [0], color=S2, linewidth=2.2),
                       Line2D([0], [0], color=INK, linestyle=(0, (4, 3)), linewidth=1.4)],
              labels=[f"collapsed to one score per pathway, centred: {auc_agg:.3f}",
                      f"all 2,016 genes together: {auc_panel:.3f}"], frameon=False, loc="upper left", handlelength=1.4,
              borderaxespad=0.3)
    ax.text(0.315, ym * 1.08, f"per pathway: median {a.median():.3f}, "
            f"IQR {a.quantile(.25):.3f}–{a.quantile(.75):.3f}",
            ha="left", va="center")
    rec("Fig4", "B", "pathways in the AUC histogram", len(a))
    rec("Fig4", "B", "AUC aggregated centred", f"{auc_agg:.3f}")
    rec("Fig4", "B", "AUC whole panel 2,016 genes", f"{auc_panel:.3f}")
    rec("Fig4", "B", "AUC per pathway median", f"{a.median():.3f}")
    rec("Fig4", "B", "AUC per pathway Q1", f"{a.quantile(.25):.3f}")
    rec("Fig4", "B", "AUC per pathway Q3", f"{a.quantile(.75):.3f}")
    ax.set_xlabel("AUC separating within- from cross-species pairs")
    ax.set_ylabel("pathways")
    ax.set_xlim(0.30, 1.02)
    ax.set_title("(B)  Pathway resolution holds it;\n"
                 "aggregation depends on centring", loc="left", pad=6)
    clean(ax, grid="y")

    # ---- (C) サイズ非依存
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(SEP.n_genes, SEP[col], s=22, marker="o", facecolor=S1, edgecolor="white",
               linewidth=0.8, zorder=4)
    ax.axhline(auc_panel, color=INK, linestyle=(0, (4, 3)), linewidth=1.2, zorder=3)
    rho = stats.spearmanr(SEP.n_genes, SEP[col])[0]
    ax.text(0.97, 0.05, f"Spearman $\\rho$ = {rho:+.3f}", transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=MUTED, linewidth=0.8))
    rec("Fig4", "C", "Spearman rho AUC vs set size", f"{rho:+.3f}")
    ax.set_xlabel("Genes in pathway")
    ax.set_ylabel("AUC")
    ax.set_title("(C)  Size explains part of it", loc="left", pad=6)
    clean(ax)
    # ---- (D) 標準化と経路平均を分けた3条件
    ax = fig.add_subplot(gs[2, :])
    cols = ["cos_raw", "cos_centred", "cos_agg_centred", "cos_agg_raw"]
    xs = [0, 1, 2, 3]
    for w, col, lab, dy in [(1, S1, "within species", 0.05), (0, S3, "cat vs mouse", -0.05)]:
        d = BA[BA.within == w]
        for _, r in d.iterrows():                       # 4列目は対照条件なので線でつながない
            ax.plot(xs[:3], [r[c] for c in cols[:3]], color=col, linewidth=0.8,
                    alpha=0.4, zorder=2)
        for x, c in zip(xs, cols):
            ax.scatter([x] * len(d), d[c], s=11, color=col, zorder=3,
                       label=lab if x == 0 else None)
            med = d[c].median()
            ax.plot([x - 0.08, x + 0.08], [med] * 2, color=col, linewidth=2.4, zorder=5,
                    solid_capstyle="butt")
            rec("Fig4", "D", f"median {col_name(c)} [{lab}]", f"{med:.3f}")
            ax.text(x + 0.11, med + (dy * 1.6 if x == 2 else 0), f"{med:.3f}", ha="left", va="center", zorder=6)
    ax.axvline(2.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=1)
    ax.set_xlim(-0.55, 3.62); ax.set_xticks(xs)
    ax.set_xticklabels(["Δ as analysed\n(2,016 genes)",
                        "Δ centred per state\n(= Pearson $r$ of Δ)",
                        "pathway means of\ncentred Δ",
                        "pathway means of\nuncentred Δ"])
    ax.set_ylabel(r"$\cos\theta$ between states")
    ax.legend( frameon=False, loc="lower left", handlelength=1.4,
              borderaxespad=0.3)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.14)
    for x, k in [(0, "cos_raw"), (1, "cos_centred"), (2, "cos_agg_centred"), (3, "cos_agg_raw")]:
        a_ = BA[BA.within == 1][k].to_numpy(); b_ = BA[BA.within == 0][k].to_numpy()
        auc_ = float(np.mean((a_[:, None] > b_[None, :]) + 0.5 * (a_[:, None] == b_[None, :])))
        rec("Fig4", "D", f"AUC {col_name(k)}", f"{auc_:.3f}")
        ax.text(x, hi + 0.10, f"AUC {auc_:.3f}",
                ha="center", va="center")
    ax.set_title("(D)  Removing each state's mean change before averaging\n"
                 "reverses what pathway means show", loc="left", pad=6)
    clean(ax, grid="y")

    fig.tight_layout(rect=(0, 0.01, 1, 1), h_pad=0.8, w_pad=0.8)
    fig.subplots_adjust(hspace=0.46, wspace=0.34)
    print("Fig4:", savefig(fig, OUT / "Fig4_pathway_ranking_and_aggregation")[0].name)
    plt.close(fig)


# ============================================================ Fig 6
def joint_rank(path, genes):
    e = pd.read_csv(HERE / path, sep="\t", index_col=0).reindex(genes)
    return e.rank(pct=True).mean(axis=1, skipna=True)


def fig6():
    ST = pd.read_csv(RES / "reliability" / "orthology_strata.tsv", sep="\t")
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    xy = D[["cat_CKD34", "mouse_2W"]].dropna()
    genes = xy.index
    gA = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    inA = genes.isin(gA)
    cov = {"Mean expression": joint_rank("control_log2cpm.tsv", genes),
           "SE of Δ\n(empirical Bayes)": joint_rank("delta_se_moderated.tsv", genes),
           "Noise of Δ\n(split-half)": joint_rank("delta_se_splithalf.tsv", genes)}

    est = []
    u = pd.read_csv(RES / "groupAB_control_log2cpm" / "check1_unmatched_delta.tsv",
                    sep="\t", index_col=0)["value"]
    est.append({"label": "unmatched", "d": float(u["point_estimate"]),
                "lo": float(u["lo95"]), "hi": float(u["hi95"]), "mk": "D"})
    for tag, lab, mk in [("control_log2cpm", "matched on mean expression", "o"),
                         ("delta_se_moderated", "matched on SE of Δ", "s"),
                         ("delta_se_splithalf", "matched on split-half noise of Δ", "^")]:
        v = pd.read_csv(RES / f"groupAB_{tag}" / "check1_nn_matched.tsv",
                        sep="\t", index_col=0)["value"]
        est.append({"label": lab, "d": -float(v["delta_GroupB_minus_GroupA"]),
                    "lo": -float(v["delta_hi95"]), "hi": -float(v["delta_lo95"]), "mk": mk})
    c = pd.read_csv(RES / "check1" / "check1_nn_matched.tsv", sep="\t", index_col=0)["value"]
    ctrl = {"label": "matched on mean expression",
            "d": float(c["delta_added_minus_core"]), "lo": float(c["delta_lo95"]),
            "hi": float(c["delta_hi95"]), "mk": "x"}

    fig = plt.figure(figsize=(FIGW, 8.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.15])

    # ---- (A) 同一性四分位（写像品質）
    axA = fig.add_subplot(gs[0])
    w = 0.34
    for i_, (cls, col, lab) in enumerate([
            ("cross_species", S3, "Cross-species"),
            ("same_species_diff_dataset", S2, "Same species, different dataset")]):
        d = ST[ST["class"] == cls].sort_values("pid_quartile")
        x = np.arange(4) + (i_ - 0.5) * w
        axA.bar(x, d.median_cos, width=w * 0.9, color=col, edgecolor="white",
                linewidth=0.8, zorder=3, label=lab)
        for xi, v in zip(x, d.median_cos):
            axA.text(xi, v + 0.014, f"{v:.3f}",
                     ha="center", va="bottom")
    axA.set_xticks(range(4))
    axA.set_xticklabels([f"Q{q}" for q in range(1, 5)])
    axA.set_xlabel("Orthologue percent-identity quartile (Q1 lowest, Q4 highest)")
    axA.set_ylabel(r"Median $\cos\theta$")
    axA.set_ylim(0, 0.95)
    axA.annotate("lowest where identity is highest", xy=(2.82, 0.345), xytext=(1.55, 0.14), ha="left", va="bottom", zorder=6,
                 arrowprops=dict(arrowstyle="->", color=S3, linewidth=1.0))
    axA.set_title("(A)  Confounder 1 — orthologue mapping quality\n"
                  "does not order the deficit", loc="left", pad=24)
    axA.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=2,
               handlelength=1.4)
    clean(axA, grid="y")

    # ---- (B) 共変量のバランス
    axB = fig.add_subplot(gs[1])
    ticks, labs = [], []
    for j_, (name, v) in enumerate(cov.items()):
        for i_, (sel, col, hatch) in enumerate([(inA, S1, ""), (~inA, NEUTRAL, "///")]):
            axB.boxplot([v[sel].dropna()], positions=[j_ * 2.4 + i_ * 0.85], widths=0.68,
                        vert=True, patch_artist=True, showfliers=False,
                        medianprops=dict(color=INK, linewidth=1.5),
                        whiskerprops=dict(color=INK2, linewidth=0.8),
                        capprops=dict(color=INK2, linewidth=0.8),
                        boxprops=dict(facecolor=col, edgecolor=INK2, linewidth=0.8,
                                      hatch=hatch))
        ticks.append(j_ * 2.4 + 0.42); labs.append(name)
    axB.set_xticks(ticks); axB.set_xticklabels(labs)
    axB.set_ylabel("Percentile rank among analysed genes")
    axB.set_ylim(0, 1.0)
    axB.axhline(0.5, color=MUTED, linestyle=(0, (4, 3)), linewidth=1.0, zorder=2)
    axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=S1, edgecolor=INK2),
                        plt.Rectangle((0, 0), 1, 1, facecolor=NEUTRAL, edgecolor=INK2,
                                      hatch="///")],
               labels=["Group A", "Group B"], frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, handlelength=1.4)
    axB.set_title("(B)  Confounders 2 and 3 — only mean expression\n"
                  "separates the gene sets", loc="left", pad=24)
    clean(axB, grid="y")

    # ---- (C) forest（Group A/B の4推定量と、手続き対照を分ける）
    axC = fig.add_subplot(gs[2])
    ys = [4.4, 3.4, 2.4, 1.4]
    for y, e in zip(ys, est):
        axC.plot([e["lo"], e["hi"]], [y, y], color=INK2, linewidth=1.4,
                 solid_capstyle="butt", zorder=3)
        axC.scatter([e["d"]], [y], s=64, marker=e["mk"], facecolor=S1,
                    edgecolor="white", linewidth=1.0, zorder=4)
        axC.text(0.398, y, f"{e['d']:+.4f}", ha="right", va="center")
    axC.axhline(0.72, color=MUTED, linewidth=1.0, zorder=2)
    axC.plot([ctrl["lo"], ctrl["hi"]], [0.1, 0.1], color=INK2, linewidth=1.4,
             solid_capstyle="butt", zorder=3)
    axC.scatter([ctrl["d"]], [0.1], s=72, marker="x", color=S2, linewidth=2.0, zorder=4)
    axC.text(0.398, 0.1, f"{ctrl['d']:+.4f}", ha="right", va="center")
    axC.axvline(0, color=INK, linewidth=1.1, zorder=2)
    axC.set_yticks(ys + [0.1])
    axC.set_yticklabels([e["label"] for e in est] + [ctrl["label"]])
    axC.text(-0.085, 5.05, "Group A vs Group B", style="italic",
             ha="left", va="center")
    axC.text(-0.085, 0.60, "Procedure control — symbol- minus orthologue-matched genes",
             style="italic", ha="left", va="center")
    axC.set_xlabel(r"$\Delta\rho$ between the two gene sets (95% gene-bootstrap interval)")
    axC.set_xlim(-0.09, 0.41)
    axC.set_ylim(-0.45, 5.35)
    axC.set_title("(C)  Matching on expression removes part of the difference;\n"
                  "matching on precision removes none", loc="left", pad=6)
    clean(axC, grid="x")
    axC.tick_params(axis="y", length=0)

    fig.tight_layout(rect=(0, 0.01, 1, 1), h_pad=1.2)
    fig.subplots_adjust(hspace=0.75)
    # A（種間の cos）と B・C（Group A/B の比較）は対象が違うので区切りを入れる。
    # 軸の枠ではなく軸ラベル・見出しを含めた外形で測る（見出しに重ねないため）。
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    yA = axA.get_tightbbox(rnd).transformed(fig.transFigure.inverted()).y0
    yB = axB.get_tightbbox(rnd).transformed(fig.transFigure.inverted()).y1
    ymid = (yA + yB) / 2
    fig.add_artist(Line2D([0.04, 0.98], [ymid, ymid], color=MUTED, linewidth=0.9,
                          linestyle=(0, (5, 4)), transform=fig.transFigure))
    print("Fig5:", savefig(fig, OUT / "Fig5_confounders")[0].name)
    plt.close(fig)


def main():
    for f in (fig1, fig2, fig3, fig4, fig5, fig6):
        f()


if __name__ == "__main__":
    main()
