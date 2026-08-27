"""Fig 2: 種間一致とその物差し。

(A) ネコ皮質晩期 vs Pod-TRECK 2W の log2FC 散布図（transcripts）
(B) 同、蛋白（vs Day21）
(C) 種間一致 vs 種内ベンチマークの横棒（モジュールレベルも併記）
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

log = get_logger("23_fig2")
apply_style()
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
FIG = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))

CROSS = "#C44E52"; WITHIN = "#4C72B0"; MODULE = "#8172B3"


def lfc(ds):
    s = DE.loc[ds, "lfc"]
    return s[np.isfinite(s)]


def to_mouse(s):
    o = s.copy(); o.index = pd.Series(s.index).map(CAT2MOUSE).values
    o = o[pd.notna(o.index)]
    o.index = [str(i).upper() for i in o.index]
    return o[~o.index.duplicated()]


def norm(s):
    o = s.copy(); o.index = [str(i).upper() for i in o.index]
    return o[~o.index.duplicated()]


def scatter(ax, cat_ds, mouse_ds, xlab, ylab, title, letter):
    x, y = norm(to_mouse(lfc(cat_ds))), norm(lfc(mouse_ds))
    i = x.index.intersection(y.index)
    xv, yv = x.loc[i].values, y.loc[i].values
    rho = stats.spearmanr(xv, yv)[0]
    ax.axhline(0, color="0.75", lw=.8, zorder=1)
    ax.axvline(0, color="0.75", lw=.8, zorder=1)
    ax.scatter(xv, yv, s=5, alpha=.22, color=CROSS, edgecolor="none", zorder=2, rasterized=True)
    z = np.polyfit(xv, yv, 1)
    xs = np.linspace(np.percentile(xv, .3), np.percentile(xv, 99.7), 50)
    ax.plot(xs, np.polyval(z, xs), color="0.15", lw=1.3, zorder=4)
    ax.set_xlabel(xlab, fontsize=9); ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(f"{letter}  {title}\nSpearman rho = {rho:.3f}   n = {len(i):,}",
                 loc="left", fontsize=10)
    ax.grid(alpha=.18, zorder=0)
    return rho, len(i)


def main():
    fig = plt.figure(figsize=(15.0, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12], wspace=.20)
    axA, axB, axC = (fig.add_subplot(gs[i]) for i in range(3))

    rA, nA = scatter(axA, "cat_rna_ctx_late", "m_rna_2w",
                     "cat cortex, IRIS 3/4 vs control (log2 FC)",
                     "Pod-TRECK 14 d vs control (log2 FC)",
                     "Transcriptome", "A")
    rB, nB = scatter(axB, "cat_prot_ctx_late", "m_prot_d21",
                     "cat cortex protein, IRIS 3/4 vs control (log2 FC)",
                     "Pod-TRECK Day 21 vs control (log2 FC)",
                     "Proteome", "B")

    # ---------------- (C) ----------------
    bench = pd.read_csv(RES / "within_species_benchmarks.csv")

    def b(lbl):
        return float(bench.loc[bench.label == lbl, "rho"].iloc[0])

    mods = pd.read_csv(RES / "module_crossspecies.csv")

    def m(pair):
        return float(mods.loc[mods.pair == pair, "rho"].iloc[0])

    # モジュールレベルはこのパネルに置かない。モデル間ペアでも同等の値が出るため
    # （3.8 参照）、ここに並べると「集約すると種間一致が種内を超える」と誤読される。
    rows = [
        ("cross-species, transcripts", rA, CROSS, "cross-species"),
        ("cross-species, protein", rB, CROSS, "cross-species"),
        ("cat RNA vs protein", b("ネコ皮質 RNA vs protein"), WITHIN, "within-species"),
        ("mouse RNA vs protein", b("マウス RNA2W vs protein D14"), WITHIN, "within-species"),
        ("cat cortex vs medulla (RNA)", b("ネコ 皮質 vs 髄質 RNA"), WITHIN, "within-species"),
    ]
    # 長いカテゴリ名を軸外の ytick に置くと隣のパネル（B）の作図領域に
    # 食い込むため、ラベルはバーの直上に軸内配置する。
    ypos = np.arange(len(rows))
    axC.barh(ypos, [r[1] for r in rows], color=[r[2] for r in rows], height=.52, zorder=3)
    axC.set_yticks([])
    for i, r in enumerate(rows):
        axC.text(0.008, i + 0.30, r[0], va="bottom", ha="left", fontsize=8.2, zorder=5)
        axC.text(r[1] + .014, i, f"{r[1]:.3f}", va="center", fontsize=8.2, zorder=4)

    # 種内ベンチマークの帯（種間がその何割かを示す）
    lo, hi = b("ネコ皮質 RNA vs protein"), b("ネコ 皮質 vs 髄質 RNA")
    # 帯は薄く敷く。ラベルは枠内の最上段（バーの無い帯状の空き）に置き、
    # 注記は枠内の最下段に確保した余白に置く。どちらも軸外に出さないので
    # タイトルとも右枠とも交差しない。
    axC.axvspan(lo, hi, color=WITHIN, alpha=.09, zorder=0)
    axC.set_ylim(-1.15, len(rows) + 0.18)
    axC.text(0.63, len(rows) - 0.06, "within-species benchmark band  0.660 - 0.760",
             ha="center", va="center", fontsize=6.9, color=WITHIN, zorder=5)

    pctA, pctB = 100 * rA / lo, 100 * rB / lo
    axC.text(0.008, -0.72,
             f"cross-species reaches {pctA:.0f}% (transcripts) and {pctB:.0f}% (protein)\n"
             f"of the lowest within-species benchmark (0.660)",
             fontsize=7.4, color="0.25", va="center", ha="left", linespacing=1.9, zorder=5)
    axC.set_xlim(0, 1.0); axC.set_xlabel("Spearman rho", fontsize=9)
    axC.set_title("C  Cross-species agreement against within-species benchmarks",
                  loc="left", fontsize=10, pad=10)
    axC.grid(axis="x", alpha=.22, zorder=0)

    _savefig(fig, FIG / "Fig2_concordance")
    log.info("saved %s", FIG / "Fig2_concordance.png")

    vals = {"A transcripts rho": round(rA, 3), "A n": nA,
            "B protein rho": round(rB, 3), "B n": nB,
            **{f"C {r[0]}": round(r[1], 3) for r in rows},
            "C cross/within(min) transcripts %": round(pctA, 1),
            "C cross/within(min) protein %": round(pctB, 1)}
    pd.Series(vals).to_csv(RES / "fig2_values.csv")
    log.info("%s", vals)
    append_summary("23_fig2 / Fig2", {"図": "results/figures/Fig2_concordance.png (300 dpi)",
                                      "数値": vals}, cfg)


if __name__ == "__main__":
    main()
