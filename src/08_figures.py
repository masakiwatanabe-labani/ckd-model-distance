"""主要図の生成。日本語フォントが無い環境を想定し、図中ラベルは英語で作る。"""
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
from lib_stats import load_config, get_logger  # noqa: E402

log = get_logger("08_figures")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
# 論文用の Fig1-6 と混ざらないよう、探索段階の図は supplementary/ に出す。
FIG = RES / "supplementary"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
C2M = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))


def series(ds, to_mouse=False):
    s = DE.loc[ds, "lfc"]
    s = s[np.isfinite(s)]
    if to_mouse:
        s = s.set_axis(pd.Series(s.index).map(C2M))
        s = s[s.index.notna()]
    s.index = [str(i).upper() for i in s.index]
    return s[~s.index.duplicated()]


def fig_crossspecies():
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    for a, (cat_ds, mouse_ds, title) in zip(ax, [
            ("cat_rna_ctx_late", "m_rna_2w", "Transcriptome"),
            ("cat_prot_ctx_late", "m_prot_d21", "Proteome")]):
        x, y = series(cat_ds, True), series(mouse_ds)
        idx = x.index.intersection(y.index)
        rho = stats.spearmanr(x[idx], y[idx])[0]
        a.hexbin(x[idx], y[idx], gridsize=70, bins="log", cmap="viridis", mincnt=1)
        a.axhline(0, color="w", lw=.5); a.axvline(0, color="w", lw=.5)
        a.set_title(f"{title}\nSpearman rho = {rho:.3f} (n={len(idx)})")
        a.set_xlabel("Cat cortex, CKD3/4 vs Ctrl")
        a.set_ylabel("Mouse Pod-TRECK vs Ctrl")
    plt.tight_layout(); plt.savefig(FIG / "fig_crossspecies.png", dpi=180); plt.close()


def fig_benchmarks():
    b = pd.read_csv(RES / "within_species_benchmarks.csv")
    x = pd.read_csv(RES / "crossspecies_concordance.csv")
    x = x[x.mapping == "ortholog"]
    key = x[x.label.isin(["cat_rna_ctx_late x m_rna_2w", "cat_prot_ctx_late x m_prot_d21"])]
    lab = list(b.label) + ["CROSS-SPECIES " + l for l in key.label]
    val = list(b.rho) + list(key.rho)
    col = ["#9aa0a6"] * len(b) + ["#e8710a"] * len(key)
    o = np.argsort(val)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh([lab[i] for i in o], [val[i] for i in o], color=[col[i] for i in o])
    ax.set_xlim(0, 1); ax.set_xlabel("Spearman rho of log2FC vectors")
    ax.set_title("Cross-species similarity vs internal benchmarks")
    plt.tight_layout(); plt.savefig(FIG / "fig_benchmarks.png", dpi=180); plt.close()


def fig_module_stage():
    e = pd.read_csv(RES / "module_eigengene_by_stage.csv")
    for tis in ["ctx", "med"]:
        d = e[e.tissue == tis]
        if d.empty:
            continue
        fig, ax = plt.subplots(figsize=(7, max(4, 0.28 * len(d))))
        y = np.arange(len(d))
        for col, c, m in [("Control", "#9aa0a6", "o"), ("CKD1/2", "#1a73e8", "o"), ("CKD3/4", "#e8710a", "o")]:
            ax.scatter(d[col], y, color=c, label=col, s=34, marker=m, zorder=3)
        for i, (_, r) in enumerate(d.iterrows()):
            ax.plot([r["Control"], r["CKD3/4"]], [i, i], color="#ccc", lw=1, zorder=0)
        ax.set_yticks(y); ax.set_yticklabels(d.module, fontsize=7.5)
        ax.axvline(0, color="k", lw=.5); ax.legend(fontsize=7.5)
        ax.set_xlabel("module eigengene"); ax.set_title(f"Stage trajectory ({tis})")
        plt.tight_layout(); plt.savefig(FIG / f"fig_module_stage_{tis}.png", dpi=180); plt.close()


def fig_progression_axis():
    p = pd.read_csv(RES / "progression_axis_crossspecies.csv")
    d = p[p.cat_axis == "ネコ皮質進行軸"] if "ネコ皮質進行軸" in set(p.cat_axis) else p
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.barh(range(len(d)), d.rho, color=["#34a853" if v > 0 else "#d93025" for v in d.rho])
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d.mouse_axis, fontsize=8)
    ax.axvline(0, color="k", lw=.6)
    ax.set_xlabel("Spearman rho with cat progression axis")
    plt.tight_layout(); plt.savefig(FIG / "fig_progression_axis.png", dpi=180); plt.close()


if __name__ == "__main__":
    for fn in [fig_crossspecies, fig_benchmarks, fig_module_stage, fig_progression_axis]:
        try:
            fn()
            log.info("ok: %s", fn.__name__)
        except Exception as exc:  # noqa: BLE001
            log.warning("skip %s: %s", fn.__name__, exc)
