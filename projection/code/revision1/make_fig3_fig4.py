"""Round 4: Figure 3 と Figure 4 だけを描き直し、印字した数値を TSV に落とす。

Figure 3 は F-1（パネルB を従来3コレクションの 74 セットに戻す）と
F-2（パネルA の総数をコレクション込みのキーで数え、422 / 244 にする）を反映する。
Figure 4 はパネルA の集計キーだけが変わる。印字値のうち IRI 28 d のリード率が
21% から 22% になり、本文 Section 2.3 の記述と一致する。

PNG は PDF のラスタ化ではなく 150 dpi で描き直す。
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))

DEST = Path(os.path.expanduser("~/Desktop/CKD_xspecies_fig3"))
(DEST / "figures").mkdir(parents=True, exist_ok=True)
(DEST / "figures_png").mkdir(parents=True, exist_ok=True)

import make_figures as mf  # noqa: E402

# --- ベクタ PDF（既定の出力先 projection/manuscript/figures/）---
mf.fig4()          # Figure 3
mf.fig5()          # Figure 4
tsv = mf.dump_printed()

# --- 150 dpi PNG は描き直し ---
import lib_figure  # noqa: E402

_orig = lib_figure.savefig


def png150(fig, path, dpi=300, formats=("png", "pdf"), **kw):
    return _orig(fig, DEST / "figures_png" / Path(path).name, dpi=150, formats=("png",), **kw)


lib_figure.savefig = png150
mf.savefig = png150
mf.PRINTED.clear()
mf.fig4()
mf.fig5()

# --- 納品フォルダへ PDF をコピー ---
import shutil  # noqa: E402

for name in ("Fig3_pathway_alignment.pdf", "Fig4_pathway_ranking_and_aggregation.pdf"):
    shutil.copy2(mf.OUT / name, DEST / "figures" / name)
shutil.copy2(tsv, DEST / "figure_printed_numbers.tsv")
print("done:", sorted(p.name for p in (DEST / "figures").glob("*")),
      sorted(p.name for p in (DEST / "figures_png").glob("*")))
