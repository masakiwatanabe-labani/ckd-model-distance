"""論文図 7 点を書き出し、納品フォルダにまとめる。

PDF はベクタ（ラスタ画像なし）。PNG は PDF のラスタ化ではなく 150 dpi で描き直す。
図に印字した数値は figure_printed_numbers.tsv に、書体と寸法は figure_typography.tsv に落ちる。
"""
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "code" / "revision1"))
sys.path.insert(0, str(HERE.parent / "src"))

DEST = Path(os.path.expanduser("~/Desktop/CKD_xspecies_fig3"))
(DEST / "figures").mkdir(parents=True, exist_ok=True)
(DEST / "figures_png").mkdir(parents=True, exist_ok=True)

import pandas as pd  # noqa: E402
import lib_figure  # noqa: E402
import make_figures as mf  # noqa: E402
import cosine_attenuation_sim as sim  # noqa: E402

SIM = pd.read_csv(HERE / "results" / "revision1" / "cosine_attenuation_sim.tsv", sep="\t")
MAIN = (mf.fig1, mf.fig2, mf.fig3, mf.fig4, mf.fig5, mf.fig6)

# ---- ベクタ PDF と 300 dpi PNG（既定の出力先 projection/manuscript/figures/）----
for f in MAIN:
    f()
sim.plot(SIM)
tsv = mf.dump_printed()

# ---- 150 dpi PNG は描き直す ----
_orig = lib_figure.savefig


def png150(fig, path, dpi=300, formats=("png", "pdf"), **kw):
    return _orig(fig, DEST / "figures_png" / Path(path).name, dpi=150, formats=("png",), **kw)


lib_figure.savefig = png150
mf.savefig = png150
sim.savefig = png150
mf.PRINTED.clear()
for f in MAIN:
    f()
sim.plot(SIM)

for p in sorted(mf.OUT.glob("Fig*.pdf")):
    shutil.copy2(p, DEST / "figures" / p.name)
shutil.copy2(tsv, DEST / "figure_printed_numbers.tsv")
print("PDF :", sorted(p.name for p in (DEST / "figures").glob("*.pdf")))
print("PNG :", sorted(p.name for p in (DEST / "figures_png").glob("*.png")))
