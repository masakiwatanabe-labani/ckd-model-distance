"""確認用の 150 dpi PNG を出力する。PDF をラスタ化せず、図を描き直して保存する。"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(HERE.parent / 'src'))

OUT = Path(os.path.expanduser('~/Desktop/CKD_xspecies_figures/figures_png'))
OUT.mkdir(parents=True, exist_ok=True)

import lib_figure  # noqa: E402

_orig = lib_figure.savefig


def png150(fig, path, dpi=300, formats=("png", "pdf"), **kw):
    p = Path(path)
    return _orig(fig, OUT / p.name, dpi=150, formats=("png",), **kw)


lib_figure.savefig = png150
import make_figures  # noqa: E402
make_figures.savefig = png150
make_figures.main()

import cosine_attenuation_sim as sim  # noqa: E402
sim.main()
print('150 dpi PNG:', sorted(p.name for p in OUT.glob('*.png')))
