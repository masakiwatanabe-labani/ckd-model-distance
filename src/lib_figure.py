"""全図に共通のスタイルと保存。

論文図 7 点（Figure 1-6 と Figure S1）はすべてこのモジュールの apply_style() だけで
書体・級数・線幅・文字色を決める。図ごとに fontsize や色を指定しないこと。
ばらつきはここを直せば 7 点同時に直る、という状態を保つための約束である。

級数（役割ごと、いずれも rcParams で与える）:

  8 pt   目盛りラベル、凡例、図中の注記      font.size / xtick.labelsize / legend.fontsize
  9 pt   軸ラベル                             axes.labelsize
  10 pt  パネル見出し (A) (B) ...、bold       axes.titlesize / axes.titleweight

8 pt を下限とし、6 pt 台・7 pt 台は使わない。数式の添字は本文級数の 0.7 倍で描かれる
（9 pt なら 6.3 pt）ので、ラベルに $x_{i}$ の形を使わないこと。太字はパネル見出しと、
キャプションが名指ししている数値ラベル（EMPH を付けたもの）に限る。図中の文字は
すべて黒で、薄いグレーの文字は使わない。

図全体の見出し（suptitle）は置かない。MDPI はキャプションを図の下に組むので、
同じ文が図の上下に二重に出てしまう。見出しの文言はキャプション側が持つ。
figure.titlesize / figure.titleweight は将来 suptitle を使う場合に備えて残してある。

書体は Arial 一本。macOS の Arial.ttf / Arial Bold.ttf を pdf.fonttype 42 で
サブセット埋め込みする。本文は Palatino（serif）だが、図が sans で本文が serif の
組み合わせは通常であり、合わせる必要はない。

  pdf.fonttype = 42
      matplotlib の既定は Type 3。Type 3 はグリフを独自演算子で埋め込むため
      Illustrator ではテキストとして扱えないことが多い。42 は TrueType を
      埋め込み、文字情報が保持されるので選択・編集できる。
  svg.fonttype = "none"
      SVG ではアウトライン化せず <text> 要素のまま出す。
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---- 書体 ----------------------------------------------------------------
# Arial のみ。日本語グリフが要る探索用の図のために Hiragino Sans を後ろに置くが、
# 論文図 7 点に日本語は無いので埋め込まれるのは Arial だけになる。
FONT_FAMILY = "Arial"
FONT_STACK = [FONT_FAMILY, "Hiragino Sans"]

# ---- 級数（本文から参照できるように名前を付けておく）----------------------
FS_ANNOT = 8      # 図中の注記、目盛りラベル、凡例
FS_LABEL = 9      # 軸ラベル
FS_PANEL = 10     # パネル見出し
FS_SUPTITLE = 11  # 図の見出し

# ---- 色 ------------------------------------------------------------------
TEXT = "#000000"        # 図中の文字はすべてこれ。グレーの文字は使わない
RULE = "#000000"        # 軸・目盛り
GRID = "#dddcd8"        # 罫のみ薄く

# キャプションが名指ししている数値ラベルに付ける強調。ほかに bold は使わない。
EMPH = {"fontweight": "bold"}

MIN_LINEWIDTH = 0.8     # 軸・目盛り・データの線とマーカー縁の下限


def apply_style() -> None:
    matplotlib.rcParams.update({
        "pdf.fonttype": 42,      # TrueType: テキストとして編集できる
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "pdf.compression": 6,

        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "axes.unicode_minus": False,
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_FAMILY,
        "mathtext.it": f"{FONT_FAMILY}:italic",
        "mathtext.bf": f"{FONT_FAMILY}:bold",
        "mathtext.default": "it",

        "font.size": FS_ANNOT,
        "axes.labelsize": FS_LABEL,
        "axes.titlesize": FS_PANEL,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "xtick.labelsize": FS_ANNOT,
        "ytick.labelsize": FS_ANNOT,
        "legend.fontsize": FS_ANNOT,
        "legend.title_fontsize": FS_ANNOT,
        "figure.titlesize": FS_SUPTITLE,
        "figure.titleweight": "bold",
        "figure.labelsize": FS_LABEL,

        "text.color": TEXT,
        "axes.labelcolor": TEXT,
        "axes.titlecolor": TEXT,
        "axes.edgecolor": RULE,
        "xtick.color": RULE,
        "ytick.color": RULE,
        "xtick.labelcolor": TEXT,
        "ytick.labelcolor": TEXT,
        "legend.labelcolor": TEXT,
        "grid.color": GRID,

        "axes.linewidth": MIN_LINEWIDTH,
        "xtick.major.width": MIN_LINEWIDTH,
        "ytick.major.width": MIN_LINEWIDTH,
        "xtick.minor.width": MIN_LINEWIDTH,
        "ytick.minor.width": MIN_LINEWIDTH,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "grid.linewidth": 0.7,
        "lines.linewidth": 1.8,
        "lines.markersize": 5.0,
        "lines.markeredgewidth": 0.9,
        "patch.linewidth": MIN_LINEWIDTH,
        "boxplot.boxprops.linewidth": MIN_LINEWIDTH,
        "boxplot.whiskerprops.linewidth": MIN_LINEWIDTH,
        "boxplot.capprops.linewidth": MIN_LINEWIDTH,
        "boxplot.medianprops.linewidth": 1.6,

        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "legend.borderaxespad": 0.3,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })


# MDPI の本文幅。図はこれを超えない。
MAX_WIDTH_IN = 6.3


def savefig(fig, path, dpi: int = 300, formats=("png", "pdf"), **kw) -> list[Path]:
    """同じ図を複数形式で保存する。path の拡張子は無視して stem を使う。

    保存後に実寸を測り、MAX_WIDTH_IN を超えていれば例外にする。文字を大きくした
    ぶんを図幅で吸収してしまわないための歯止め。
    """
    path = Path(path)
    bbox = kw.pop("bbox_inches", "tight")
    pad = kw.pop("pad_inches", 0.02)
    facecolor = kw.pop("facecolor", "white")
    out = []
    for ext in formats:
        p = path.with_suffix(f".{ext}")
        fig.savefig(p, dpi=dpi, bbox_inches=bbox, pad_inches=pad,
                    facecolor=facecolor, **kw)
        out.append(p)
    w, h = measure(fig, bbox, pad)
    if w > MAX_WIDTH_IN + 1e-3:
        raise ValueError(f"{path.name}: 幅 {w:.2f} in が上限 {MAX_WIDTH_IN} in を超える。"
                         "余白とパネル間隔で詰めること（図の物理寸法を広げない）")
    return out


def measure(fig, bbox="tight", pad: float = 0.02) -> tuple[float, float]:
    """保存される PDF の実寸（インチ）を返す。"""
    if bbox != "tight":
        return tuple(fig.get_size_inches())
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    return bb.width + 2 * pad, bb.height + 2 * pad
