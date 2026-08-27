"""図の共通スタイルと保存。

Illustrator でテキストを編集できる PDF を出すための設定をまとめる。

  pdf.fonttype = 42
      matplotlib の既定は Type 3。Type 3 はグリフを独自演算子で埋め込むため
      Illustrator ではテキストとして扱えないことが多い。42 は TrueType を
      埋め込み、文字情報が保持されるので選択・編集できる。
  svg.fonttype = "none"
      SVG ではアウトライン化せず <text> 要素のまま出す。

フォントは macOS に標準である Helvetica を第一候補にする。埋め込まれた
フォントが手元に無いと Illustrator が代替フォントに置き換えるため、
一般的な書体にしておくほうが扱いやすい。日本語グリフは Hiragino Sans に
フォールバックする（探索用の図に日本語ラベルが残っているため）。
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FONT_STACK = ["Helvetica", "Arial", "Hiragino Sans", "DejaVu Sans"]


def apply_style() -> None:
    matplotlib.rcParams.update({
        "pdf.fonttype": 42,      # TrueType: テキストとして編集できる
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "pdf.compression": 6,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "axes.unicode_minus": False,
    })


def savefig(fig, path, dpi: int = 300, formats=("png", "pdf"), **kw) -> list[Path]:
    """同じ図を複数形式で保存する。path の拡張子は無視して stem を使う。"""
    path = Path(path)
    out = []
    for ext in formats:
        p = path.with_suffix(f".{ext}")
        fig.savefig(p, dpi=dpi, bbox_inches=kw.pop("bbox_inches", "tight"),
                    facecolor=kw.pop("facecolor", "white"), **kw)
        out.append(p)
    return out
