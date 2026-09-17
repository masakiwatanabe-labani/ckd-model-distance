"""Round 4: 図 3 と本文の食い違いのうち、本文側で直すべき 1 か所。

Section 2.3 は 422 セット x 16 状態 = 6,752 セルを Figure 3B の指示で書いていたが、
パネルB が描くのは従来 3 コレクションの 74 セット分 1,184 セルである。
解析のセル数と、パネルが表示するセル数を分けて書く。
数値そのものは変えていない（422 を正とする判断は F-2 のとおり）。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from edit_helper import rp

R = 'manuscript/RESULTS.md'

rp(R, "Each of the 16 states was projected onto the reference axis using only that pathway's "
      "genes, giving 6,752 cells (Figure 3B).",
      "Each of the 16 states was projected onto the reference axis using only that pathway's "
      "genes, giving 6,752 pathway-by-state values. Figure 3B displays the 1,184 of them that "
      "fall in the 74 sets it shows.")
print('Round 4 applied')
