# -*- coding: utf-8 -*-
"""AW2: Table S10 の段階表を結果ファイルから作り直す（手で直さない）。"""
from __future__ import annotations
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
SUP = HERE / "manuscript_C" / "SUPPLEMENTARY.md"
GS = pd.read_csv(HERE / "results" / "roundR" / "groupA_detection_summary.tsv", sep="\t")


def main() -> int:
    want = "\n".join(["| Step | Genes |", "|---|---|"] +
                     [f"| {r.step} | {int(r.genes):,} |" for r in GS.itertuples()])
    t = SUP.read_text()
    # 表の行だけを、空行に当たるまで取る。以降を飲み込まないよう 1 行ずつ数える。
    lines = t.splitlines(keepends=True)
    start = next((i for i, l in enumerate(lines) if l.startswith("| Step | Genes |")), None)
    assert start is not None, "段階表が見つからない"
    end = start
    while end < len(lines) and lines[end].startswith("|"):
        end += 1
    cur = "".join(lines[start:end]).rstrip("\n")

    class _M:
        pass
    m = _M()
    m.start = lambda: sum(len(x) for x in lines[:start])
    m.end = lambda: sum(len(x) for x in lines[:end])
    if cur == want:
        print("段階表は既に最新")
        return 0
    SUP.write_text(t[:m.start()] + want + "\n" + t[m.end():])
    print("段階表を作り直した:")
    for line in want.splitlines()[2:]:
        print("  " + line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
