# -*- coding: utf-8 -*-
"""AT1: 「Table S#」の引用が、その番号のキャプションの主題と合っているか検査する。

番号だけの照合では、§2.2 が S5 を S6 と引いていた種類の誤りが通ってしまう。
引用のある文と、引用先のキャプションが、特徴語を共有しているかを見る。
人が確認して残すと決めた箇所は ALLOW に登録する（トーン検査と同じ運用）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
LAB = HERE / "results" / "roundR" / "supp_labels.tsv"
FILES = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md", "RESPONSE.md"]
STOP = set("""a an and are as at be been before between by can do does for from has have in into is
it its not of on or over that the their them there these this to under was were when which with
each every all both one two three four five same such than then they we our us also only its'
table tables figure figures supplementary section sections analysis analyses results result
values value given used using use reported report describe described shown show data set sets""".split())
# 人が確認して残すと決めた引用（ファイル:番号:引用先）
ALLOW = {
    # 一覧・凡例など、主題語を共有しない位置での引用
    ("BACKMATTER.md", "Table S1"), ("BACKMATTER.md", "Figure S1"),
    ("SUPPLEMENTARY.md", "Table S1"),
    # 以下は AT1 で 1 件ずつ内容を確認した。引用先は正しいが、その段落が
    # キャプションと語を共有していないだけ。
    ("METHODS.md", "Table S1"),        # 75 ペアの部分集合を指す（§4.11）
    ("METHODS.md", "Figure S2"),       # 経路レベルの整列の図（§4.13）
    ("RESPONSE.md", "Table S9"),       # R3 major 3 の併記（per-gene standardization）
    ("RESPONSE.md", "Table S1"),       # R3 major 5 の段落（ペア部分集合の言及）
}


def terms(s: str) -> set:
    return {w for w in re.findall(r"[A-Za-z][A-Za-z0-9\-]{3,}", s.lower()) if w not in STOP}


def sentences(text: str):
    """段落単位で見る。引用は短い従属節に置かれることが多く、文単位では狭すぎる。"""
    for para in text.split("\n\n"):
        flat = " ".join(para.split())
        if flat.strip():
            yield flat.strip()


def main() -> int:
    T = pd.read_csv(LAB, sep="\t")
    cap = {r.citation: r.caption for r in T.itertuples()}
    known = set(cap)
    bad, n = [], 0
    for f in FILES:
        p = M / f
        if not p.exists():
            continue
        for s in sentences(p.read_text()):
            # 他誌論文の補足への引用は対象外
            s = re.sub(r"Supplementary (Table|Figure) S\d+ of (that paper|the companion study)",
                       "", s)
            # letter では査読票の引用と「旧番号（新番号 in this revision）」が投稿版の番号。
            # 新番号だけを検査対象にする。
            if f == "RESPONSE.md":
                s = re.sub(r"> \*.*?\*", "", s)   # 段落は 1 行に潰してある
                s = re.sub(r"\b(Table|Figure) S\d+ \((Table|Figure) (S\d+) in this revision\)",
                           r"\1 \3", s)
            for m in re.finditer(r"\b(Table|Figure) S(\d+)[A-D]?\b", s):
                cite = f"{m.group(1)} S{m.group(2)}"
                if cite not in known:
                    bad.append(f"{f}: 実在しない {cite}")
                    continue
                n += 1
                if (f, cite) in ALLOW:
                    continue
                # 引用先のキャプション自身の中の引用は数えない
                if re.match(rf"^{cite}\.", s):
                    continue
                shared = terms(s) & terms(cap[cite])
                if not shared:
                    bad.append(f"{f}: {cite} の主題と合わない引用: {s[:96]}")
    print(f"補足への引用 {n} 件を、キャプションの主題と突き合わせた")
    if bad:
        for b in bad[:12]:
            print("  ✗ " + b)
        print(f"  合わない引用 {len(bad)} 件")
        return 1
    print("  すべて主題を共有している")
    return 0


if __name__ == "__main__":
    sys.exit(main())
