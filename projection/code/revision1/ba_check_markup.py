# -*- coding: utf-8 -*-
"""BA3: 出来上がった docx / xlsx に Markdown の記号が残っていないか検査する。

docx は文字だけを見る（XML のタグではなく、読者が見るテキスト）。
"""
from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
# 記号 → 説明。太字は「**」が本文に出てはいけない
MARKS = [("`", "バッククォート"),
         ("**", "太字の記号"),
         ("## ", "見出しの記号"),
         ("# ", "見出しの記号"),
         ("](", "リンクの記号"),
         ("|---", "表の区切り")]


def docx_text(p: Path) -> str:
    x = zipfile.ZipFile(p).read("word/document.xml").decode("utf-8")
    t = re.sub(r"<[^>]+>", " ", x)
    return html.unescape(re.sub(r"\s+", " ", t))


def xlsx_text(p: Path) -> str:
    import openpyxl
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    out = []
    for s in wb.sheetnames:
        ws = wb[s]
        for row in ws.iter_rows(max_row=4, values_only=True):
            out += [str(c) for c in row if c is not None]
    return " ".join(out)


def main() -> int:
    bad, n = [], 0
    # 添付の定量表は Markdown から作っていない（受領したままの受託解析の出力）。
    # 列見出しに「# of Grouped proteins」があるので走査の対象外にする。
    SKIP = {"Supplementary_Data_S1_PodTRECK_proteome.xlsx"}
    for p in sorted(DEST.glob("*")):
        if p.name.startswith("~$") or p.name in SKIP:
            continue
        if p.suffix == ".docx":
            t = docx_text(p)
        elif p.suffix == ".xlsx":
            t = xlsx_text(p)
        else:
            continue
        n += 1
        for sym, what in MARKS:
            if sym in t:
                i = t.index(sym)
                bad.append(f"{p.name}: {what} 「{sym}」 …{t[max(0,i-40):i+40]}…")
    print(f"提出ファイル {n} 件の文字を、Markdown の記号 {len(MARKS)} 種で走査")
    if bad:
        for b in bad[:10]:
            print("  ✗ " + b)
        return 1
    print("  記号は残っていない")
    return 0


if __name__ == "__main__":
    sys.exit(main())
