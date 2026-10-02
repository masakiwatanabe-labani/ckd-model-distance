# -*- coding: utf-8 -*-
"""AW: 切り落とした SUPPLEMENTARY.md の末尾（Table S11・S12 と図の legend）を復元する。

aw_cascade_table.py の正規表現が段階表以降を飲み込んでしまったため、
その直前に組み上げた Supplementary_Notes.docx（キャプションと legend）と
Supplementary_Tables.xlsx（表本体）から組み直す。どちらも同じ Markdown から作られている。
"""
from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parents[2]
SUP = HERE / "manuscript_C" / "SUPPLEMENTARY.md"
DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"


def notes_paras() -> list[str]:
    x = zipfile.ZipFile(DEST / "Supplementary_Notes.docx").read("word/document.xml").decode()
    out = []
    for b in re.findall(r"<w:p[ >].*?</w:p>", x, re.S):
        t = "".join(re.findall(r"<w:t(?![^>]*/>)[^>]*>(.*?)</w:t>", b, re.S))
        t = html.unescape(re.sub(r"\s+", " ", t)).strip()
        if t:
            out.append(t)
    return out


def sheet_rows(name: str) -> list[str]:
    wb = openpyxl.load_workbook(DEST / "Supplementary_Tables.xlsx", read_only=True,
                               data_only=True)
    ws = wb[name]
    raw = []
    for r in ws.iter_rows(min_row=3, values_only=True):
        cells = ["" if c is None else str(c) for c in r]
        while cells and cells[-1] == "":
            cells.pop()
        if cells:
            raw.append(cells)
    if not raw:
        return []
    # 見出しの列数にそろえる。元の表には全幅の小見出し行（「Within dataset (n = 45)」）が
    # あり、空セルを落とすと列数が合わなくなる。
    width = len(raw[0])
    rows = ["| " + " | ".join(raw[0]) + " |", "|" + "---|" * width]
    for cells in raw[1:]:
        cells = cells[:width] + [""] * (width - len(cells))
        rows.append("| " + " | ".join(cells) + " |")
    return rows


def main() -> int:
    t = SUP.read_text().rstrip("\n")
    if "Table S11." in t:
        print("既に復元されている")
        return 0
    paras = notes_paras()
    cap = {}
    for p in paras:
        m = re.match(r"^(Table S\d+|Figure S\d+)\.", p)
        if m:
            cap[m.group(1)] = p
    parts = [t, ""]
    for n in ("Table S11", "Table S12"):
        parts += [cap[n], "", "\n".join(sheet_rows(n)), ""]
    parts += ["## Supplementary Figure legends", ""]
    for n in ("Figure S1", "Figure S2", "Figure S3", "Figure S4", "Figure S5"):
        parts += [cap[n], ""]
    SUP.write_text("\n".join(parts).rstrip("\n") + "\n")
    print(f"復元した: Table S11・S12 と図の legend 5 件（{len(SUP.read_text().splitlines())} 行）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
