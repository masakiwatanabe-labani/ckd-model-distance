# -*- coding: utf-8 -*-
"""C案の docx を Markdown 一式に変換する。

C案はこれまで docx だけで存在していて Markdown が無かった。以降の改訂を
Markdown 側で管理し、verify_numbers.py に載せられるようにするための変換。

段落と表を本文の順に取り、見出し（1. / 2.1. など）で節ファイルに分ける。
表は <w:tbl> から行と列を取って pipe 表にする（段落の並びから復元しない）。
"""
from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]

T = re.compile(r"<w:t(?: [^>]*)?>(.*?)</w:t>", re.S)
BLOCK = re.compile(r"<w:tbl>.*?</w:tbl>|<w:p\b[^>]*/>|<w:p\b.*?</w:p>", re.S)
ROW = re.compile(r"<w:tr\b.*?</w:tr>", re.S)
CELL = re.compile(r"<w:tc>.*?</w:tc>", re.S)

# 見出し → 出力ファイル
FILES = [("FRONTMATTER", None), ("INTRODUCTION", "1."), ("RESULTS", "2."),
         ("DISCUSSION", "3."), ("METHODS", "4."), ("CONCLUSIONS", "5.")]


def text_of(node: str) -> str:
    return html.unescape("".join(T.findall(node))).replace(" ", " ").strip()


def table_md(tbl: str) -> str:
    rows = []
    for tr in ROW.findall(tbl):
        cells = [text_of(tc).replace("|", r"\|") for tc in CELL.findall(tr)]
        if any(c for c in cells):
            rows.append(cells)
    if not rows:
        return ""
    n = max(len(r) for r in rows)
    rows = [r + [""] * (n - len(r)) for r in rows]
    out = ["| " + " | ".join(rows[0]) + " |", "|" + "---|" * n]
    for r in rows[1:]:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "Downloads" / "C案 MS.docx"
    dest = HERE / "manuscript_C"
    dest.mkdir(exist_ok=True)

    xml = zipfile.ZipFile(src).read("word/document.xml").decode("utf-8")
    body = re.search(r"<w:body>(.*)</w:body>", xml, re.S).group(1)

    blocks = []
    for b in BLOCK.findall(body):
        if b.startswith("<w:tbl>"):
            md = table_md(b)
            if md:
                blocks.append(("table", md))
        else:
            t = text_of(b)
            if t:
                blocks.append(("text", t))

    # 見出しの検出と節分け
    head = re.compile(r"^(\d+(?:\.\d+)*)\.\s+(.+)$")
    parts: dict[str, list[str]] = {f: [] for f, _ in FILES}
    parts["BACKMATTER"] = []
    parts["SUPPLEMENTARY"] = []
    parts["REFERENCES"] = []
    parts["FIGURES"] = []
    cur = "FRONTMATTER"
    in_refs = seen_back = False
    for kind, t in blocks:
        if kind == "text":
            if t.strip() == "References":
                cur, in_refs = "REFERENCES", True
                parts[cur].append("# References")
                continue
            if in_refs and re.match(r"^\d{1,2}\.\s+[A-Z\[]", t):
                parts["REFERENCES"].append(t)
                continue
            if re.match(r"^Figure \d+\.", t):
                cur, in_refs = "FIGURES", False
                parts[cur].append(t)
                continue
            # 補足へ回すのは、後付け資料の始まり以降だけ。Methods 本文にも
            # 「Figure S2B displays …」のような言及があるので、そこで切り替えない。
            if seen_back and re.match(r"^(Supplementary Results|Table S\d|Figure S\d)", t):
                cur, in_refs = "SUPPLEMENTARY", False
            m = head.match(t)
            if m and not in_refs and cur not in ("SUPPLEMENTARY", "FIGURES") \
                    and len(m.group(1).split(".")) <= 2:
                top = m.group(1).split(".")[0]
                for f, pre in FILES:
                    if pre and pre.startswith(top + "."):
                        cur = f
                        break
                level = "#" if "." not in m.group(1) else "##"
                parts[cur].append(f"{level} {m.group(1)}. {m.group(2)}")
                continue
            if re.match(r"^(Author Contributions|Funding|Institutional|Informed|Data Availability|"
                        r"Acknowledgments|Conflicts of Interest|Supplementary Materials):", t):
                cur, in_refs, seen_back = "BACKMATTER", False, True
            parts[cur].append(t)
        else:
            parts[cur].append(t)

    for f, blocks_ in parts.items():
        if blocks_:
            (dest / f"{f}.md").write_text("\n\n".join(blocks_) + "\n")
    print(f"書き出し: {dest}")
    for p in sorted(dest.glob("*.md")):
        print(f"  {p.name:20s} {len(p.read_text().split()):6d} 語")
    return 0


if __name__ == "__main__":
    sys.exit(main())
