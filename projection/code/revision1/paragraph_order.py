"""Part 7. 本文の段落順マニフェストを書き出す。

docx 側で段落が別の節に紛れ込んでいる箇所を安全に直せるよう、Markdown 側の正しい
出現順を 1 行 1 ブロックで出す。文章そのものは変更しない。

列:
  seq      通し番号
  section  所属する節番号（2.2、4.13 など。節の外は 0）
  kind     heading / text / equation / table / figure / caption
  first60  空白と記号を除いた先頭 60 文字
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
MD = HERE / "manuscript"
OUT = HERE / "extracted"
FILES = [("1. Introduction", "INTRODUCTION.md"),
         ("2. Results", "RESULTS.md"),
         ("3. Discussion", "DISCUSSION.md"),
         ("4. Materials and Methods", "METHODS.md"),
         ("5. Conclusions", "CONCLUSIONS.md")]


def norm60(s: str) -> str:
    """空白と記号を除いた先頭 60 文字。docx 側の段落と照合するための鍵。"""
    t = re.sub(r"[^0-9A-Za-zÀ-ɏͰ-Ͽ]", "", s)
    return t[:60]


def kind_of(block: str) -> str:
    b = block.lstrip()
    if b.startswith("#"):
        return "heading"
    if b.startswith("|"):
        return "table"
    if b.startswith("$$") or b.startswith(r"\("):
        return "equation"
    if re.match(r"\*\*(Table|Figure) ", b):
        return "caption"
    if b.startswith("---") or b.startswith("###"):
        return "heading"
    return "text"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [("seq", "section", "kind", "first60")]
    seq = 0
    for top, fn in FILES:
        section = top.split(".")[0]
        for block in (MD / fn).read_text().split("\n\n"):
            b = block.strip()
            if not b:
                continue
            m = re.match(r"^##\s+(\d+\.\d+)\.", b)
            if m:
                section = m.group(1)
            elif re.match(r"^#\s+(\d+)\.", b):
                section = re.match(r"^#\s+(\d+)\.", b).group(1)
            seq += 1
            rows.append((str(seq), section, kind_of(b), norm60(b)))
    p = OUT / "paragraph_order.tsv"
    p.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    print(p, f"({len(rows) - 1} blocks)")
    from collections import Counter
    c = Counter(r[2] for r in rows[1:])
    print("  kinds:", dict(c))
    secs = [r[1] for r in rows[1:]]
    order = []
    for s in secs:
        if not order or order[-1] != s:
            order.append(s)
    dup = [s for s in set(order) if order.count(s) > 1]
    print("  sections in order:", " ".join(order))
    if dup:
        print("  ** 同じ節番号が離れた位置に二度現れる:", sorted(dup))
    return 0


if __name__ == "__main__":
    sys.exit(main())
