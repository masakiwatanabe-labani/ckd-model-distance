# -*- coding: utf-8 -*-
"""AR1: 節順を投稿先の書式に戻す（AG R1-5 の撤回）。

  Results 3 → 2 / Discussion 4 → 3 / Materials and Methods 2 → 4 / Conclusions 5 そのまま
apply_AG_sections.py の逆写像。見出しと「Section N.M」形式の参照だけを書き換える。
R1-4 で Methods へ移した 2 段落は Methods に残す。
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
FILES = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md", "RESPONSE.md"]
MAP = {"3": "2", "4": "3", "2": "4"}          # AG の {"2":"3","3":"4","4":"2"} の逆
REF = re.compile(r"(Sections?\s+)([2-4]\.\d{1,2})((?:\s*(?:,|and)\s*[2-4]\.\d{1,2})*)")
HEAD = re.compile(r"(?m)^(#{1,3}\s*)([2-4])(\.\d{1,2})?\.")
NUMREF = re.compile(r"\b([2-4])\.(\d{1,2})\b")


def main() -> int:
    n_head = n_ref = 0
    for f in FILES:
        p = M / f
        if not p.exists():
            continue
        t = p.read_text()

        def head_sub(m):
            nonlocal n_head
            n_head += 1
            return f"{m.group(1)}\x01{MAP[m.group(2)]}\x01{m.group(3) or ''}."

        def ref_sub(m):
            nonlocal n_ref
            body = m.group(2) + (m.group(3) or "")
            n_ref += len(NUMREF.findall(body))
            return m.group(1) + NUMREF.sub(
                lambda x: f"\x01{MAP[x.group(1)]}\x01.{x.group(2)}", body)

        t = HEAD.sub(head_sub, t)
        t = REF.sub(ref_sub, t)
        p.write_text(re.sub(r"\x01(\d)\x01", r"\1", t))
    print(f"  見出し {n_head} 件 / 参照 {n_ref} 件を戻した")
    for f in ("RESULTS.md", "DISCUSSION.md", "METHODS.md", "CONCLUSIONS.md"):
        print(f"  {f}: {(M / f).read_text().splitlines()[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
