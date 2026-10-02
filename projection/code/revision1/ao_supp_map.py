# -*- coding: utf-8 -*-
"""AO2: 投稿版の補足番号と改訂版の補足番号の対応を、キャプションの内容で突き合わせる。

査読者は投稿版の番号で書いているので、回答では併記する必要がある。
番号の付け替えではなく内容で対応づけるので、途中の振り直しの経緯に依存しない。
"""
from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
BASE = HERE / "data" / "base" / "base_supplementary_notes.docx"
MC = HERE / "manuscript_C"
OUT = HERE / "results" / "roundR" / "supplementary_number_map.tsv"


def sim(a: str, b: str) -> float:
    import difflib
    return difflib.SequenceMatcher(None, a.lower(), b.lower(), autojunk=False).ratio()


def base_caps() -> dict:
    x = zipfile.ZipFile(BASE).read("word/document.xml").decode("utf-8")
    out = {}
    for b in re.findall(r"<w:p[ >].*?</w:p>", x, re.S):
        t = html.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", b, re.S)))
        t = re.sub(r"\s+", " ", t).strip()
        m = re.match(r"^(Table|Figure) S(\d+)\.\s*(.*)$", t)
        if m:
            out[(m.group(1), int(m.group(2)))] = m.group(3)
    return out


def new_caps() -> dict:
    out = {}
    for b in (MC / "SUPPLEMENTARY.md").read_text().split("\n\n"):
        s = " ".join(b.split())
        m = re.match(r"^(Table|Figure) S(\d+)\.\s*(.*)$", s)
        if m:
            out[(m.group(1), int(m.group(2)))] = m.group(3)
    return out


def main() -> int:
    old, new = base_caps(), new_caps()
    rows, used = [], set()
    for (kind, n), txt in sorted(old.items()):
        best, hit = 0.0, None
        for (k2, n2), t2 in new.items():
            if k2 != kind:
                continue
            s = sim(txt[:300], t2[:300])
            if s > best:
                best, hit = s, n2
        rows.append({"kind": kind, "submitted": n, "revised": hit if best >= 0.5 else "",
                     "similarity": round(best, 3),
                     "caption": txt[:88]})
        if best >= 0.5:
            used.add((kind, hit))
    for (kind, n), txt in sorted(new.items()):
        if (kind, n) not in used:
            rows.append({"kind": kind, "submitted": "", "revised": n, "similarity": "",
                         "caption": txt[:88]})
    T = pd.DataFrame(rows).sort_values(["kind", "revised", "submitted"], kind="stable")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_colwidth", 90)
    print(T.to_string(index=False))
    miss = T[(T.submitted != "") & (T.revised == "")]
    print(f"\n投稿版にあって改訂版に対応が無いもの: {len(miss)}")
    print(f"改訂版で新設されたもの: {int((T.submitted == '').sum())}")
    print(f"書き出し: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
