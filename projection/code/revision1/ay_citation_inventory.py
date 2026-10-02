# -*- coding: utf-8 -*-
"""AY1: 投稿版以降に追加・書き換えた文に現れる [n] を、第一著者と並べて一覧にする（目視用）。

an_block_map.tsv の rewritten / new の段落だけを対象にする（unchanged は投稿版のまま）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
OUT = HERE / "results" / "roundR" / "citation_inventory.tsv"


def refs() -> dict:
    out = {}
    for line in (M / "REFERENCES.md").read_text().splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m:
            first = m.group(2).split(";")[0].strip()
            title = re.split(r"\.\s", m.group(2), maxsplit=1)
            out[m.group(1)] = (first, (title[1] if len(title) > 1 else "")[:66])
    return out


def main() -> int:
    R = refs()
    BM = pd.read_csv(HERE / "results" / "roundR" / "an_block_map.tsv", sep="\t",
                     dtype=str, keep_default_na=False)
    changed = BM[BM.action.str.startswith(("rewritten", "new"))]
    rows, seen = [], set()
    for r in changed.itertuples():
        for m in re.finditer(r"\[(\d+(?:\s*[,–-]\s*\d+)*)\]", str(r.text)):
            for tok in re.split(r"\s*[,–-]\s*", m.group(1)):
                if tok in R and (tok, r.section) not in seen:
                    seen.add((tok, r.section))
                    rows.append({"ref": int(tok), "section": r.section,
                                 "first_author": R[tok][0], "title": R[tok][1],
                                 "context": str(r.text)[:56]})
    # letter に新しく書いた参照も加える
    lt = (M / "RESPONSE.md").read_text()
    for m in re.finditer(r"reference (\d+)", lt):
        tok = m.group(1)
        if tok in R and (tok, "RESPONSE") not in seen:
            seen.add((tok, "RESPONSE"))
            rows.append({"ref": int(tok), "section": "RESPONSE",
                         "first_author": R[tok][0], "title": R[tok][1],
                         "context": " ".join(lt[max(0, m.start() - 50):m.end()].split())[-52:]})
    T = pd.DataFrame(rows).sort_values(["ref", "section"], kind="stable")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 210)
    pd.set_option("display.max_colwidth", 52)
    print(T.to_string(index=False))
    print(f"\n{len(T)} 件 → {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
