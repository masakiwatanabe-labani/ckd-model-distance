# -*- coding: utf-8 -*-
"""全表の列見出しと中身の対応を検査する（AB）。

査読者 3 が見つけた Table S5 の列ずれは、セルの中に `|` を含む名前
（`elapsed time |Δlog10 h|`、`time | onset compartment`）が Markdown の
pipe 表を通ったときに列区切りとして解釈されたために起きた。
同じ壊れ方が他の表にないかを機械的に調べる。

検査項目:
  1. 見出しの列数と各行の列数が一致するか
  2. セルの中に未エスケープの `|` が残っていないか（変換で列がずれる原因）
  3. `\\` の残骸が入っていないか
  4. 数値であるべき列に文字列が混ざっていないか（列ずれの兆候）
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
MD = HERE / "manuscript_C"
FILES = ["RESULTS.md", "METHODS.md", "SUPPLEMENTARY.md"]
SEP = re.compile(r"\|[\s\-|:]+\|")
NUM = re.compile(r"^[+\-−]?[\d,]*\.?\d+(\s*\(p\s*=\s*[\d.]+\))?$")


def tables(text: str):
    """(キャプション, 行のリスト) を返す。"""
    out, cur, cap, last = [], [], None, ""
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("|"):
            cur.append(s)
        else:
            if cur:
                out.append((cap or last, cur))
                cur = []
            if re.match(r"^(Table|Figure) S?\d+\.", s):
                cap = s[:70]
            elif s:
                last = s[:70]
    if cur:
        out.append((cap or last, cur))
    return [(c, r) for c, r in out if any(SEP.fullmatch(x) for x in r)]


def cells(row: str):
    return [c.strip() for c in re.split(r"(?<!\\)\|", row.strip().strip("|"))]


def main() -> int:
    bad, n_tab, n_cell = [], 0, 0
    for f in FILES:
        for cap, rows in tables((MD / f).read_text()):
            n_tab += 1
            head = cells(rows[0])
            body = [r for r in rows[1:] if not SEP.fullmatch(r)]
            for i, r in enumerate(body):
                c = cells(r)
                n_cell += len(c)
                if len(c) != len(head):
                    bad.append(f"{f} / {cap[:44]}: 行 {i+1} が {len(c)} 列（見出しは {len(head)}）")
                for x in c:
                    if "\\|" in x or re.search(r"(?<!\\)\|", x):
                        bad.append(f"{f} / {cap[:44]}: セルに | が残る: {x[:40]!r}")
                    if x.endswith("\\") or " \\ " in x:
                        bad.append(f"{f} / {cap[:44]}: セルに \\ の残骸: {x[:40]!r}")
            # 数値列らしい列に文字が混ざっていないか
            for j in range(1, len(head)):
                col = [cells(r)[j] for r in body if len(cells(r)) == len(head)]
                if not col:
                    continue
                numlike = sum(bool(NUM.match(v)) for v in col if v not in ("", "—", "-"))
                nonempty = sum(1 for v in col if v not in ("", "—", "-"))
                odd = [v for v in col if v not in ("", "—", "-") and not NUM.match(v)]
                # yes/no や「3 each」のように、設計上そう書く列は誤検出しない
                VOCAB = {"yes", "no", "n/a", "na", "none"}
                designed = all(v.lower() in VOCAB or re.fullmatch(r"\d+ each", v) for v in odd)
                if nonempty >= 3 and odd and not designed and numlike >= nonempty * 0.5:
                    bad.append(f"{f} / {cap[:44]}: 列 {j+1}「{head[j][:24]}」に非数値が混在 "
                               f"{odd[:2]}")
    print(f"表 {n_tab} 件 / セル {n_cell} 件を検査")
    if bad:
        print(f"問題 {len(bad)} 件:")
        for b in bad:
            print(f"  ✗ {b}")
    else:
        print("  問題なし")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
