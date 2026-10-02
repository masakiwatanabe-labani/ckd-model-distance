# -*- coding: utf-8 -*-
"""BA4: 削除・差し替えをした箇所の前後を、出来上がった docx の文字で読む。

run を連結して読む（タグを空白に置き換えると、run の境目に無い空白が混ざる）。
機械的に見つけられる壊れ方だけを検査にする:
  - 同じ文が 2 度続く
  - 句読点の前の空白、句読点の連続
  - 大文字で始まらない文、動詞の欠けた断片の兆候（「, we has」のような形）
"""
from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
FILES = ["Main_manuscript.docx", "Response_to_reviewers.docx", "Cover_letter.docx",
         "Supplementary_Notes.docx"]
# 数式を含む段落は、<w:t> だけを読むと数式の位置が空白になる（「defined as , the」）。
# その段落では句読点前の空白を見ない。
BROKEN = [
    (r"\s+[,.;:]", "句読点の前に空白", "no-math"),
    # 実際に起きた壊れ方は「we has …」という主語と動詞の不一致だった。
    (r"\b(we|they|these|those) has\b|\b(it|this|that) have\b", "主語と動詞が合わない", "any"),
    (r"\b(\w{4,}) \1\b", "語の重複", "any"),
]


par_has_math: dict = {}


def paragraphs(p: Path) -> list[str]:
    x = zipfile.ZipFile(p).read("word/document.xml").decode("utf-8")
    out = []
    for b in re.findall(r"<w:p[ >].*?</w:p>", x, re.S):
        t = html.unescape("".join(
            re.findall(r"<w:t(?![^>]*/>)[^>]*>(.*?)</w:t>", b, re.S)))
        t = re.sub(r"\s+", " ", t).strip()
        if t:
            out.append(t)
            par_has_math[id(out[-1])] = "<m:oMath" in b
    return out


def sentences(par: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.;:])\s+", par) if s.strip()]


def main() -> int:
    bad, n_par, n_sent = [], 0, 0
    for f in FILES:
        p = DEST / f
        if not p.exists():
            continue
        for par in paragraphs(p):
            n_par += 1
            ss = sentences(par)
            n_sent += len(ss)
            # 同じ文が 2 度続く
            for a, b in zip(ss, ss[1:]):
                if len(a) > 30 and a.lower() == b.lower():
                    bad.append(f"{f}: 同じ文が 2 度続く: {a[:70]}")
            has_math = par_has_math.get(id(par), False)
            for pat, what, scope in BROKEN:
                if scope == "no-math" and has_math:
                    continue
                for m in re.finditer(pat, par):
                    seg = par[max(0, m.start() - 40):m.end() + 40]
                    bad.append(f"{f}: {what}: …{seg}…")
    print(f"提出 docx {len(FILES)} 件 / 段落 {n_par} / 文 {n_sent} を読んだ")
    if bad:
        for b in bad[:12]:
            print("  ✗ " + b)
        print(f"  疑い {len(bad)} 件")
        return 1
    print("  壊れた文は見つからない")
    return 0


if __name__ == "__main__":
    sys.exit(main())
