# -*- coding: utf-8 -*-
"""BK: 回答書の段落が Markdown と同じ順で docx に入っているか。

build_letters_docx.py が地の文を溜めてから吐いていたため、「**Major 1.**」のような
見出しが、同じ段落にある引用より後ろに出ていた。見出しと引用の並びを突き合わせる。
"""
from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
PAIRS = [("RESPONSE.md", "Response_to_reviewers.docx"),
         ("COVER_LETTER.md", "Cover_letter.docx")]


def key(s: str) -> str:
    return " ".join(s.split())[:60]


def expected(md: str) -> list[tuple[str, str]]:
    """Markdown から、見出し（**…**／#）と引用（> …）を出てくる順に。"""
    out = []
    for line in md.splitlines():
        s = line.strip()
        if s.startswith("> "):
            out.append(("quote", key(re.sub(r"^> \*?|\*$", "", s))))
        elif s.startswith("#"):
            out.append(("head", key(s.lstrip("# "))))
        elif re.fullmatch(r"\*\*(.+?)\*\*", s):
            out.append(("head", key(s.strip("*"))))
    return out


def got(path: Path) -> list[str]:
    x = zipfile.ZipFile(path).read("word/document.xml").decode()
    return [key(re.sub(r"<[^>]+>", "", b))
            for b in re.findall(r"<w:p\b.*?</w:p>", x, re.S)
            if re.sub(r"<[^>]+>", "", b).strip()]


def main() -> int:
    bad, n = [], 0
    for md_name, dx_name in PAIRS:
        md, dx = M / md_name, DEST / dx_name
        if not (md.exists() and dx.exists()):
            continue
        want = expected(md.read_text())
        paras = got(dx)
        pos, last = [], -1
        for kind, k in want:
            n += 1
            hit = next((i for i, p in enumerate(paras) if i > last and p.startswith(k[:40])), None)
            if hit is None:
                bad.append(f"{dx_name}: {kind} が見つからないか順が違う: {k}")
                continue
            pos.append((kind, k, hit))
            last = hit
        # 見出しの直後が、その見出しに属する引用であること
        for i, (kind, k, at) in enumerate(pos[:-1]):
            if kind == "head" and pos[i + 1][0] == "quote" and pos[i + 1][2] < at:
                bad.append(f"{dx_name}: 引用が見出しより前に出ている: {k}")
        print(f"  {dx_name}: 見出しと引用 {len(pos)} 件の並びを確認")
    print(f"\n突き合わせた段落 {n} 件")
    if bad:
        for b in bad[:8]:
            print("  ✗ " + b)
        return 1
    print("  Markdown と同じ順に入っている")
    return 0


if __name__ == "__main__":
    sys.exit(main())
