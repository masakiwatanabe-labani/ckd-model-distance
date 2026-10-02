# -*- coding: utf-8 -*-
"""BH: 青字版（Main_manuscript_blue.docx）の検査。

  1. 文字が Main_manuscript.docx とブロックごとに厳密一致すること
  2. 青の文字列が、変更履歴版の w:ins の文字列と完全一致すること
  3. 変更履歴の印（ins / del / moveFrom / moveTo / delText）が残っていないこと
  4. 本文の文字色が黒と青だけで、青に themeColor が付いていないこと
"""
from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path.home() / "Desktop" / "tracked_tools"))
from blk import NS, q, text_of  # noqa: E402

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
BLUE_VAL = "0000FF"
MARKS = ("ins", "del", "moveFrom", "moveTo", "delText",
         "moveFromRangeStart", "moveFromRangeEnd", "moveToRangeStart", "moveToRangeEnd")

FAIL: list[str] = []
OK = 0


def chk(what, ok, detail=""):
    global OK
    if ok:
        OK += 1
    else:
        FAIL.append(f"{what}{(' — ' + detail) if detail else ''}")
    print(f"  {'ok ' if ok else '✗  '} {what}{(' — ' + detail) if detail else ''}")


def body_of(p: Path):
    x = zipfile.ZipFile(p).read("word/document.xml")
    return etree.fromstring(x).find(q("w:body"))


def lname(e) -> str:
    return etree.QName(e).localname


def blocks(body):
    return [c for c in body if lname(c) != "sectPr"]


def norm(s: str) -> str:
    return " ".join(s.split())


def run_texts(body, pred) -> list[str]:
    """pred を満たす run の文字を、出てくる順に並べる。"""
    out = []
    for r in body.iter():
        if lname(r) != "r":
            continue
        t = text_of(r)
        if t and pred(r):
            out.append(t)
    return out


def is_blue(r) -> bool:
    c = r.find("w:rPr/w:color", NS)
    return c is not None and (c.get(q("w:val")) or "").upper() == BLUE_VAL


def inside_ins(r) -> bool:
    for a in r.iterancestors():
        if lname(a) == "ins":
            return True
        if lname(a) in ("del", "moveFrom"):
            return False
    return False


def main() -> int:
    blue, clean, tracked = (DEST / "Main_manuscript_blue.docx",
                            DEST / "Main_manuscript.docx",
                            DEST / "Main_manuscript_tracked.docx")
    for p in (blue, clean, tracked):
        if not p.exists():
            print(f"無い: {p}", file=sys.stderr)
            return 2
    B, C, T = body_of(blue), body_of(clean), body_of(tracked)

    # 1. 文字がブロックごとに厳密一致
    bb, cb = blocks(B), blocks(C)
    same_n = len(bb) == len(cb)
    chk("1. ブロック数が修正稿と同じ", same_n, f"{len(bb)} / {len(cb)}")
    diff = [i for i, (a, b) in enumerate(zip(bb, cb))
            if lname(a) != lname(b) or norm(text_of(a)) != norm(text_of(b))]
    chk("1. 全ブロックの文字が修正稿と一致", same_n and not diff,
        "" if not diff else f"{len(diff)} ブロックが違う（最初は {diff[0]}）")
    if diff:
        for i in diff[:3]:
            print(f"      blue : {norm(text_of(bb[i]))[:120]}")
            print(f"      clean: {norm(text_of(cb[i]))[:120]}")

    # 2. 青の文字列 == 変更履歴版の w:ins の文字列
    got = run_texts(B, is_blue)
    want = run_texts(T, inside_ins)
    chk("2. 青の run の数が w:ins の run の数と同じ", len(got) == len(want),
        f"{len(got)} / {len(want)}")
    chk("2. 青の文字列が w:ins の文字列と完全一致", got == want,
        "" if got == want else
        f"最初の相違: 青『{(got + [''])[next((i for i in range(max(len(got), len(want))) if (got[i:i+1] or [None]) != (want[i:i+1] or [None])), 0)][:60] if got else ''}』")
    chk("2. つないだ文字列も一致", "".join(got) == "".join(want),
        f"{len(''.join(got)):,} 文字 / {len(''.join(want)):,} 文字")

    # 3. 変更履歴の印が残っていない
    x = zipfile.ZipFile(blue).read("word/document.xml").decode()
    left = {m: len(re.findall(rf"<w:{m}[ />]", x)) for m in MARKS}
    chk("3. 変更履歴の印が残っていない", not any(left.values()),
        ", ".join(f"{k}={v}" for k, v in left.items() if v))
    chk("3. 著者名つきの改訂属性が残っていない",
        not re.search(r'<w:(ins|del|moveFrom|moveTo)\b[^>]*w:author', x))

    # 4. 文字色は黒と青だけ
    vals = set(re.findall(r'<w:color [^>]*w:val="([^"]*)"', x))
    chk("4. 本文の文字色が黒と青だけ", vals <= {"000000", BLUE_VAL}, f"{sorted(vals)}")
    theme = re.findall(rf'<w:color [^>]*w:val="{BLUE_VAL}"[^>]*theme', x)
    chk("4. 青の文字に themeColor が付いていない", not theme, f"{len(theme)} 件")
    n_blue = len(re.findall(rf'<w:color w:val="{BLUE_VAL}"', x))
    print(f"\n  青の run {len(got)} 件 / 青の色指定 {n_blue} 件 / "
          f"青の文字 {len(''.join(got)):,} 文字")

    print(f"\n合格 {OK} / 不合格 {len(FAIL)}")
    for f in FAIL:
        print("  ✗", f)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
