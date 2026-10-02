"""Resolve a tracked document to its accept-all view and paint the inserted text blue.

    python bh_build_blue.py <tracked document.xml> <output document.xml>

All changes are accepted; only the characters that were inside w:ins become blue
(0000FF). Deletions and move sources are dropped, move destinations keep their own
color. Theme colors are removed from the blue runs: a themeColor attribute wins over
w:val in Word and would paint them black again.

Masaki の tracked_tools/build_blue.py が置かれたらそちらを使う。これはその仕様どおりの
実装で、呼び出し方と引数は同じ。
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path.home() / "Desktop" / "tracked_tools"))
from blk import NS, q, load  # noqa: E402

BLUE = "0000FF"
# CT_RPr の並びで w:color より前に来る要素。ここの直後に color を挿す。
BEFORE_COLOR = ["rStyle", "rFonts", "b", "bCs", "i", "iCs", "caps", "smallCaps",
                "strike", "dstrike", "outline", "shadow", "emboss", "imprint",
                "noProof", "snapToGrid", "vanish", "webHidden"]
THEME_ATTRS = ("themeColor", "themeTint", "themeShade")


def lname(e) -> str:
    return etree.QName(e).localname


def rpr_of(run):
    """run の rPr を返す。無ければ先頭に作る（w:r も m:r も rPr が先頭）。"""
    rpr = run.find("w:rPr", NS)
    if rpr is None:
        rpr = etree.Element(q("w:rPr"))
        run.insert(0, rpr)
    return rpr


def paint(run) -> bool:
    """run の文字色を青にする。themeColor は外す。色を付けたら True。"""
    if run.find("w:t", NS) is None and run.find("m:t", NS) is None \
            and run.find("w:delText", NS) is None:
        return False                      # 文字を持たない run（改ページなど）は触らない
    rpr = rpr_of(run)
    color = rpr.find("w:color", NS)
    if color is None:
        color = etree.Element(q("w:color"))
        i = 0
        for c in rpr:
            if lname(c) in BEFORE_COLOR:
                i = list(rpr).index(c) + 1
        rpr.insert(i, color)
    color.set(q("w:val"), BLUE)
    for a in THEME_ATTRS:
        if color.get(q("w:" + a)) is not None:
            del color.attrib[q("w:" + a)]
    return True


def accept(body):
    """変更をすべて承諾した本文にする（check_tracked.resolve の accept と同じ規則）。"""
    b = copy.deepcopy(body)
    drop, keep = {"del", "moveFrom"}, {"ins", "moveTo"}

    for tr in list(b.iter(q("w:tr"))):
        trpr = tr.find("w:trPr", NS)
        if trpr is not None and "del" in {lname(c) for c in trpr}:
            tr.getparent().remove(tr)
    for p in list(b.iter(q("w:p"))):
        rpr = p.find("w:pPr/w:rPr", NS)
        if rpr is not None and len(rpr) and lname(rpr[0]) in drop:
            p.getparent().remove(p)
    for t in list(b.iter(q("w:tbl"))):
        if t.find("w:tr", NS) is None:
            t.getparent().remove(t)

    # w:ins の中の文字だけを青にする。moveTo は元の色のまま。
    painted = 0
    for ins in b.iter(q("w:ins")):
        if lname(ins.getparent()) in ("rPr", "trPr"):
            continue                      # 段落記号・行の印には色が無い
        for run in ins.iter():
            if lname(run) == "r" and paint(run):
                painted += 1

    for tag in ("ins", "del", "moveFrom", "moveTo"):
        for w in list(b.iter(q("w:" + tag))):
            par = w.getparent()
            if par is None:
                continue
            if lname(par) in ("rPr", "trPr"):
                par.remove(w)
                continue
            if tag in drop:
                par.remove(w)
            else:
                i = par.index(w)
                for c in list(w):
                    par.insert(i, c)
                    i += 1
                par.remove(w)
    for tag in ("moveFromRangeStart", "moveFromRangeEnd",
                "moveToRangeStart", "moveToRangeEnd"):
        for e in list(b.iter(q("w:" + tag))):
            e.getparent().remove(e)
    assert keep and not list(b.iter(q("w:delText"))), "承諾後に delText が残っている"
    return b, painted


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    tree, body, _ = load(src)
    new_body, painted = accept(body)
    body.getparent().replace(body, new_body)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    tree.write(dst, xml_declaration=True, encoding="UTF-8", standalone=True)
    print(f"blue runs: {painted}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
