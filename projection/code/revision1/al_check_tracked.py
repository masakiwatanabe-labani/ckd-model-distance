# -*- coding: utf-8 -*-
"""AL1: 変更履歴版 Main_manuscript_tracked.docx を機械で検査する。

Word の「文書の比較」は元文書と変更後文書を取り違えやすいので、
承諾版・却下版のテキストを自分で作って突き合わせる。
変更履歴版そのものは編集しない。読むだけ。

  承諾（accept all）: <w:del>…</w:del> を落とし、<w:ins> の中身を残す
  却下（reject all）: <w:ins>…</w:ins> を落とし、<w:del> の中身（w:delText）を残す
"""
from __future__ import annotations

import difflib
import re
import sys
import zipfile
from pathlib import Path

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
TRACKED = DEST / "Main_manuscript_tracked.docx"
CLEAN = DEST / "Main_manuscript.docx"
ORIGINAL = Path(__file__).resolve().parents[2] / "data" / "base" / "base_manuscript.docx"


def doc_xml(p: Path) -> str:
    return zipfile.ZipFile(p).read("word/document.xml").decode("utf-8")


def plain(xml: str) -> str:
    """段落ごとのテキスト。比較のため空白を潰す。"""
    out = []
    for m in re.finditer(r"<w:p[ >].*?</w:p>", xml, re.S):
        # <w:t/> のような自己終了タグを開始タグと誤認すると、次の </w:t> までの
        # XML をまるごと本文として拾ってしまう。除外する。
        t = "".join(re.findall(
            r"<w:(?:t|delText)(?![^>]*/>)[^>]*>(.*?)</w:(?:t|delText)>", m.group(0), re.S))
        t = re.sub(r"\s+", " ", t).strip()
        if t:
            out.append(t)
    return "\n".join(out)


# 範囲マーカーは文字を持たないので落とす（移動の記録に使われる）
RANGE = re.compile(r"</?w:(?:moveFrom|moveTo)Range(?:Start|End)[^>]*/?>")


def _unwrap(x: str, name: str) -> str:
    """<w:name …> と </w:name> だけを外す。

    要素名の直後が空白・`>`・`/` であることを要求する。これを入れないと
    `<w:delText …>` が `<w:del…>` として巻き込まれ、削除された文字が消える。
    """
    return re.sub(rf"</?w:{name}(?=[ >/])[^>]*>", "", x)


def _drop(x: str, name: str) -> str:
    """<w:name …>…</w:name> を中身ごと落とす。"""
    return re.sub(rf"<w:{name}(?=[ >])[^>]*>.*?</w:{name}>", "", x, flags=re.S)


def accept(xml: str) -> str:
    """すべて承諾。削除と移動元を落とし、挿入と移動先の中身を残す。"""
    x = _drop(xml, "del")
    x = _drop(x, "moveFrom")
    x = _unwrap(x, "ins")
    x = _unwrap(x, "moveTo")
    return RANGE.sub("", x)


def reject(xml: str) -> str:
    """すべて却下。挿入と移動先を落とし、削除と移動元の中身を残す。"""
    x = _drop(xml, "ins")
    x = _drop(x, "moveTo")
    # 先に delText を通常の w:t に直してから、del の殻を外す
    x = x.replace("<w:delText", "<w:t").replace("</w:delText>", "</w:t>")
    x = _unwrap(x, "del")
    x = _unwrap(x, "moveFrom")
    return RANGE.sub("", x)


def diff(a: str, b: str, n=6) -> list[str]:
    d = [l for l in difflib.unified_diff(a.splitlines(), b.splitlines(), lineterm="", n=0)
         if l[:1] in "+-" and l[:3] not in ("+++", "---")]
    return d[:n]


def main() -> int:
    if not TRACKED.exists():
        print(f"変更履歴版がまだ無い: {TRACKED}\n"
              "  Word の「文書の比較」で作ってから、もう一度これを実行する", file=sys.stderr)
        return 2
    bad = 0
    xml = doc_xml(TRACKED)
    acc, rej = plain(accept(xml)), plain(reject(xml))
    cln, org = plain(doc_xml(CLEAN)), plain(doc_xml(ORIGINAL))

    ok1 = acc == cln
    print(f"1. 承諾版 == Main_manuscript.docx : {'一致' if ok1 else '不一致'}")
    if not ok1:
        bad += 1
        print(f"   段落 {len(acc.splitlines())} 対 {len(cln.splitlines())}")
        for l in diff(acc, cln):
            print("   ", l[:150])
    ok2 = rej == org
    print(f"2. 却下版 == 投稿稿 : {'一致' if ok2 else '不一致'}")
    if not ok2:
        bad += 1
        print(f"   段落 {len(rej.splitlines())} 対 {len(org.splitlines())}")
        for l in diff(rej, org):
            print("   ", l[:150])
    if not ok1 and not ok2:
        # 向きの取り違えを直接調べる
        if plain(reject(xml)) == cln and plain(accept(xml)) == org:
            print("   → 元文書と変更後文書が逆になっている。比較をやり直すこと")

    a_x = accept(xml)
    n_math = len(re.findall(r"<m:oMath[ >]", a_x))
    n_block = len(re.findall(r"<m:oMathPara[ >]", a_x))
    ok3 = (n_math, n_block) == (33, 13)
    print(f"3. 承諾版の数式 : oMath {n_math} / ブロック {n_block} "
          f"{'（33 / 13 で一致）' if ok3 else '✗ 33 / 13 のはず'}")
    bad += 0 if ok3 else 1

    ok4 = "<w:lnNumType" in xml
    print(f"4. 通し行番号 lnNumType : {'あり' if ok4 else '✗ 無い'}")
    bad += 0 if ok4 else 1

    # 6. 黒以外の文字色は、削除・移動元の中（元稿の書式を保つ箇所）だけであること
    stray = []
    for m in re.finditer(r'<w:r\b[^>]*>(?:(?!</w:r>).)*?<w:color w:val="(?!000000)'
                         r'[0-9A-Fa-f]{6}"[^>]*/>.*?</w:r>', xml, re.S):
        before = xml[:m.start()]
        depth = {}
        for o in re.findall(r"<(w:ins|w:del|w:moveFrom|w:moveTo)\b", before):
            depth[o] = depth.get(o, 0) + 1
        for c in re.findall(r"</(w:ins|w:del|w:moveFrom|w:moveTo)>", before):
            depth[c] = depth.get(c, 0) - 1
        if not {k for k, v in depth.items() if v > 0} & {"w:del", "w:moveFrom"}:
            txt = "".join(re.findall(r"<w:(?:t|delText)[^>]*>(.*?)</w:(?:t|delText)>",
                                     m.group(0), re.S))
            stray.append(txt[:40])
    print(f"6. 黒以外の文字色 : {'削除・移動元の中だけ' if not stray else '✗ ' + str(stray[:3])}")
    bad += 0 if not stray else 1

    authors = sorted(set(re.findall(r'w:author="([^"]*)"', xml)))
    print(f"7. 変更履歴の作成者 : {authors if authors else '（変更が無い）'}")
    if len(authors) != 1:
        print("   ✗ 作成者が 1 名でない。別の名前か空欄が混じっている")
        bad += 1
    n_ins = len(re.findall(r"<w:ins[ >]", xml))
    n_del = len(re.findall(r"<w:del[ >]", xml))
    n_mf = len(re.findall(r"<w:moveFrom[ >]", xml))
    n_mt = len(re.findall(r"<w:moveTo[ >]", xml))
    print(f"   挿入 {n_ins} / 削除 {n_del} / 移動元 {n_mf} / 移動先 {n_mt}")
    moved = re.findall(r"<w:moveTo[ >].*?</w:moveTo>", xml, re.S)
    for m in moved:
        tx = re.sub(r"\s+", " ", "".join(
            re.findall(r"<w:t[^>]*>(.*?)</w:t>", m, re.S))).strip()
        if tx:
            print(f"   移動された内容: {tx[:88]}")
    print(f"\n要対応 {bad} 件")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
