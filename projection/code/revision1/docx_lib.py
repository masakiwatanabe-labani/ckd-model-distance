# -*- coding: utf-8 -*-
"""docx の XML を段落単位で読み書きするための最小限のヘルパー。

数式（m:oMath）は run の一種として段落に入っているので、段落を作り直すときは
その要素をそのまま持ち越す。python-docx は OMML を扱えないので XML を直接触る。
"""
from __future__ import annotations

import html
import re
import shutil
import zipfile
from pathlib import Path

T = re.compile(r"<w:t(?: [^>]*)?>(.*?)</w:t>", re.S)
BLOCK = re.compile(r"<w:tbl>.*?</w:tbl>|<w:p\b[^>]*/>|<w:p\b.*?</w:p>", re.S)
OMATH = re.compile(r"<m:oMathPara[ >].*?</m:oMathPara>|<m:oMath[ >].*?</m:oMath>", re.S)


def text_of(node: str) -> str:
    """段落の可読テキスト。数式は [MATH] に置き換えて位置が分かるようにする。"""
    s = OMATH.sub("∅", node)
    return html.unescape("".join(T.findall(s))).replace(" ", " ").strip()


def raw_text(node: str) -> str:
    return html.unescape("".join(T.findall(node))).replace(" ", " ").strip()


def blocks_of(xml: str) -> list[str]:
    body = re.search(r"<w:body>(.*)</w:body>", xml, re.S).group(1)
    return BLOCK.findall(body)


def n_math(xml: str) -> tuple[int, int]:
    return (len(re.findall(r"<m:oMath[ >]", xml)),
            len(re.findall(r"<m:oMathPara[ >]", xml)))


def esc(s: str) -> str:
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def para(text: str, style: str | None = None, *, bold=False, italic=False,
         align: str | None = None) -> str:
    """素の段落を1つ作る。"""
    pr = []
    if style:
        pr.append(f'<w:pStyle w:val="{style}"/>')
    if align:
        pr.append(f'<w:jc w:val="{align}"/>')
    ppr = f"<w:pPr>{''.join(pr)}</w:pPr>" if pr else ""
    rpr = ""
    if bold or italic:
        rpr = "<w:rPr>" + ("<w:b/>" if bold else "") + ("<w:i/>" if italic else "") + "</w:rPr>"
    runs = ""
    for i, line in enumerate(text.split("\n")):
        if i:
            runs += "<w:r><w:br/></w:r>"
        runs += f'<w:r>{rpr}<w:t xml:space="preserve">{esc(line)}</w:t></w:r>'
    return f"<w:p>{ppr}{runs}</w:p>"


def set_text(node: str, new: str) -> str:
    """数式を含まない段落のテキストだけを差し替える。書式は最初の run のものを使う。"""
    assert "<m:oMath" not in node, "数式を含む段落にこれを使わない"
    m = re.search(r"<w:r\b[^>]*>(?:\s*<w:rPr>.*?</w:rPr>)?", node, re.S)
    rpr = ""
    if m:
        mm = re.search(r"<w:rPr>.*?</w:rPr>", m.group(0), re.S)
        rpr = mm.group(0) if mm else ""
    ppr = re.search(r"<w:pPr>.*?</w:pPr>", node, re.S)
    ppr = ppr.group(0) if ppr else ""
    runs = ""
    for i, line in enumerate(new.split("\n")):
        if i:
            runs += "<w:r><w:br/></w:r>"
        runs += f'<w:r>{rpr}<w:t xml:space="preserve">{esc(line)}</w:t></w:r>'
    return f"<w:p>{ppr}{runs}</w:p>"


def strip_comments(xml: str) -> str:
    """コメント参照を本文から取り除く。"""
    xml = re.sub(r"<w:commentRangeStart[^>]*/>", "", xml)
    xml = re.sub(r"<w:commentRangeEnd[^>]*/>", "", xml)
    xml = re.sub(r"<w:r\b[^>]*>(?:(?!</w:r>).)*?<w:commentReference[^>]*/>.*?</w:r>", "", xml,
                 flags=re.S)
    return xml


def strip_external_refs(zin, drop: set, extra: dict) -> None:
    """雛形から引き継がれる外部参照と第三者のトラッキング ID を落とす（AL2）。

    投稿する docx には、テンプレートを置いていた別の利用者のローカルパス
    （word/_rels/settings.xml.rels の attachedTemplate）と、
    Grammarly の文書 ID（docProps/custom.xml）が残る。どちらも著者と無関係の
    第三者情報なので、本文に触れずに除く。
    """
    names = set(zin.namelist())
    if "word/_rels/settings.xml.rels" in names:
        rels = zin.read("word/_rels/settings.xml.rels").decode("utf-8")
        if "attachedTemplate" in rels:
            rels = re.sub(r"<Relationship\b[^>]*attachedTemplate[^>]*/>", "", rels)
            extra["word/_rels/settings.xml.rels"] = rels
    if "word/settings.xml" in names:
        s = zin.read("word/settings.xml").decode("utf-8")
        if "attachedTemplate" in s:
            extra["word/settings.xml"] = re.sub(r"<w:attachedTemplate\b[^>]*/>", "", s)
    if "docProps/custom.xml" in names:
        drop.add("docProps/custom.xml")
        if "[Content_Types].xml" in names:
            ct = zin.read("[Content_Types].xml").decode("utf-8")
            extra["[Content_Types].xml"] = re.sub(
                r'<Override[^>]*PartName="/docProps/custom\.xml"[^>]*/>', "", ct)
        if "_rels/.rels" in names:
            rr = zin.read("_rels/.rels").decode("utf-8")
            extra["_rels/.rels"] = re.sub(
                r"<Relationship\b[^>]*docProps/custom\.xml[^>]*/>", "", rr)


def write_docx(src: Path, dest: Path, document_xml: str, *, drop=(), extra=None,
               add=None) -> Path:
    """src を雛形に、document.xml を差し替えた docx を書く。

    add は雛形に無いパートを足すときに使う（name -> bytes）。図を埋め込む場合の
    word/media/*.png がこれにあたる。
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    zin = zipfile.ZipFile(src)
    drop, extra, add = set(drop), dict(extra or {}), dict(add or {})
    strip_external_refs(zin, drop, extra)
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename in drop:
                continue
            data = (document_xml.encode("utf-8") if item.filename == "word/document.xml"
                    else zin.read(item.filename))
            if item.filename in extra:
                data = extra[item.filename].encode("utf-8")
            zout.writestr(item, data)
        for name, data in add.items():
            assert name not in zin.namelist(), f"{name} は雛形にすでにある"
            zout.writestr(name, data)
    return dest


EMU_PER_IN = 914400


def image_para(rid: str, px: tuple[int, int], dpi: float, max_in: float,
               *, doc_id: int, name: str, descr: str) -> str:
    """画像 1 点を本文幅に収めて中央に置く段落を作る。

    幅は実寸（px / dpi）と max_in の小さいほうに合わせ、縦横比は保つ。
    """
    w_in = min(px[0] / dpi, max_in)
    h_in = w_in * px[1] / px[0]
    cx, cy = int(round(w_in * EMU_PER_IN)), int(round(h_in * EMU_PER_IN))
    A = 'http://schemas.openxmlformats.org/drawingml/2006/main'
    P = 'http://schemas.openxmlformats.org/drawingml/2006/picture'
    return (
        '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>'
        '<wp:inline distT="0" distB="0" distL="0" distR="0">'
        f'<wp:extent cx="{cx}" cy="{cy}"/>'
        '<wp:effectExtent l="0" t="0" r="0" b="0"/>'
        f'<wp:docPr id="{doc_id}" name="{esc(name)}" descr="{esc(descr)}"/>'
        f'<wp:cNvGraphicFramePr><a:graphicFrameLocks xmlns:a="{A}" noChangeAspect="1"/>'
        '</wp:cNvGraphicFramePr>'
        f'<a:graphic xmlns:a="{A}"><a:graphicData uri="{P}">'
        f'<pic:pic xmlns:pic="{P}">'
        f'<pic:nvPicPr><pic:cNvPr id="0" name="{esc(name)}"/><pic:cNvPicPr/></pic:nvPicPr>'
        f'<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch>'
        '</pic:blipFill>'
        '<pic:spPr><a:xfrm><a:off x="0" y="0"/>'
        f'<a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
        '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
        '</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')
