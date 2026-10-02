# -*- coding: utf-8 -*-
"""文献を「鍵」で引く。番号を書き込んだ箇所を無くすための共通部品（BJ2）。

番号は REFERENCES.md から題名の手がかりで引く。文献が増減して番号が動いても、
補足ノート・補足表・回答書が追随する。本文の Markdown では renumber_refs_C.py が
番号を振り直すので、鍵を使うのは本文の外（補足・回答書）である。
"""
from __future__ import annotations

import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
REFS = HERE / "manuscript_C" / "REFERENCES.md"

# 鍵 → 書誌に必ず含まれる手がかり
NEEDLE = {
    "li": "Integrated multi-omics analysis of renal metabolism in domestic cats",
    "companion": "Integrated proteomic–transcriptomic profiling identifies",
    "iri": "Molecular characterization of the transition from acute to chronic",
    "ensembl": "Harrison",
    "podtreck": "An inducible and podocyte-specific diphtheria toxin receptor",
    "ercb": "Tissue transcriptome-driven identification of epidermal growth factor",
}
TOKEN = re.compile(r"\[REF:([a-z0-9_]+)\]")


def table() -> dict[str, str]:
    """鍵 → 現在の番号。"""
    txt = REFS.read_text()
    out = {}
    for key, needle in NEEDLE.items():
        for line in txt.splitlines():
            m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
            if m and needle in m.group(2):
                out[key] = m.group(1)
                break
        else:
            raise SystemExit(f"REFERENCES.md に見つからない鍵: {key}（{needle}）")
    return out


def entry(key: str) -> str:
    """鍵に対応する書誌の 1 行（番号を除く）。検査で中身を照合するのに使う。"""
    needle = NEEDLE[key]
    for line in REFS.read_text().splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m and needle in m.group(2):
            return m.group(2)
    raise SystemExit(f"REFERENCES.md に見つからない鍵: {key}")


def resolve(text: str) -> str:
    """文中の [REF:key] を現在の番号に置き換える。"""
    t = table()

    def sub(m):
        k = m.group(1)
        if k not in t:
            raise SystemExit(f"未知の文献キー: [REF:{k}]")
        return f"[{t[k]}]"

    return TOKEN.sub(sub, text)


# --------------------------------------------------------------- BK: 投稿時の番号
BASE_DOCX = HERE / "data" / "base" / "base_manuscript.docx"
DOI = re.compile(r"https?://doi\.org/(\S+?)\.?$")


def submitted_table() -> dict[str, str]:
    """投稿稿（土台 docx）の文献一覧。番号 → 書誌。

    投稿先の雛形は番号を自動採番するので、本文中に数字は無い。参考文献スタイルの
    段落を出てくる順に数えて番号とする。
    """
    import zipfile
    x = zipfile.ZipFile(BASE_DOCX).read("word/document.xml").decode()
    out, n = {}, 0
    for b in re.findall(r"<w:p\b.*?</w:p>", x, re.S):
        st = re.search(r'<w:pStyle w:val="([^"]+)"', b)
        t = " ".join(re.sub(r"<[^>]+>", "", b).split())
        if st and "reference" in st.group(1).lower() and t:
            n += 1
            out[str(n)] = t
    return out


def _doi(line: str) -> str | None:
    m = DOI.search(line.strip())
    return m.group(1).lower().rstrip(".") if m else None


def current_number_of(submitted_no: str) -> tuple[str | None, str]:
    """投稿時の番号 → （改訂版の番号, その書誌）。DOI で突き合わせる。"""
    sub = submitted_table()
    if submitted_no not in sub:
        return None, ""
    d = _doi(sub[submitted_no])
    if not d:
        return None, sub[submitted_no]
    for line in REFS.read_text().splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m and _doi(m.group(2)) == d:
            return m.group(1), m.group(2)
    return None, sub[submitted_no]


def quote_renumbering(quote_text: str) -> list[tuple[str, str, str]]:
    """査読票の原文に出る [n] を（投稿時の番号, 改訂版の番号, 書誌）に開く。"""
    out = []
    for m in re.finditer(r"\[(\d+)\]", quote_text):
        old = m.group(1)
        new, entry = current_number_of(old)
        if new and new != old:
            out.append((old, new, entry))
    return out
