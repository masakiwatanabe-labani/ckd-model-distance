# -*- coding: utf-8 -*-
"""C_MS_clean.docx を土台に、最終 Markdown の本文を当てた Main_manuscript.docx を作る。

数式（m:oMath）は Word のオブジェクトで Markdown 側に持てないので、原典の段落 XML を
そのまま持ち越す。原典のブロック列を順に歩き、

  数式を含む段落  → そのまま（<w:t> の綴りだけ直す）
  表              → Markdown の表から作り直す
  本文の段落      → 対応する Markdown 段落の文へ差し替え、対応が無ければ落とす
  空段落          → そのまま

としたうえで、Markdown にあって原典に無い段落を、順序どおりに差し込む。
"""
from __future__ import annotations

import difflib
import html
import re
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import docx_lib as D  # noqa: E402

SRC = Path.home() / "Downloads" / "C_MS_clean.docx"
MD = HERE / "manuscript_C"
# 投稿先の書式は Introduction → Results → Discussion → Materials and Methods → Conclusions。
# 投稿した元稿がすでにこの順序だったので、AG R1-5 の移動は撤回した（AR1）。
ORDER = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "REFERENCES.md", "FIGURES.md"]
# 数式を含む段落は XML ごと持ち越すので、その中のテキストだけをここで直す。
# 数式段落と Markdown の差はこの4種類だけであることを確認済み。
SPELL = [("normalised", "normalized"), ("standardisation", "standardization"),
         ("standardised", "standardized"), ("summarised", "summarized"),
         ("characterised", "characterized"), ("centred", "centered"),
         ("uncentred", "uncentered"), ("analysed", "analyzed"), ("ischaemia", "ischemia"),
         ("organised", "organized"), ("recognised", "recognized"),
         ("Supplementary Results S", "Supplementary Note S"),
         ("(GEO accession GSE104954)", "(GEO accession GSE104954) [37]")]


def norm(s: str) -> str:
    s = re.sub(r"[^a-z0-9 ]", " ", s.lower())
    return " ".join(s.split())


def md_blocks():
    """Markdown を本文順の段落列にする。$$ 数式は原典側が持つので落とす。"""
    out = []
    for f in ORDER:
        raw_text = (MD / f).read_text()
        # 文献は 1 件 1 行で書いてあるので、空行ではなく行で切る
        chunks = ([ln for ln in raw_text.splitlines() if ln.strip()]
                  if f == "REFERENCES.md" else raw_text.split("\n\n"))
        for raw in chunks:
            s = raw.strip()
            if not s or s.startswith("$$"):
                continue
            if s.startswith("#"):
                s = re.sub(r"^#+\s*", "", s)
            out.append(s)
    return out


def table_xml(md_rows: list[str], template: str) -> str:
    """Markdown の pipe 表を、原典の表の書式を雛形にして組み立てる。"""
    rows = [r for r in md_rows if not re.fullmatch(r"\|[\s\-|:]+\|", r.strip())]
    cells = [[c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", r.strip().strip("|"))]
             for r in rows]
    n = max(len(c) for c in cells)
    pr = re.search(r"<w:tblPr>.*?</w:tblPr>", template, re.S)
    pr = pr.group(0) if pr else "<w:tblPr><w:tblW w:w=\"0\" w:type=\"auto\"/></w:tblPr>"
    grid = "<w:tblGrid>" + "<w:gridCol/>" * n + "</w:tblGrid>"
    body = ""
    for r, row in enumerate(cells):
        row = row + [""] * (n - len(row))
        tcs = ""
        for c in row:
            rpr = "<w:rPr><w:b/></w:rPr>" if r == 0 else ""
            tcs += ('<w:tc><w:tcPr><w:tcW w:w="0" w:type="auto"/></w:tcPr>'
                    f'<w:p><w:r>{rpr}<w:t xml:space="preserve">{D.esc(c)}</w:t></w:r></w:p></w:tc>')
        body += f"<w:tr>{tcs}</w:tr>"
    return f"<w:tbl>{pr}{grid}{body}</w:tbl>"


def main() -> int:
    xml = zipfile.ZipFile(SRC).read("word/document.xml").decode("utf-8")
    xml = D.strip_comments(xml)
    blocks = D.blocks_of(xml)
    end = next(i for i, b in enumerate(blocks) if D.text_of(b).strip() == "Supplementary Materials")
    main_blocks = blocks[:end]

    md = md_blocks()
    md_tables = {}                       # キャプション先頭 → 表の行
    md_text = []
    i = 0
    while i < len(md):
        s = md[i]
        if s.lstrip().startswith("|"):
            key = norm(md_text[-1])[:40] if md_text else f"t{i}"
            md_tables.setdefault(key, []).extend(s.splitlines())
        else:
            md_text.append(s)
        i += 1

    # 原典の本文段落と Markdown 段落を対応づける
    orig_idx = [i for i, b in enumerate(main_blocks)
                if not b.startswith("<w:tbl") and D.text_of(b) and "<m:oMath" not in b]
    a_seq = [norm(D.text_of(main_blocks[i]))[:120] for i in orig_idx]
    b_seq = [norm(s)[:120] for s in md_text]
    sm = difflib.SequenceMatcher(None, a_seq, b_seq, autojunk=False)
    pair, used = {}, set()
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ("equal", "replace"):
            for k in range(min(i2 - i1, j2 - j1)):
                pair[orig_idx[i1 + k]] = j1 + k
                used.add(j1 + k)

    out, cursor = [], 0
    for bi, blk in enumerate(main_blocks):
        if blk.startswith("<w:tbl"):
            prev = next((D.text_of(main_blocks[j]) for j in range(bi - 1, -1, -1)
                         if D.text_of(main_blocks[j])), "")
            key = norm(prev)[:40]
            if key not in md_tables:
                # キャプションが変わっている場合は、見出し行の先頭セルで探す
                m1 = re.search(r"<w:tr\b.*?</w:tr>", blk, re.S)
                cells0 = ([D.text_of(c) for c in
                           re.findall(r"<w:tc>.*?</w:tc>", m1.group(0), re.S)] if m1 else [])
                head0 = norm(cells0[0]) if cells0 else ""
                for k, rows in md_tables.items():
                    if not rows:
                        continue
                    mdcells = [norm(c) for c in rows[0].strip().strip("|").split("|")]
                    if head0 and mdcells and head0 == mdcells[0]:
                        key = k
                        break
            if key in md_tables:
                out.append(table_xml(md_tables.pop(key), blk))
            else:
                out.append(blk)
            continue
        if "<m:oMath" in blk:
            fixed = blk
            for a, b in SPELL:
                fixed = re.sub(rf"(<w:t(?: [^>]*)?>)([^<]*)",
                               lambda m: m.group(1) + m.group(2).replace(a, b), fixed)
            out.append(fixed)
            continue
        if "<w:drawing>" in blk:
            # 原典に貼られていた図は古い版なので外す。図は別ファイルで提出し、
            # キャプションは本文の該当箇所に置く。
            continue
        if not D.text_of(blk):
            out.append(blk)
            continue
        if bi in pair:
            j = pair[bi]
            # 追いついていない Markdown 段落を先に差し込む
            while cursor < j:
                if cursor not in used:
                    out.append(D.para(md_text[cursor]))
                cursor += 1
            out.append(D.set_text(blk, md_text[j]))
            cursor = j + 1
    while cursor < len(md_text):
        if cursor not in used:
            out.append(D.para(md_text[cursor]))
        cursor += 1

    out = place_figures(out)

    sect = re.search(r"<w:sectPr[ >].*?</w:sectPr>", xml, re.S)
    sect_xml = sect.group(0) if sect else ""
    # 原典の sectPr には既に lnNumType が入っている。無い原典から作り直した場合に備えた保険。
    if sect_xml and "<w:lnNumType" not in sect_xml:
        LN = '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="360"/>'
        m_cols = re.search(r"<w:cols[^>]*/>", sect_xml)
        if m_cols:                       # cols の直後が正しい並び順
            sect_xml = sect_xml.replace(m_cols.group(0), m_cols.group(0) + LN, 1)
        else:
            sect_xml = sect_xml.replace("</w:sectPr>", LN + "</w:sectPr>")
    body = "".join(out) + sect_xml
    head = xml[:xml.index("<w:body>") + len("<w:body>")]
    new = head + body + "</w:body></w:document>"

    dest = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Main_manuscript.docx"
    z = zipfile.ZipFile(SRC)
    media = tuple(n for n in z.namelist() if n.startswith("word/media/"))
    rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
    rels = re.sub(r'<Relationship[^>]*Target="media/[^"]*"[^>]*/>', "", rels)
    D.write_docx(SRC, dest, new,
                 drop=("word/comments.xml", "word/commentsExtended.xml",
                       "word/commentsIds.xml") + media,
                 extra={"docProps/core.xml": clean_core(z),
                        "word/_rels/document.xml.rels": rels})
    om, omp = D.n_math(new)
    print(f"  {dest}")
    print(f"  数式 oMath {om}（うちブロック {omp}）／ ブロック {len(out)}")
    print(f"  表 {sum(1 for b in out if b.startswith('<w:tbl'))} 件 "
          f"／ 未使用の Markdown 表 {list(md_tables)[:3]}")
    return 0


def place_figures(out: list[str]) -> list[str]:
    """図のキャプションを、その図への最後の言及の直後に移す。

    投稿時は図を別ファイルで出すが、キャプションは本文の該当箇所に置く。
    """
    caps = {}
    for i, b in enumerate(out):
        m = re.match(r"Figure (\d)\.", D.text_of(b) or "")
        if m:
            caps[int(m.group(1))] = i
    if not caps:
        return out
    keep = [b for i, b in enumerate(out) if i not in set(caps.values())]
    for n in sorted(caps, reverse=True):
        cap = out[caps[n]]
        # 図は Results のものなので、探す範囲は 3. Discussion の手前まで。
        stop = next((i for i, b in enumerate(keep)
                     if (D.text_of(b) or "").startswith("3. Discussion")), len(keep))

        # パネル付きの言及（Figure 1D など）を優先する。無ければ素の言及。
        def hits(pat):
            return [i for i, b in enumerate(keep[:stop])
                    if re.search(pat, D.text_of(b) or "")
                    and not re.match(rf"Figure {n}\.", D.text_of(b) or "")]
        panel = hits(rf"Figure {n}[A-D]")
        last = max(panel) if panel else max(hits(rf"Figure {n}\b"), default=None)
        if last is None:
            keep.append(cap)
        else:
            keep.insert(last + 1, cap)
    return keep


def clean_core(z) -> str:
    core = z.read("docProps/core.xml").decode("utf-8")
    core = re.sub(r"<dc:creator>.*?</dc:creator>", "<dc:creator></dc:creator>", core)
    core = re.sub(r"<cp:lastModifiedBy>.*?</cp:lastModifiedBy>",
                  "<cp:lastModifiedBy></cp:lastModifiedBy>", core)
    core = re.sub(r"<dc:subject>.*?</dc:subject>", "<dc:subject></dc:subject>", core)
    core = re.sub(r"<dc:description>.*?</dc:description>", "<dc:description></dc:description>", core)
    return core


if __name__ == "__main__":
    sys.exit(main())
