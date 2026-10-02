# -*- coding: utf-8 -*-
"""C案の提出物を組み立てる。"""
from __future__ import annotations
import os
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
MD = HERE / "manuscript_C"
DEST = Path(os.path.expanduser("~/Desktop/CKD_xspecies_C_final"))
BODY = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
        "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md", "REFERENCES.md"]


def pipe_tables(text):
    out, cur = [], []
    for line in text.splitlines():
        if line.lstrip().startswith("|"):
            cur.append(line)
        elif cur:
            out.append(cur); cur = []
    if cur:
        out.append(cur)
    sep = re.compile(r"\|[\s\-|:]+\|")
    return [b for b in out if any(sep.fullmatch(x.strip()) for x in b)]


def to_tsv(block):
    rows = []
    for line in block:
        if re.fullmatch(r"\|[\s\-|:]+\|", line.strip()):
            continue
        cells = re.split(r"(?<!\\)\|", line.strip().strip("|"))
        rows.append("\t".join(c.strip().replace("\\|", "|") for c in cells))
    return "\n".join(rows) + "\n"


def main() -> int:
    for sub in ("manuscript", "extracted", "figures", "figures_png", "results", "code"):
        (DEST / sub).mkdir(parents=True, exist_ok=True)

    for f in BODY:
        shutil.copy2(MD / f, DEST / "manuscript" / f)

    # ---- abstract
    m = re.search(r"^Abstract:(.*?)$", (MD / "FRONTMATTER.md").read_text(), re.M)
    abst = " ".join(m.group(1).split())
    (DEST / "extracted" / "abstract.txt").write_text(abst + "\n")

    # ---- captions（本文の出現順）
    caps = []
    for f in ("RESULTS.md", "METHODS.md", "FIGURES.md", "SUPPLEMENTARY.md"):
        for p in (MD / f).read_text().split("\n\n"):
            if re.match(r"^(Table|Figure) S?\d+\.", p.strip()):
                caps.append(" ".join(p.split()))
    (DEST / "extracted" / "captions.md").write_text(
        "# Figure and table captions, C draft final\n\n---\n\n" + "\n\n---\n\n".join(caps) + "\n")

    # ---- table1..N
    n = 0
    for f in ("RESULTS.md", "METHODS.md"):
        text = (MD / f).read_text()
        for mm in re.finditer(r"^(Table (\d+))\.", text, re.M):
            num = int(mm.group(2))
            rest = text[mm.end():]
            stop = re.search(r"\n\n(?:Table|Figure) \d+\.", rest)
            blocks = pipe_tables(rest[:stop.start() if stop else len(rest)])
            if blocks:
                (DEST / "extracted" / f"table{num}.tsv").write_text(
                    "".join(to_tsv(b) for b in blocks))
                n += 1
    assert n >= 3, f"表が {n} 件"

    # ---- paragraph_order
    rows = ["seq\tfile\tsection\tkind\tfirst60"]
    seq = 0
    sec = "0"
    for f in BODY:
        for blk in (MD / f).read_text().split("\n\n"):
            s = blk.strip()
            if not s:
                continue
            seq += 1
            hm = re.match(r"^#+\s+([\d.]+)\.", s)
            if hm:
                sec = hm.group(1)
            kind = ("heading" if s.startswith("#") else "table" if s.startswith("|")
                    else "equation" if s.startswith("$$")
                    else "caption" if re.match(r"^(Table|Figure) S?\d+\.", s) else "text")
            flat = re.sub(r"\s+", " ", re.sub(r"[#|$*]", "", s)).strip()[:60]
            rows.append(f"{seq}\t{f}\t{sec}\t{kind}\t{flat}")
    (DEST / "extracted" / "paragraph_order.tsv").write_text("\n".join(rows) + "\n")

    # ---- figures
    for p in sorted((MD / "figures").glob("*.pdf")):
        shutil.copy2(p, DEST / "figures" / p.name)
    for p in sorted((MD / "figures").glob("*.png")):
        shutil.copy2(p, DEST / "figures_png" / p.name)

    # ---- results（このラウンドで新しく作ったもの）
    shutil.copy2(HERE / "results" / "roundC" / "centering_controls.tsv", DEST / "results")
    shutil.copy2(HERE / "results" / "revision1" / "pathway_reactome"
                 / "figure_printed_numbers_C.tsv", DEST / "results")
    shutil.copy2(HERE / "results" / "provenance.tsv", DEST / "results")

    # ---- code
    for f in ("centering_controls.py", "make_figures_C.py", "docx_to_markdown_C.py",
              "renumber_refs_C.py", "build_manuscript_C.sh", "apply_C_part2.py",
              "apply_C_part2b.py", "apply_C_parts.py", "apply_C_part5b.py",
              "apply_C_refs_eq.py", "apply_C_final_polish.py", "deliver_C.py"):
        shutil.copy2(HERE / "code" / "revision1" / f, DEST / "code" / f)
    shutil.copy2(HERE / "verify_numbers.py", DEST / "code" / "verify_numbers.py")

    print(f"提出物: {DEST}")
    for sub in ("manuscript", "extracted", "figures", "figures_png", "results", "code"):
        items = sorted(p.name for p in (DEST / sub).iterdir())
        print(f"  {sub}/ ({len(items)}): {', '.join(items[:6])}{' …' if len(items) > 6 else ''}")
    print(f"  Abstract {len(abst.split())} 語 / キャプション {len(caps)} 件 / 表 {n} 件")
    return 0


if __name__ == "__main__":
    sys.exit(main())
