# -*- coding: utf-8 -*-
"""第9ラウンドの提出物を組み立てる。

~/Desktop/CKD_xspecies_round9/
  manuscript/       更新後の Markdown 一式
  extracted/        abstract.txt, captions.md, table1〜5.tsv, paragraph_order.tsv
"""
from __future__ import annotations
import os
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
MD = HERE / "manuscript"
DEST = Path(os.path.expanduser("~/Desktop/CKD_xspecies_round9"))


def main() -> int:
    (DEST / "extracted").mkdir(parents=True, exist_ok=True)
    if (DEST / "manuscript").exists():
        shutil.rmtree(DEST / "manuscript")
    shutil.copytree(MD, DEST / "manuscript",
                    ignore=shutil.ignore_patterns("figures", "NUMBERS.md", "*.docx"))

    # ---- abstract.txt
    fm = (MD / "FRONTMATTER.md").read_text()
    m = re.search(r"## Abstract\n+(.+?)\n\n\*\*Keywords", fm, re.S)
    assert m, "Abstract が見つからない"
    abst = " ".join(m.group(1).split())
    (DEST / "extracted" / "abstract.txt").write_text(abst + "\n")

    # ---- captions.md（本文の出現順に、キャプション段落だけを拾う）
    caps = []
    for f in ("RESULTS.md", "METHODS.md"):
        for p in (MD / f).read_text().split("\n\n"):
            if re.match(r"\*\*(Figure|Table) S?\d", p.strip()):
                caps.append(" ".join(p.split()))
    for p in (MD / "SUPPLEMENTARY.md").read_text().split("\n\n"):
        if re.match(r"\*\*(Figure|Table) S?\d", p.strip()):
            caps.append(" ".join(p.split()))
    (DEST / "extracted" / "captions.md").write_text(
        "# Figure and table captions, final text after round 9\n\n---\n\n"
        + "\n\n---\n\n".join(caps) + "\n")

    # ---- table1〜5.tsv（pipe 表をそのまま TSV に）
    def pipe_tables(text):
        """キャプション直後に続く pipe 表のブロックだけを返す。

        キャプション本文を読み飛ばし、最初の pipe 行から数える。地の文にも `|` で
        囲まれた数式（|log2(...)| など）が出るので、表が始まったあとで pipe 以外の
        非空行に当たったらそこで打ち切る。
        """
        out, cur = [], []
        for line in text.splitlines():
            if line.lstrip().startswith("|"):
                cur.append(line)
            elif cur:
                out.append(cur); cur = []
        if cur:
            out.append(cur)
        # 区切り行 |---| を持つものだけが本物の表。地の文の |…| はここで落ちる。
        sep = re.compile(r"\|[\s\-|:]+\|")
        return [b for b in out if any(sep.fullmatch(x.strip()) for x in b)]

    def split_cells(line):
        """`\|` でエスケープされた縦棒はセルの区切りにしない。"""
        s = line.strip().strip("|")
        parts = re.split(r"(?<!\\)\|", s)
        return [c.strip().replace("\\|", "|") for c in parts]

    def to_tsv(block):
        rows = []
        for line in block:
            if re.fullmatch(r"\|[\s\-|:]+\|", line.strip()):
                continue
            rows.append("\t".join(split_cells(line)))
        return "\n".join(rows) + "\n"

    res, met = (MD / "RESULTS.md").read_text(), (MD / "METHODS.md").read_text()
    # Table 1〜4 は Results、Table 5 は Methods。キャプションの直後のブロックを取る。
    n = 0
    for text in (res, met):
        for m in re.finditer(r"\*\*Table (\d+)\.\*\*", text):
            num = int(m.group(1))
            rest = text[m.end():]
            blocks = pipe_tables(rest[:rest.find("\n\n**") if "\n\n**" in rest else len(rest)])
            if not blocks:
                continue
            body = "".join(to_tsv(b) for b in blocks)
            (DEST / "extracted" / f"table{num}.tsv").write_text(body)
            n += 1
    assert n >= 5, f"表が {n} 件しか取れていない"

    # ---- paragraph_order.tsv
    sys.path.insert(0, str(HERE / "code" / "revision1"))
    import paragraph_order  # noqa: E402
    paragraph_order.main()
    shutil.copy2(HERE / "extracted" / "paragraph_order.tsv",
                 DEST / "extracted" / "paragraph_order.tsv")

    print(f"提出物: {DEST}")
    for p in sorted((DEST / "extracted").iterdir()):
        print(f"  extracted/{p.name} ({p.stat().st_size} bytes)")
    print(f"  manuscript/ ({len(list((DEST/'manuscript').iterdir()))} files)")
    print(f"  Abstract {len(abst.split())} 語 / キャプション {len(caps)} 件")
    return 0


if __name__ == "__main__":
    sys.exit(main())
