"""Table 1-4 の全セルを、図の印字値と同じ形式で TSV に落とす。

図は make_figures.py が描画時に rec() で記録するが、表は Markdown の pipe 表が
そのまま版面になるので、pipe 表から読む。列は図側の figure_printed_numbers.tsv に
合わせて (table, row, column, printed) とし、verify_numbers.py が本文と突き合わせる。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
MD = HERE / "manuscript" / "RESULTS.md"
OUT = HERE / "results" / "round5" / "table_printed_numbers.tsv"


def cells(text: str):
    """**Table N.** のキャプションから次のキャプションまでの pipe 表をすべて返す。

    Table 2 は一つのキャプションの下に 2 ブロックあるので、ブロックごとに
    先頭行を見出しとして扱い、(表番号, ブロック番号, 行) を返す。
    """
    marks = [(m.start(), m.group(1)) for m in re.finditer(r"\*\*Table (\d+)\.\*\*", text)]
    stops = sorted(m.start() for m in re.finditer(r"\*\*(?:Table|Figure) ", text))
    for pos, n in marks:
        end = next((s for s in stops if s > pos), len(text))
        span = text[pos:end]
        block, rows = 0, []
        for line in span.split("\n"):
            st = line.strip()
            if st.startswith("|"):
                rows.append(st)
            elif rows:
                block += 1
                yield n, block, rows
                rows = []
        if rows:
            yield n, block + 1, rows


def main() -> int:
    text = MD.read_text()
    out = [("table", "row", "column", "printed")]
    ndata = {}
    for n, block, rows in cells(text):
        parsed = []
        for raw in rows:
            cols = [c.strip() for c in raw.strip("|").split("|")]
            if all(set(c) <= set("-: ") for c in cols):      # 区切り行
                continue
            parsed.append(cols)
        if not parsed:
            continue
        head, body = parsed[0], parsed[1:]
        tag = f"Table{n}" if block == 1 else f"Table{n}b{block}"
        for j, c in enumerate(head):
            out.append((tag, "header", str(j), c))
        for cols in body:
            for j, c in enumerate(cols[1:], start=1):
                colname = head[j] if j < len(head) else str(j)
                out.append((tag, cols[0], colname, c))
        ndata[tag] = (len(body), len(head))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join("\t".join(r) for r in out) + "\n")
    print(OUT, f"({len(out) - 1} rows)")
    for tag in sorted(ndata):
        r, c = ndata[tag]
        print(f"  {tag}: {r} data rows x {c} columns")
    return 0


if __name__ == "__main__":
    sys.exit(main())
