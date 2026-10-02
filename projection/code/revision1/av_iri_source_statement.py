# -*- coding: utf-8 -*-
"""AV2: GSE98622 の単位について、原著が何と書いているかを記録する。

査読コメントや第三者の要約ではなく、原著（Liu et al., JCI Insight 2017;2:e94716,
GSE98622 の元論文）の本文から逐語で取る。Europe PMC の全文 XML を使う。
"""
from __future__ import annotations

import html
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
OUT = HERE / "results" / "roundR" / "iri_source_statement.tsv"
CACHE = HERE / "data" / "base" / "liu2017_jciinsight_fulltext.xml"
URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5612583/fullTextXML"
DOI = "10.1172/jci.insight.94716"


def fulltext() -> str:
    if CACHE.exists():
        return CACHE.read_text(errors="replace")
    with urllib.request.urlopen(URL, timeout=120) as r:
        x = r.read().decode("utf-8", "replace")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(x)
    return x


def main() -> int:
    flat = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", fulltext())))
    m = re.search(r"we analyzed RNA-seq data at the gene level.*?GSE98622\s*\)", flat)
    if not m:
        print("原著の該当文が見つからない", file=sys.stderr)
        return 1
    quote = re.sub(r"\s+", " ", m.group(0)).replace("f ragments p er k ilobase of transcript "
                                                    "per m illion mapped reads",
                                                    "fragments per kilobase of transcript per "
                                                    "million mapped reads")
    n_genes = re.search(r"from ([\d,]+) genes", quote)
    rows = [
        {"quantity": "source paper", "value": "Liu et al., JCI Insight 2017;2:e94716"},
        {"quantity": "doi", "value": DOI},
        {"quantity": "full text", "value": "Europe PMC PMC5612583"},
        {"quantity": "genes stated in the paper", "value": n_genes.group(1) if n_genes else ""},
        {"quantity": "verbatim statement", "value": quote},
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    for r in rows:
        print(f"  {r['quantity']:26s} {r['value'][:150]}")
    print(f"\n書き出し: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
