# -*- coding: utf-8 -*-
"""BI: 手法の原著などを文献に加える。

書誌は **すべて Crossref の記録から作る**（記憶で書かない）。DOI だけを手で与え、
著者・題名・巻・号・頁・年は api.crossref.org の返り値をそのまま使う。
取得した生の記録は results/roundR/crossref_records.json に残す。

本文には [SSGSEA] のような名前のキーで引用を置き、renumber_refs_C.py が
初出順に番号を振り直す。
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
CACHE = HERE / "results" / "roundR" / "crossref_records.json"
MAIL = "vm15139.jitudou@gmail.com"

# キー → (DOI, 投稿先の書式の雑誌略記)。年・著者・題名・巻・頁は Crossref から取る。
REFS = {
    # BI1 必須: 手法の原著
    "SSGSEA":    ("10.1038/nature08460", "Nature"),
    "GSVA":      ("10.1186/1471-2105-14-7", "BMC Bioinform."),
    "SINGSCORE": ("10.1186/s12859-018-2435-4", "BMC Bioinform."),
    "PLAGE":     ("10.1186/1471-2105-6-225", "BMC Bioinform."),
    "GSEAPY":    ("10.1093/bioinformatics/btac757", "Bioinformatics"),
    # BI2 推奨
    "PAGE":      ("10.1186/1471-2105-6-144", "BMC Bioinform."),
}
# 照会はするが本文には入れない。BJ1 で Lin / Gilad を外し、BI3 は採用しないと決まった。
OPTIONAL = {
    "LIN":    ("10.1073/pnas.1413624111", "Proc. Natl. Acad. Sci. USA"),
    "GILAD":  ("10.12688/f1000research.6536.1", "F1000Research"),
    "GOEMAN": ("10.1093/bioinformatics/btm051", "Bioinformatics"),
    "TEUFEL": ("10.1053/j.gastro.2016.05.051", "Gastroenterology"),
    "ZHOU":   ("10.1681/ASN.0000000000000217", "J. Am. Soc. Nephrol."),
}
ADD = list(REFS)          # 本文に入れるキー（BI1 + PAGE）


def crossref(doi: str) -> dict:
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    req = urllib.request.Request(url, headers={"User-Agent": f"ckd-xspecies/1.0 (mailto:{MAIL})"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["message"]


def records(refresh=False) -> dict:
    all_refs = {**REFS, **OPTIONAL}
    cache = json.loads(CACHE.read_text()) if CACHE.exists() and not refresh else {}
    for k, (doi, _ab) in all_refs.items():
        if k not in cache or cache[k].get("DOI", "").lower() != doi.lower():
            cache[k] = crossref(doi)
            print(f"  Crossref 照会: {k} {doi}")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=1))
    return cache


def initials(given: str) -> str:
    """'David A.' → 'D.A.' / 'Seon-Young' → 'S.-Y.'"""
    out = []
    for part in re.split(r"[ ]+", given.strip()):
        for sub in re.split(r"(-)", part):
            if sub == "-":
                out.append("-")
            elif sub:
                out.append(sub[0].upper() + ".")
    return "".join(out)


def authors(msg: dict, limit=10) -> str:
    au = msg.get("author", [])
    names = [f"{a.get('family','').strip()}, {initials(a.get('given',''))}".rstrip(", ")
             for a in au]
    if len(names) > limit:
        return "; ".join(names[:limit]) + "; et al."
    return "; ".join(names)


def year(msg: dict) -> int:
    """印刷版の年を優先する（オンライン先行で年がずれるため）。"""
    for f in ("published-print", "published", "published-online", "issued"):
        d = msg.get(f, {}).get("date-parts") or []
        if d and d[0] and d[0][0]:
            return int(d[0][0])
    raise SystemExit(f"年が取れない: {msg.get('DOI')}")


def pages(msg: dict) -> str:
    p = msg.get("page") or msg.get("article-number") or ""
    return p.replace("-", "–")


def entry(key: str, msg: dict) -> str:
    abbrev = {**REFS, **OPTIONAL}[key][1]
    title = (msg.get("title") or [""])[0].strip().rstrip(".")
    vol = msg.get("volume", "")
    bits = [authors(msg), f"{title}.", f"{abbrev} {year(msg)},", f"{vol}," if vol else "",
            f"{pages(msg)}." if pages(msg) else ""]
    s = " ".join(b for b in bits if b)
    return f"{s} https://doi.org/{msg['DOI']}."


# ----------------------------------------------------------------- 本文の編集
EDITS = [
    # BI1: 4 手法の初出（Results 2.1）
    ("RESULTS.md",
     "with the Bioconductor implementations of ssGSEA, GSVA, singscore and PLAGE "
     "(Section 4.14, Table S4)",
     "with the Bioconductor implementations of ssGSEA [SSGSEA], GSVA [GSVA], "
     "singscore [SINGSCORE] and PLAGE [PLAGE] (Section 4.14, Table S4)"),
    # BI1: §4.14 でもう一度、GSEApy の初出もここ
    ("METHODS.md",
     "GSVA 2.6.6 for ssGSEA, GSVA and PLAGE, and singscore 1.32.0 for singscore, "
     "with GSEABase 1.74.0 reading the gene sets.",
     "GSVA 2.6.6 for ssGSEA [SSGSEA], GSVA [GSVA] and PLAGE [PLAGE], and singscore "
     "1.32.0 for singscore [SINGSCORE], with GSEABase 1.74.0 reading the gene sets."),
    ("METHODS.md",
     "SciPy 1.13.1 and GSEApy 1.3.1.",
     "SciPy 1.13.1 and GSEApy 1.3.1 [GSEAPY]."),
    # BI2: PAGE を、中心化を定義する Methods（4.13）で
    ("METHODS.md",
     "The centered version standardized Δ within each state across the common gene "
     "universe before taking those means:",
     "The centered version standardized Δ within each state across the common gene "
     "universe before taking those means, so that each pathway mean is taken after "
     "the state-wide mean has been removed, as in the numerator of the PAGE "
     "statistic [PAGE]:"),
    # BI2: PAGE を §3.4 で
    ("DISCUSSION.md",
     "The arithmetic pathway means of the main analysis gave the same 0.356;",
     "A pathway mean taken after the state-wide mean has been removed is the "
     "numerator of the PAGE statistic [PAGE], so the centered form is the one that "
     "parametric gene-set scores of that family already use. The arithmetic pathway "
     "means of the main analysis gave the same 0.356;"),
]

# BJ1: Limitations は増やさない。BI2 で入れた Lin / Gilad の文を元に戻す。
# 第 2 要素が本文に残っていれば第 3 要素へ戻す（番号に振り直されていても拾えるよう
# 正規表現で書く）。
REVERTS = [
    ("DISCUSSION.md",
     r"Species and dataset are confounded, so the retained component cannot be "
     r"assigned to disease biology rather than to measurement; cross-species "
     r"expression comparisons have previously grouped by species in a way that was "
     r"later attributed to the confounding of species with processing batch "
     r"\[(?:LIN,GILAD|\d+[,–-]\d+)\],",
     "Species and dataset are confounded, so the retained component cannot be "
     "assigned to disease biology rather than to measurement,"),
]


def main() -> int:
    msgs = records("--refresh" in sys.argv)

    # 1. 書誌を作る（Crossref の記録からのみ）
    print("### Crossref から作った書誌")
    built = {}
    for k in {**REFS, **OPTIONAL}:
        built[k] = entry(k, msgs[k])
        mark = "  " if k in ADD else "（任意）"
        print(f"{mark}[{k}] {built[k]}")

    # 2. 本文に引用を置く（まず BJ1 の取り消し）
    print("\n### 本文の編集")
    for f, rx, back in REVERTS:
        q = M / f
        t = q.read_text()
        n = len(re.findall(rx, t))
        if n:
            q.write_text(re.sub(rx, back, t))
            print(f"  戻 {f}: Lin / Gilad の文を取り消した（{n} 箇所）")
    for f, old, new in EDITS:
        p = M / f
        t = p.read_text()
        # すでに当たっているか。番号に振り直されたあとでも拾えるように、
        # [KEY] の部分だけ「キーまたは数字」にゆるめた正規表現で見る。
        rx = re.escape(new)
        rx = re.sub(r"\\\[([A-Z,]+)\\\]",
                    lambda m: r"\[(?:" + re.escape(m.group(1)) + r"|[\d,–-]+)\]", rx)
        if re.search(rx, t):
            print(f"  済 {f}: {new[:52]}…")
            continue
        if t.count(old) != 1:
            raise SystemExit(f"✗ {f} に {t.count(old)} 件: {old[:70]}")
        p.write_text(t.replace(old, new))
        print(f"  入 {f}: {old[:52]}…")

    # 3. 文献一覧にキーで足す
    R = M / "REFERENCES.md"
    txt = R.read_text()
    # 採用しない文献が一覧に残っていれば外す（キーのままでも番号つきでも）
    dropped = []
    for k in OPTIONAL:
        doi = OPTIONAL[k][0].lower()
        keep = []
        for line in txt.splitlines():
            if doi in line.lower() and re.match(r"^(?:\d+|[A-Z]+)\.\s", line.strip()):
                dropped.append(k)
                continue
            keep.append(line)
        txt = "\n".join(keep)
    if dropped:
        print(f"文献一覧から外した（採用しない）: {dropped}")
    added = []
    for k in ADD:
        doi = REFS[k][0].lower()
        if doi in txt.lower():          # 番号に振り直されたあとは DOI で見る
            continue
        txt = txt.rstrip() + f"\n{k}. {built[k]}\n"
        added.append(k)
    R.write_text(txt)
    print(f"\n文献一覧に追加: {added or 'なし（すでにある）'}")
    print("  renumber_refs_C.py を回して番号を振り直すこと")
    return 0


if __name__ == "__main__":
    sys.exit(main())
