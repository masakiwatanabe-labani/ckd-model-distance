# -*- coding: utf-8 -*-
"""AY1: データセット・著者名・アクセッションを含む文の [n] が、その文献と合っているか。

投稿版以降に書いた文で「ネコのデータの出典」を [25]（Kriegeskorte）と引いていた種類の
誤りを機械的に捕まえる。文中の手がかり語と、引かれた文献の書誌を突き合わせる。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
FILES = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md", "RESPONSE.md", "COVER_LETTER.md"]
# 文中の手がかり → 引かれた文献の書誌に必ず含まれるべき語
CLUES = [
    (r"\bLi et al\.|42003_2025_9164|MOESM3", ("Li", "domestic cats")),
    (r"\bLiu et al\.|GSE98622", ("Liu", "ischemia")),
    (r"\bEnsembl\b|BioMart", ("Ensembl", "Harrison")),
    (r"\bPXD066590\b", ("Li", "domestic cats")),
    (r"\bGSE303653\b", ("Li", "domestic cats")),
    (r"\bGSE299326\b|Pod-TRECK model", ("Watanabe", "diphtheria")),
    (r"\bGSE104954\b", ("Ju", "Nair", "tubulointerstitial", "cDNA")),
    # BI: 手法の原著。手法名を出す文が文献を引くなら、その原著であること
    (r"\bssGSEA\b", ("Barbie", "TBK1")),
    (r"\bGSVA\b", ("Hänzelmann", "gene set variation")),
    (r"\bsingscore\b", ("Foroutan", "Single sample scoring")),
    (r"\bPLAGE\b", ("Tomfohr", "singular value decomposition")),
    (r"\bGSEApy\b", ("Fang", "GSEApy")),
    (r"\bPAGE\b", ("Kim", "Parametric Analysis")),
    (r"confounding of species with processing batch|mouse ENCODE",
     ("Lin", "Gilad", "ENCODE")),
]


def refs() -> dict:
    out = {}
    for line in (M / "REFERENCES.md").read_text().splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m:
            out[m.group(1)] = m.group(2)
    return out


def main() -> int:
    R = refs()
    bad, nocite, n = [], [], 0
    for f in FILES:
        p = M / f
        if not p.exists():
            continue
        text = re.sub(r"> \*.*?\*", "", p.read_text())     # 査読票の引用は対象外
        # 段落ではなく文単位で見る。段落単位だと、同じ段落の別の文の引用で通ってしまう。
        sents = []
        for para in text.split("\n\n"):
            flat = " ".join(para.split())
            if flat.startswith("|"):
                continue
            sents += [s.strip() for s in re.split(r"(?<=[.;])\s+", flat) if s.strip()]
        for flat in sents:
            for pat, need in CLUES:
                if not re.search(pat, flat):
                    continue
                # その文（引用のある節）の中の [n]
                for mm in re.finditer(r"\[(\d+(?:[,–-]\d+)*)\]", flat):
                    for tok in re.split(r"[,–-]", mm.group(1)):
                        if tok not in R:
                            continue
                        n += 1
                        if any(w.lower() in R[tok].lower() for w in need):
                            break
                    else:
                        continue
                    break
                else:
                    if re.search(r"\[\d", flat):
                        cited = re.findall(r"\[(\d+(?:[,–-]\d+)*)\]", flat)
                        bad.append(f"{f}: {need[0]} ら を指すべき文の引用が {cited}: "
                                   f"{flat[:80]}")
                    else:
                        nocite.append(f"{f}: {need[0]} ら に触れて引用が無い: {flat[:76]}")
    print(f"データセット・著者名・アクセッションを含む文の引用 {n} 件を突き合わせた（文単位）")
    if nocite:
        print(f"  参考: 出典に触れているが引用の無い文 {len(nocite)} 件（誤りではない。目視用）")
        for s in nocite[:10]:
            print("    - " + s)
    if bad:
        for b in bad[:8]:
            print("  ✗ " + b)
        return 1
    print("  すべて書誌と合っている")
    return 0


if __name__ == "__main__":
    sys.exit(main())
