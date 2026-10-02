# -*- coding: utf-8 -*-
"""AH4: 主張のトーンを全文で検査する（査読で求められたトーンの点検）。

機械的に削らない。ヒットを文脈つきで出し、判断した結果を ALLOW に登録していく。
ALLOW は「この用法は問題ない」と人が確認した箇所だけを、文の断片で持つ。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
FILES = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md"]

# 語彙。§3.4 で使っていたものに今回の点検向けを足す
TERMS = [
    r"recommend\w*", r"should be preferred", r"we suggest", r"best practice", r"guidelines?",
    r"superior\w*", r"outperform\w*", r"better reflects?\w*", r"more accurate\w*",
    r"improved?", r"improvement\w*", r"optimal\w*", r"validate[ds]? performance", r"fidelity",
    r"preferable", r"advantage\w*", r"gold standard", r"more faithful\w*", r"faithfully",
    r"better model", r"best model", r"model quality", r"translational value",
]
RX = re.compile("|".join(f"(?:{t})" for t in TERMS), re.I)

# 人が確認して残すと決めた用法。前後 40 字ていどの一意な断片で指定する
ALLOW = [
    # 「〜とは言えない」形で、主張を否定するために語を使っている箇所
    "does not establish that the uncentered score is biologically preferable",
    "do not establish a focal-adhesion-specific model advantage",
    "does not identify a preferred model time point",
    "no implication that either representation is preferable",
    "not a recommendation",
    "We draw no recommendation about choosing an earlier or later model time point",
    # 引用元の設計や既存文献の記述で、本研究の主張ではない箇所
    "improved genome assembly",
    # 以下は AH4 で一件ずつ文脈を確認した。いずれも語を否定の中で使っている
    "rather than a preferred species or a validated measure of disease fidelity",
    "Neither has been established here as a superior predictor of biological fidelity",
    "do not rank biological fidelity",
    "Neither result indicates which representation better reflects biological similarity",
    "use this analysis to rank species by human disease fidelity",
    "neither identifies a representation that is biologically more faithful",
    "They cannot be used to identify an overall best model",
]


def main() -> int:
    hits, allowed = [], 0
    for f in FILES:
        for i, line in enumerate((M / f).read_text().splitlines(), 1):
            for m in RX.finditer(line):
                a, b = max(0, m.start() - 90), min(len(line), m.end() + 90)
                ctx = line[a:b]
                if any(s in line for s in ALLOW):
                    if any(s in ctx for s in ALLOW):
                        allowed += 1
                        continue
                hits.append((f, i, m.group(0), ctx))
    print(f"語彙 {len(TERMS)} 件で全文を走査。確認済みとして許可 {allowed} 件 / 未確認 {len(hits)} 件")
    for f, i, w, ctx in hits:
        print(f"\n  ✗ {f}:{i}  「{w}」")
        print(f"     …{ctx}…")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
