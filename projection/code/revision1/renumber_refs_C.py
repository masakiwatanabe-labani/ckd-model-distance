# -*- coding: utf-8 -*-
"""C案の文献番号を本文の初出順に振り直す。

新しい文献を 2 件加えたので、[ENSEMBL] / [ERCB] という仮のキーで本文に置いてある。
本文（フロント → 序論 → 結果 → 考察 → 方法 → 結論 → 背表紙 → 図 → 補足）の順に
引用を走査し、最初に出た順で 1 から番号を振り直す。範囲表記 [28-30] も展開して扱う。
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
ORDER = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
         "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md", "SUPPLEMENTARY.md"]
# BI で足した手法の原著などは名前のキーで本文に置く。長いものを先に並べる。
NAMED = ("ENSEMBL", "ERCB", "MODELS[123]", "SSGSEA", "SINGSCORE", "GSEAPY",
         "GSVA", "PLAGE", "PAGE", "GILAD", "GOEMAN", "TEUFEL", "ZHOU", "LIN")
KEYS = r"\d{1,2}|" + "|".join(NAMED)
CITE = re.compile(rf"\[((?:{KEYS})(?:\s*[,–-]\s*(?:{KEYS}))*)\]")
DASH = "–—-"


def expand(token: str) -> list[str]:
    """[21-23,26] のような中身をキーの列に開く。"""
    out = []
    for part in re.split(r"\s*,\s*", token.strip()):
        m = re.fullmatch(rf"(\d+)\s*[{DASH}]\s*(\d+)", part)
        if m:
            out += [str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1)]
        else:
            out.append(part.strip())
    return out


def collapse(nums: list[int]) -> str:
    """連番は範囲にまとめる。"""
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:                      # 3 件以上の連番だけ範囲にまとめる
            out.append(f"{nums[i]}–{nums[j]}")
        else:                              # 1 件または 2 件はそのまま並べる
            out += [str(n) for n in nums[i:j + 1]]
        i = j + 1
    return ",".join(out)


def main() -> int:
    # 1. 初出順にキーを集める
    order, seen = [], set()
    for f in ORDER:
        for m in CITE.finditer((M / f).read_text()):
            for k in expand(m.group(1)):
                if k not in seen:
                    seen.add(k); order.append(k)
    new = {k: i + 1 for i, k in enumerate(order)}
    print(f"本文に現れる文献 {len(order)} 件")

    # 2. 旧リストを読む
    refs = {}
    reftxt = (M / "REFERENCES.md").read_text()
    for line in reftxt.splitlines():
        m = re.match(rf"^(\d{{1,2}}|{'|'.join(NAMED)})\.\s+(.*)$", line.strip())
        if m:
            refs[m.group(1)] = m.group(2)
    missing = [k for k in order if k not in refs]
    unused = [k for k in refs if k not in new]
    assert not missing, f"本文にあってリストに無い: {missing}"
    if unused:
        print(f"  一度も引用されない文献（そのまま残す）: {sorted(unused)}")

    # 3. 本文の引用を書き換える
    def sub(m):
        return "[" + collapse([new[k] for k in expand(m.group(1))]) + "]"
    for f in ORDER:
        p = M / f
        p.write_text(CITE.sub(sub, p.read_text()))

    # 4. リストを新しい順に書き直す
    lines = ["# References", ""]
    for k in order:
        lines.append(f"{new[k]}. {refs[k]}")
    (M / "REFERENCES.md").write_text("\n".join(lines) + "\n")

    for k in order:
        if not k.isdigit() or int(k) >= 30:
            print(f"  {k:>8s} → [{new[k]}]  {refs[k][:56]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
