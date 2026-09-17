"""Part 4. マッチした遺伝子空間を、保存された 1 つの対応表から定義し直す。

これまでの食い違いの原因は、公表済みの aa_groupB_matched.txt を作った 1:1 マッチングが
乱数順で走っており、**対応表そのものが保存されていなかった**ことにある。そのため
  - Section 4.14 の「2,016 のうち 1,618 が相手を見つけた」は Group B 側の大きさであって
    Group A 側の本数ではない
  - 方向を変えて作り直すと件数が変わる（A→B で 1,632、B→A で 1,578、対称制限で 1,509）
という状態になっていた。

ここでは乱数を使わない 1 回のマッチングを対応表ごと保存し、matched Group A と
matched Group B の両方をその対応表から定義する。両空間は同じ大きさになり、SMD も
その対応表から計算される。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
from check1_geneset import joint_expression_rank  # noqa: E402

OUT = HERE / "results" / "round7"
OUT.mkdir(parents=True, exist_ok=True)
CALIPER = 0.02


def greedy(src, tgt, rank, caliper=CALIPER):
    """src の各遺伝子に、未使用の tgt から最も近いものを 1 対 1 で割り当てる（決定的）。"""
    tgt_sorted = sorted(tgt, key=lambda g: (float(rank[g]), g))
    tv = np.array([float(rank[g]) for g in tgt_sorted])
    used = np.zeros(len(tv), bool)
    pairs = []
    for g in sorted(src, key=lambda x: (float(rank[x]), x)):
        t = float(rank[g])
        pos = int(np.searchsorted(tv, t))
        best, bestd = -1, caliper
        l, r = pos - 1, pos
        while True:
            moved = False
            if r < len(tv) and (tv[r] - t) <= bestd:
                moved = True
                if not used[r] and (tv[r] - t) < bestd:
                    best, bestd = r, tv[r] - t
                r += 1
            if l >= 0 and (t - tv[l]) <= bestd:
                moved = True
                if not used[l] and (t - tv[l]) < bestd:
                    best, bestd = l, t - tv[l]
                l -= 1
            if not moved:
                break
        if best >= 0:
            used[best] = True
            pairs.append((g, tgt_sorted[best], float(tv[best])))
    return pairs


def smd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float((a.mean() - b.mean()) / np.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1))))


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0).dropna(axis=0, how="any")
    expr = pd.read_csv(HERE / "control_log2cpm.tsv", sep="\t", index_col=0)
    gA = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    jr = joint_expression_rank(expr, D.index).dropna()
    inA = sorted(g for g in jr.index if g in gA)
    inB = sorted(g for g in jr.index if g not in gA)

    pr = greedy(inA, inB, jr)
    P = pd.DataFrame(pr, columns=["groupA_gene", "groupB_gene", "groupB_expr_rank"])
    P["groupA_expr_rank"] = [float(jr[g]) for g in P.groupA_gene]
    P = P[["groupA_gene", "groupB_gene", "groupA_expr_rank", "groupB_expr_rank"]]
    P.round(6).to_csv(OUT / "matched_pairs.tsv", sep="\t", index=False)

    selA = sorted(P.groupA_gene)
    selB = sorted(P.groupB_gene)
    (HERE / "aa_groupA_matched.txt").write_text("\n".join(selA) + "\n")
    (HERE / "aa_groupB_matched2.txt").write_text("\n".join(selB) + "\n")

    un = [g for g in inA if g not in set(selA)]
    before = smd(jr[inA], jr[inB])
    after = smd(jr[selA], jr[selB])
    lines = [
        ("complete case genes", len(jr)),
        ("Group A genes", len(inA)),
        ("Group B genes", len(inB)),
        ("matched pairs", len(P)),
        ("Group A genes with no partner", len(un)),
        ("median expression rank of those", round(float(np.median([float(jr[g]) for g in un])), 4)),
        ("median expression rank, matched Group A", round(float(P.groupA_expr_rank.median()), 4)),
        ("median expression rank, matched Group B", round(float(P.groupB_expr_rank.median()), 4)),
        ("SMD, all Group A vs all Group B", round(before, 4)),
        ("SMD, matched Group A vs matched Group B", round(after, 4)),
        ("largest absolute rank difference within a pair",
         round(float((P.groupA_expr_rank - P.groupB_expr_rank).abs().max()), 6)),
    ]
    S = pd.DataFrame(lines, columns=["quantity", "value"])
    S.to_csv(OUT / "matched_spaces_summary.tsv", sep="\t", index=False)
    print(S.to_string(index=False))
    print(f"\n書き出し: {OUT / 'matched_pairs.tsv'}, {HERE / 'aa_groupA_matched.txt'}, "
          f"{HERE / 'aa_groupB_matched2.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
