"""Part 4. Table 4 の matched 空間の遺伝子数と SMD が食い違う原因を突き止める。

食い違い:
  - §4.14 は「2,016 の Group A のうち 1,618 が相手を見つけた」(除外 398)
  - Table 4 の matched Group A は 1,578 遺伝子、除外 438
両者は別方向のマッチングである可能性が高い。ここで両方向を決定的に作り直し、
公表済みの aa_groupB_matched.txt（1,618）と突き合わせる。
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
    """src の各遺伝子に、未使用の tgt から最も近いものを 1 対 1 で割り当てる。

    処理順は src の順位の小さい順（決定的）。同点は遺伝子名の辞書順。
    """
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
            pairs.append((g, tgt_sorted[best]))
    return pairs


def smd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float((a.mean() - b.mean()) / np.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1))))


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0).dropna(axis=0, how="any")
    expr = pd.read_csv(HERE / "control_log2cpm.tsv", sep="\t", index_col=0)
    gA = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    pubB = [g.strip() for g in (HERE / "aa_groupB_matched.txt").read_text().split() if g.strip()]

    jr = joint_expression_rank(expr, D.index).dropna()
    inA = sorted(g for g in jr.index if g in gA)
    inB = sorted(g for g in jr.index if g not in gA)
    pubB = sorted(g for g in pubB if g in jr.index)
    print(f"complete case {len(jr)}: Group A {len(inA)}, Group B {len(inB)}, "
          f"published matched B {len(pubB)}")

    rows = []
    a2b = greedy(inA, inB, jr)
    b2a = greedy(pubB, inA, jr)
    sym = greedy(inA, pubB, jr)
    for lab, pr, sa, sb in [
            ("A -> all B (the direction Section 4.14 describes)", a2b, "Group A", "Group B"),
            ("published matched B -> A (used for Table 4 in round 5)", b2a, "matched B", "Group A"),
            ("A -> published matched B (symmetric restriction)", sym, "Group A", "matched B")]:
        left = [p[0] for p in pr]
        right = [p[1] for p in pr]
        rows.append(dict(matching=lab, n_pairs=len(pr),
                         left_set=sa, right_set=sb,
                         left_unmatched=(len(inA) if sa.startswith("Group A") else len(pubB)) - len(pr),
                         smd_before=round(smd(jr[inA], jr[pubB]), 4),
                         smd_after=round(smd(jr[left], jr[right]), 4),
                         right_equals_published_B=(set(right) == set(pubB))))
        print(f"\n{lab}\n  pairs {len(pr)}  SMD after {rows[-1]['smd_after']:+.4f}"
              f"  right side == published matched B: {rows[-1]['right_equals_published_B']}"
              f"  overlap with published B: {len(set(right) & set(pubB))}")

    T = pd.DataFrame(rows)
    T.to_csv(OUT / "matching_audit.tsv", sep="\t", index=False)

    # 採用する対応表: 両側を公表済み matched B に対して対称に制限したもの
    P = pd.DataFrame(sym, columns=["groupA_gene", "groupB_gene"])
    P["groupA_expr_rank"] = [float(jr[g]) for g in P.groupA_gene]
    P["groupB_expr_rank"] = [float(jr[g]) for g in P.groupB_gene]
    P.round(6).to_csv(OUT / "matched_pairs.tsv", sep="\t", index=False)
    (HERE / "aa_groupA_matched.txt").write_text("\n".join(sorted(P.groupA_gene)) + "\n")
    (HERE / "aa_groupB_matched_sym.txt").write_text("\n".join(sorted(P.groupB_gene)) + "\n")
    un = [g for g in inA if g not in set(P.groupA_gene)]
    print(f"\n採用: 対称制限 {len(P)} 組")
    print(f"  Group A 側で対応がつかない {len(un)} 遺伝子の発現量順位 中央 "
          f"{np.median([float(jr[g]) for g in un]):.4f}")
    print(f"  matched B 側で対応がつかない {len(pubB) - len(P)} 遺伝子")
    print(f"  SMD 全 Group A vs 公表 matched B {smd(jr[inA], jr[pubB]):+.4f} "
          f"→ 対称制限後 {smd(jr[P.groupA_gene], jr[P.groupB_gene]):+.4f}")
    print(f"\n書き出し: {OUT / 'matching_audit.tsv'}, {OUT / 'matched_pairs.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
