"""Part F-1. 公表済みの matched Group B（1,618 遺伝子）に 1:1 で対応する Group A の
部分集合を作る。

Table 4 は Group A（2,016）と matched Group B（1,618）を並べているが、Group A 側には
対応相手のない 398 遺伝子（高発現側に偏る）が残っており、比較が対称でない。ここでは
matched Group B の各遺伝子に、対照群の発現量パーセンタイル順位が最も近い Group A 遺伝子を
caliper 0.02 以内で 1 対 1 に割り当て、同じ大きさの Group A 部分集合を作る。

割り当ては乱数を使わない。B 側の遺伝子を順位の小さい順に処理し、未使用の Group A 遺伝子の
うち最も近いものを取る（同点は遺伝子名の辞書順）。したがって何度走らせても同じ集合になる。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
from check1_geneset import joint_expression_rank  # noqa: E402

OUT_LIST = HERE / "aa_groupA_matched.txt"
OUT = HERE / "results" / "round5"
CALIPER = 0.02


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0).dropna(axis=0, how="any")
    expr = pd.read_csv(HERE / "control_log2cpm.tsv", sep="\t", index_col=0)
    gA = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    mB = [g.strip() for g in (HERE / "aa_groupB_matched.txt").read_text().split() if g.strip()]

    jr = joint_expression_rank(expr, D.index).dropna()
    inA = sorted([g for g in jr.index if g in gA])
    mB = sorted([g for g in mB if g in jr.index])
    print(f"complete case {len(jr)} 遺伝子 / Group A {len(inA)} / matched Group B {len(mB)}")

    a_rank = jr[inA].to_numpy()
    order = np.argsort(a_rank, kind="stable")          # 名前順に並べた上で順位で安定ソート
    a_sorted = [inA[i] for i in order]
    a_vals = a_rank[order]
    used = np.zeros(len(a_sorted), bool)

    pairs = []
    for g in sorted(mB, key=lambda x: (float(jr[x]), x)):
        t = float(jr[g])
        pos = int(np.searchsorted(a_vals, t))
        best, bestd = -1, CALIPER
        l, r = pos - 1, pos
        while True:
            moved = False
            if r < len(a_vals) and (a_vals[r] - t) <= bestd:
                moved = True
                if not used[r] and (a_vals[r] - t) < bestd:
                    best, bestd = r, a_vals[r] - t
                r += 1
            if l >= 0 and (t - a_vals[l]) <= bestd:
                moved = True
                if not used[l] and (t - a_vals[l]) < bestd:
                    best, bestd = l, t - a_vals[l]
                l -= 1
            if not moved:
                break
        if best >= 0:
            used[best] = True
            pairs.append((a_sorted[best], g, float(a_vals[best]), t))

    P = pd.DataFrame(pairs, columns=["groupA_gene", "groupB_gene",
                                     "groupA_expr_rank", "groupB_expr_rank"])
    selA = sorted(P.groupA_gene)
    OUT_LIST.write_text("\n".join(selA) + "\n")
    OUT.mkdir(parents=True, exist_ok=True)
    P.round(6).to_csv(OUT / "groupA_groupB_matched_pairs.tsv", sep="\t", index=False)

    da, db = P.groupA_expr_rank.to_numpy(), P.groupB_expr_rank.to_numpy()
    smd = (da.mean() - db.mean()) / np.sqrt(0.5 * (da.var(ddof=1) + db.var(ddof=1)))
    allA = jr[inA].to_numpy()
    allB = jr[mB].to_numpy()
    smd_before = (allA.mean() - allB.mean()) / np.sqrt(
        0.5 * (allA.var(ddof=1) + allB.var(ddof=1)))
    print(f"対応がついた組: {len(P)}")
    print(f"発現量順位 中央値: Group A {np.median(da):.4f} / Group B {np.median(db):.4f}")
    print(f"標準化平均差 SMD: 全 Group A vs matched Group B {smd_before:+.4f}"
          f"  →  対応した部分集合どうし {smd:+.4f}")
    print(f"対応がつかなかった Group A 遺伝子: {len(inA) - len(P)}"
          f"（発現量順位 中央値 "
          f"{np.median([float(jr[g]) for g in inA if g not in set(P.groupA_gene)]):.4f}）")
    print(f"\n書き出し: {OUT_LIST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
