"""Task 1. 天井までの距離 gap を3群で集計し、群間差を状態ラベル置換で検定する。

群の定義は Figure 2D と同一:
  A 同種・対照非共有 (39)
  B 異種・信頼性の高いネコ状態（cat_CKD34 / cat_med_CKD12）を含む (24)
  C 異種・信頼性の低いネコ状態（cat_CKD12 / cat_med_CKD34）を含む (24)

天井は補正なしと Spearman-Brown の両方。検定はペアが独立でないため、
状態ラベル（どの4状態をネコとするか）の網羅置換 C(16,4)=1820 で帰無分布を作る。
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
OUT = HERE / "results" / "revision1"
OUT.mkdir(parents=True, exist_ok=True)
RELIABLE = {"cat_CKD34", "cat_med_CKD12"}


def main():
    U = pd.read_csv(OUT / "pair_uncertainty_both_ceilings.tsv", sep="\t")
    U = U[(U.cls == "cross_species") | (~U.shared_control)].copy()
    U["cat"] = U.a.where(U.a.str.startswith("cat"), U.b)
    U["group"] = np.where(U.cls != "cross_species", "A same-species, no shared controls",
                          np.where(U.cat.isin(RELIABLE),
                                   "B cross-species, reliable feline state",
                                   "C cross-species, less reliable feline state"))
    rows = []
    for g, d in U.groupby("group"):
        for tag, med, lo, hi, p in [("uncorrected", "gap_median", "gap_lo", "gap_hi", "p_gap_le0"),
                                    ("Spearman-Brown", "gap_sb_median", "gap_sb_lo",
                                     "gap_sb_hi", "p_gap_sb_le0")]:
            rows.append({"group": g, "ceiling": tag, "n_pairs": len(d),
                         "gap_median": d[med].median(),
                         "gap_q25": d[med].quantile(.25), "gap_q75": d[med].quantile(.75),
                         "lower_bound_median": d[lo].median(),
                         "n_interval_excludes_zero": int((d[lo] > 0).sum()),
                         "max_tail_proportion": d[p].max()})
    G = pd.DataFrame(rows).sort_values(["group", "ceiling"])
    G.round(4).to_csv(OUT / "gap_by_group.tsv", sep="\t", index=False)
    print(G.round(3).to_string(index=False))

    # ---- 状態ラベル置換で「同種 vs 異種の gap 中央値差」を検定 ----
    states = sorted(set(U.a) | set(U.b))
    obs_cat = {s for s in states if s.startswith("cat")}
    n_cat = len(obs_cat)

    def stat(labels, col):
        cross = np.array([labels[r.a] != labels[r.b] for r in U.itertuples()])
        # 対照共有ペアは同種側からすでに除いてある
        a, b = U[col].to_numpy()[cross], U[col].to_numpy()[~cross]
        if len(a) < 5 or len(b) < 5:
            return np.nan
        return float(np.median(a) - np.median(b))

    res = []
    for col, tag in [("gap_median", "uncorrected"), ("gap_sb_median", "Spearman-Brown")]:
        obs = stat({s: ("cat" if s in obs_cat else "mouse") for s in states}, col)
        null = []
        for combo in combinations(range(len(states)), n_cat):
            lab = {s: "mouse" for s in states}
            for i in combo:
                lab[states[i]] = "cat"
            v = stat(lab, col)
            if np.isfinite(v):
                null.append(v)
        null = np.array(null)
        res.append({"ceiling": tag, "statistic": "median gap, cross minus same species",
                    "observed": obs, "n_perm": int(null.size),
                    "null_median": float(np.median(null)),
                    "null_p95": float(np.percentile(null, 95)),
                    "p_one_sided": float((null >= obs).sum() / null.size),
                    "p_two_sided": float((np.abs(null) >= abs(obs)).sum() / null.size)})
    R = pd.DataFrame(res)
    R.round(4).to_csv(OUT / "gap_group_permutation.tsv", sep="\t", index=False)
    print()
    print(R.round(4).to_string(index=False))
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
