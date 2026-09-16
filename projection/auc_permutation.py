"""E: 種差 AUC の検定を状態レベル並べ替えに変える。

ペア単位の並べ替えや Mann-Whitney の漸近 p は、120 ペアを独立扱いしてしまう。
実際には 16 状態から作った依存ペアなので、実質の自由度は状態数に近い。
状態に貼られた species ラベルのほうを並べ替えて AUC の帰無分布を作る。

3つの水準で出す:
  1. 状態レベル網羅並べ替え  C(16,4) = 1820 通りを全列挙（厳密 p）
  2. 状態レベル・データセット構成保存  各データセットから抜く状態数を観測と揃える
  3. データセットレベル       どのデータセットを「種 A」とするかの 3 通り（最も厳しい）

3 は species が dataset と完全に交絡している（ネコのデータセットは1つしかない）ため
p >= 1/3 にしかならない。それが本質的な限界であることを明示するために出す。
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def auc(x, y):
    """P(X > Y) + 0.5 P(X = Y)。"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if x.size == 0 or y.size == 0:
        return np.nan
    gt = (x[:, None] > y[None, :]).sum()
    eq = (x[:, None] == y[None, :]).sum()
    return float((gt + 0.5 * eq) / (x.size * y.size))


def auc_for_labels(pairs, species_of, col):
    """与えられた state->species 割り当てで、異データセットペアの AUC。"""
    sa = pairs["a"].map(species_of)
    sb = pairs["b"].map(species_of)
    same = pairs[col][(sa == sb).to_numpy()]
    cross = pairs[col][(sa != sb).to_numpy()]
    if len(same) < 3 or len(cross) < 3:
        return np.nan
    return auc(same, cross)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="results/control_analysis/pairs.tsv")
    ap.add_argument("--col", default="cos")
    ap.add_argument("--out", default="results/control_analysis")
    a = ap.parse_args()

    P = pd.read_csv(HERE / a.pairs, sep="\t")
    P = P[P.dataset_a != P.dataset_b].reset_index(drop=True)   # 異データセットペアのみ
    states = sorted(set(P.a) | set(P.b))
    ds = {}
    for _, r in P.iterrows():
        ds[r["a"]] = r["dataset_a"]
        ds[r["b"]] = r["dataset_b"]
    obs_species = {s: ("cat" if ds[s] == "cat" else "mouse") for s in states}
    n_cat = sum(v == "cat" for v in obs_species.values())

    obs = auc_for_labels(P, obs_species, a.col)

    # ---- 1. 状態レベル網羅並べ替え -----------------------------------
    null1 = []
    for combo in combinations(range(len(states)), n_cat):
        lab = {s: "mouse" for s in states}
        for i in combo:
            lab[states[i]] = "cat"
        v = auc_for_labels(P, lab, a.col)
        if np.isfinite(v):
            null1.append(v)
    null1 = np.array(null1)
    p1 = float((null1 >= obs).sum() / null1.size)

    # ---- 2. データセット構成を保存した並べ替え ------------------------
    # 各データセットから「cat」に割り当てる状態数を観測と同じにする
    by_ds = {}
    for s in states:
        by_ds.setdefault(ds[s], []).append(s)
    want = {d: sum(1 for s in v if obs_species[s] == "cat") for d, v in by_ds.items()}
    opts = {d: list(combinations(range(len(v)), want[d])) for d, v in by_ds.items()}
    null2 = []
    import itertools
    for pick in itertools.product(*[opts[d] for d in by_ds]):
        lab = {s: "mouse" for s in states}
        for d, idxs in zip(by_ds, pick):
            for i in idxs:
                lab[by_ds[d][i]] = "cat"
        v = auc_for_labels(P, lab, a.col)
        if np.isfinite(v):
            null2.append(v)
    null2 = np.array(null2)
    p2 = float((null2 >= obs).sum() / null2.size) if null2.size else np.nan

    # ---- 3. データセットレベル ---------------------------------------
    null3 = []
    for d in by_ds:
        lab = {s: ("cat" if ds[s] == d else "mouse") for s in states}
        v = auc_for_labels(P, lab, a.col)
        null3.append(v if np.isfinite(v) else np.nan)
    null3 = np.array([v for v in null3 if np.isfinite(v)])
    p3 = float((null3 >= obs).sum() / null3.size) if null3.size else np.nan

    res = pd.DataFrame([
        {"level": "1. 状態レベル網羅 C(16,4)", "n_perm": null1.size, "auc_obs": obs,
         "null_median": float(np.median(null1)), "null_p95": float(np.percentile(null1, 95)),
         "null_max": float(null1.max()), "p": p1, "note": ""},
        {"level": "2. データセット構成保存", "n_perm": null2.size, "auc_obs": obs,
         "null_median": float(np.median(null2)) if null2.size else np.nan,
         "null_p95": float(np.percentile(null2, 95)) if null2.size else np.nan,
         "null_max": float(null2.max()) if null2.size else np.nan,
         "p": p2 if null2.size > 1 else np.nan,
         "note": ("退化（並べ替えが観測の1通りしかない）" if null2.size <= 1 else "")},
        {"level": "3. データセットレベル（最厳）", "n_perm": null3.size, "auc_obs": obs,
         "null_median": float(np.median(null3)), "null_p95": float(np.max(null3)),
         "null_max": float(null3.max()), "p": p3,
         "note": "species が dataset と完全交絡のため p >= 1/3 が下限"},
    ])
    out = HERE / a.out
    out.mkdir(parents=True, exist_ok=True)
    tag = "" if a.col == "cos" else f"_{a.col}"
    res.round(4).to_csv(out / f"auc_permutation{tag}.tsv", sep="\t", index=False)
    pd.DataFrame({"auc_null_state_level": null1}).to_csv(
        out / f"auc_null_state_level{tag}.tsv", sep="\t", index=False)

    print(f"指標: {a.col}   異データセットペア {len(P)} 組 / 状態 {len(states)}")
    print(res.round(4).to_string(index=False))
    print(f"\n観測 AUC = {obs:.4f}")
    print(f"状態レベル網羅並べ替えの帰無分布: 中央値 {np.median(null1):.3f}, "
          f"95パーセンタイル {np.percentile(null1, 95):.3f}, 最大 {null1.max():.3f}")
    print("水準3 は species が dataset と完全交絡（ネコのデータセットは1つ）なので "
          "p >= 1/3 が下限。これは検定の限界であって効果の不在ではない。")


if __name__ == "__main__":
    main()
