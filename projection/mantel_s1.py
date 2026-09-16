"""L: 旧稿の Mantel 解析を S1 空間・新16状態で再計算する。

条件を新稿に統一する:
    遺伝子空間 Group A intersection（cortex ∩ medulla proteome, ortholog）
    Δ 定義     log2FC, 群平均 FPKM>=1 フィルタあり（S1）
    状態       新稿の16状態
旧稿の公表値（Mantel rho 0.553、入口 −0.019 等）とは状態集合も Δ 定義も違うので、
食い違って構わない。食い違いは NUMBERS.md に記録する。

統合の核（task 4）:
    旧稿の距離は D = 1 − Spearman rho。**Spearman は単調変換に不変なので、
    この距離は既に「方向のみ」の量であり、応答振幅を見ていない。**
    そこで3つの非類似度を並べ、時間との相関が方向由来か振幅由来かを分ける。

        D_rho = 1 − Spearman(Δi, Δj)        順位ベースの方向（旧稿の定義）
        D_cos = 1 − cos θ(Δi, Δj)           線形の方向
        D_nrm = |log2(‖Δi‖ / ‖Δj‖)|         振幅のみ

    D_rho と D_cos が時間と相関し D_nrm が相関しなければ、
    旧稿の Mantel は「時間的に近い状態は向きが似る」を測っていたことになる。
    これは §2.2 の「各状態の cos は基準に対してほぼ平坦」と矛盾しない
    （固定の基準軸に対する角度が一定でも、軌道は連続的に回転しうる）。
"""
from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

# 遺伝子空間の差し替え（AA. 感度解析）。既定は Group A intersection で、
# 環境変数を与えない限り従来と完全に同じ動作になる。
TAG = os.environ.get("XSP_TAG", "")
GENES = Path(os.environ.get("XSP_GENES", str(HERE / "groupA_intersection.txt")))
OUT = HERE / "results" / ("mantel_s1" + TAG)
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260826
N_PERM = 9999
ENTRY = {"cat_CKD12": "tubular", "cat_CKD34": "tubular",
         "cat_med_CKD12": "tubular", "cat_med_CKD34": "tubular",
         "mouse_5D": "glomerular", "mouse_2W": "glomerular", "mouse_3W": "glomerular"}
ENTRY.update({f"IRI_{t}": "tubular" for t in
              ["2h", "4h", "24h", "48h", "72h", "7d", "14d", "28d", "12mo"]})
MODEL = {k: ("cat_natural" if k.startswith("cat") else
             "PodTRECK" if k.startswith("mouse") else "IRI") for k in ENTRY}
SPECIES = {k: ("cat" if k.startswith("cat") else "mouse") for k in ENTRY}


def dissimilarities(D: pd.DataFrame):
    """3つの非類似度行列。"""
    states = list(D.columns)
    X = D.to_numpy(float)
    n = len(states)
    R = np.zeros((n, n)); C = np.zeros((n, n)); N = np.zeros((n, n))
    nrm = np.linalg.norm(X, axis=0)
    for i, j in combinations(range(n), 2):
        r = 1 - stats.spearmanr(X[:, i], X[:, j])[0]
        c = 1 - float(X[:, i] @ X[:, j] / (nrm[i] * nrm[j]))
        a = abs(np.log2(nrm[i] / nrm[j]))
        R[i, j] = R[j, i] = r
        C[i, j] = C[j, i] = c
        N[i, j] = N[j, i] = a
    return states, {"D_rho (1-Spearman)": R, "D_cos (1-cos)": C, "D_nrm (|log2 norm ratio|)": N}


def mantel(D, X, rng, partial=None, n_perm=N_PERM):
    """Mantel / 偏Mantel。状態ラベルの並べ替えで p を取る（14_time_axis_full.py と同じ）。"""
    k = D.shape[0]
    iu = np.triu_indices(k, 1)

    def resid(y, z):
        A = np.column_stack([np.ones_like(z), z])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        return y - A @ beta

    d, x = D[iu], X[iu]
    if partial is not None:
        z = partial[iu]
        d, x = resid(d, z), resid(x, z)
    r_obs = stats.spearmanr(d, x)[0]
    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(k)
        dp = D[np.ix_(p, p)][iu]
        if partial is not None:
            dp = resid(dp, partial[iu])
        if abs(stats.spearmanr(dp, x)[0]) >= abs(r_obs) - 1e-12:
            cnt += 1
    return float(r_obs), (cnt + 1) / (n_perm + 1)


def run(D, states_use, label, hours, rng):
    sub = D[states_use]
    states, mats = dissimilarities(sub)
    n = len(states)
    T = np.zeros((n, n)); E = np.zeros((n, n)); M = np.zeros((n, n))
    have_t = all(np.isfinite(hours.get(s, np.nan)) for s in states)
    for i, j in combinations(range(n), 2):
        a, b = states[i], states[j]
        if have_t:
            T[i, j] = T[j, i] = abs(np.log10(hours[a]) - np.log10(hours[b]))
        E[i, j] = E[j, i] = float(ENTRY[a] != ENTRY[b])
        M[i, j] = M[j, i] = float(MODEL[a] != MODEL[b])
    iu = np.triu_indices(n, 1)
    conf = stats.spearmanr(T[iu], E[iu])[0] if have_t else np.nan

    rows = []
    for mname, DM in mats.items():
        preds = [("onset compartment", E), ("model identity", M)]
        if have_t:
            preds.append(("elapsed time |Δlog10 h|", T))
        for pname, P in preds:
            r, p = mantel(DM, P, rng)
            rows.append({"subset": label, "n_states": n, "dissimilarity": mname,
                         "predictor": pname, "mantel_rho": r, "p_perm": p})
        if have_t:
            r, p = mantel(DM, T, rng, partial=E)
            rows.append({"subset": label, "n_states": n, "dissimilarity": mname,
                         "predictor": "time | onset compartment", "mantel_rho": r, "p_perm": p})
            r, p = mantel(DM, E, rng, partial=T)
            rows.append({"subset": label, "n_states": n, "dissimilarity": mname,
                         "predictor": "onset compartment | time", "mantel_rho": r, "p_perm": p})
    return pd.DataFrame(rows), conf, states, mats


def main():
    rng = np.random.default_rng(SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    print(f"遺伝子 {len(D)}（Group A intersection の complete case。除外なし）")

    tm = pd.read_csv(HERE / "time_map.tsv", sep="\t")
    hours = tm.set_index("state")["hours"].astype(float).to_dict()
    model_states = [s for s in D.columns if np.isfinite(hours.get(s, np.nan))]

    res = []
    r16, conf16, states16, mats16 = run(D, list(D.columns), "16 states (all)", hours, rng)
    res.append(r16)
    r12, conf12, _, _ = run(D, model_states, "12 model states (time defined)", hours, rng)
    res.append(r12)
    R = pd.concat(res, ignore_index=True)
    R.round(4).to_csv(OUT / "mantel_s1.tsv", sep="\t", index=False)

    print(f"\n交絡 rho(時間, 入口) = {conf12:+.3f}（12状態）")
    for lab in ["16 states (all)", "12 model states (time defined)"]:
        d = R[R.subset == lab]
        if not len(d):
            continue
        print(f"\n=== {lab} ===")
        piv = d.pivot(index="predictor", columns="dissimilarity", values="mantel_rho")
        pv = d.pivot(index="predictor", columns="dissimilarity", values="p_perm")
        for pred in piv.index:
            cells = "  ".join(f"{c.split(' ')[0]:6s} {piv.loc[pred,c]:+.3f} (p={pv.loc[pred,c]:.4f})"
                              for c in piv.columns)
            print(f"  {pred:28s} {cells}")

    # ---------------- task 3: 種内 vs 種間の距離分布
    rows = []
    for i, j in combinations(range(len(states16)), 2):
        a, b = states16[i], states16[j]
        cls = "within species" if SPECIES[a] == SPECIES[b] else "cross species"
        rec = {"a": a, "b": b, "class": cls,
               "same_dataset": MODEL[a] == MODEL[b]}
        for mname, DM in mats16.items():
            rec[mname] = DM[i, j]
        rows.append(rec)
    P = pd.DataFrame(rows)
    P.round(4).to_csv(OUT / "distance_pairs.tsv", sep="\t", index=False)
    summ = []
    for mname in mats16:
        for cls in ["within species", "cross species"]:
            v = P.loc[P["class"] == cls, mname]
            summ.append({"dissimilarity": mname, "class": cls, "n": len(v),
                         "median": v.median(), "min": v.min(), "max": v.max()})
    S = pd.DataFrame(summ)
    S.round(4).to_csv(OUT / "distance_distributions.tsv", sep="\t", index=False)
    print("\n=== 種内 vs 種間の距離分布（16状態、120ペア） ===")
    print(S.round(3).to_string(index=False))
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
