"""Q追加: 蛋白層の split-half r を分割サイズの関数として確認する。

猫の信頼性 0.51–0.81 が検体数由来の過小評価かどうかを、
RNA と同じ手続き（r(n) 曲線 + Spearman-Brown の内部検証）で確かめる。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from protein_decomp import SPECS, load, cosine, DET  # noqa: E402

OUT = HERE / "results" / "protein"
N_SPLIT = 500
SEED = 20260826


def sb(r, k=2):
    r = np.clip(r, 0, 0.999999)
    return k * r / (1 + (k - 1) * r)


def main():
    rng = np.random.default_rng(SEED)
    D = pd.read_csv(OUT / "delta_matrix_protein.tsv", sep="\t", index_col=0)
    genes = D.index
    rows = []
    for state, mname, case, ctrl, species in SPECS:
        mat, grp = load(mname)
        cs = list(grp[grp == case].index)
        ks = list(grp[grp == ctrl].index)
        nmax = min(len(cs), len(ks)) // 2
        for n in range(1, nmax + 1):
            vals = []
            for _ in range(N_SPLIT):
                c, k = rng.permutation(cs), rng.permutation(ks)
                d1 = mat[list(c[:n])].mean(axis=1) - mat[list(k[:n])].mean(axis=1)
                d2 = mat[list(c[n:2*n])].mean(axis=1) - mat[list(k[n:2*n])].mean(axis=1)
                h1, _ = B.to_human(d1[np.isfinite(d1)], species)
                h2, _ = B.to_human(d2[np.isfinite(d2)], species)
                u, v = h1.reindex(genes), h2.reindex(genes)
                m = u.notna() & v.notna()
                if m.sum() > 50:
                    vals.append(cosine(u[m].to_numpy(), v[m].to_numpy()))
            rows.append({"state": state, "n_case": len(cs), "n_ctrl": len(ks),
                         "n_per_half": n, "n_max": nmax, "r_half": float(np.median(vals))})
            print(f"{state:20s} n={n} (max {nmax})  r_half={np.median(vals):+.3f}")
    C = pd.DataFrame(rows)
    C.round(4).to_csv(OUT / "protein_r_curve.tsv", sep="\t", index=False)

    print("\n=== Spearman-Brown の内部検証（曲線が取れる状態） ===")
    val = []
    for state, d in C.groupby("state"):
        d = d.set_index("n_per_half")["r_half"]
        if 1 not in d.index or len(d) < 2:
            continue
        for k in [2, 3]:
            if k in d.index:
                pred = float(sb(d[1], k))
                val.append({"state": state, "k": k, "r_1": d[1], "r_k_obs": d[k],
                            "r_k_pred": pred, "obs_minus_pred": d[k] - pred})
                print(f"  {state:20s} k={k}  r(1)={d[1]:+.3f}  実測 r(k)={d[k]:+.3f}  "
                      f"SB予言={pred:.3f}  差={d[k]-pred:+.3f}")
    V = pd.DataFrame(val)
    if len(V):
        V.round(4).to_csv(OUT / "protein_sb_validation.tsv", sep="\t", index=False)
        print(f"\n  実測−予言: 範囲 [{V.obs_minus_pred.min():+.3f}, {V.obs_minus_pred.max():+.3f}], "
              f"中央 {V.obs_minus_pred.median():+.3f}, 実測>予言 {int((V.obs_minus_pred>0).sum())}/{len(V)}")

    # 全検体まで外挿した信頼性（n_max の実測から）
    print("\n=== 全検体版への外挿 ===")
    ext = []
    for state, d in C.groupby("state"):
        d = d.sort_values("n_per_half")
        r = d.r_half.iloc[-1]; nmax = int(d.n_max.iloc[0])
        kfull = min(d.n_case.iloc[0], d.n_ctrl.iloc[0]) / nmax
        ext.append({"state": state, "n_max": nmax, "r_at_nmax": r,
                    "k_to_full": kfull, "reliability_extrapolated": float(sb(r, kfull))})
        print(f"  {state:20s} r({nmax})={r:+.3f} → 全検体 {float(sb(r,kfull)):.3f}")
    pd.DataFrame(ext).round(4).to_csv(OUT / "protein_reliability_extrapolated.tsv",
                                      sep="\t", index=False)


if __name__ == "__main__":
    main()
