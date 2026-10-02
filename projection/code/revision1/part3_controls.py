# -*- coding: utf-8 -*-
"""Part 3（任意）の2件。

3-1 座標数の効果: Group A の 422 集合から 163 集合を無作為に抜いて集約 AUC を再計算する。
     matched Group B の 0.996 が「経路座標の本数が少ないから高い」のかを切り分ける。
3-2 IRI 12 か月状態: 唯一 age-matched control を使い 11 ペアに関与する状態を外したときの
     主結果の AUC。

入力は Δ 行列・遺伝子リスト・ペア表・GMT だけ。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402

OUT = HERE / "results" / "roundF"
SEED = 20260826
N_DRAW = 300


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def main() -> int:
    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    wi = P.wi.to_numpy()
    base = (~P.shared_control).to_numpy()

    def agg(M, groups, mask):
        A = np.vstack([M[g].mean(axis=0) for g in groups])
        An = A / np.linalg.norm(A, axis=0, keepdims=True)
        c = (An.T @ An)[ia, ib]
        return auc(c[mask & wi], c[mask & ~wi])

    def gene(M, mask):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        c = (Mn.T @ Mn)[ia, ib]
        return auc(c[mask & wi], c[mask & ~wi])

    real = [gi for _c, _n, gi in paths]
    rows = []

    # ---- 3-1 座標数の効果
    full_c, full_r = agg(C, real, base), agg(X, real, base)
    rng = np.random.default_rng(SEED)
    sub = []
    for _ in range(N_DRAW):
        pick = rng.choice(len(real), size=163, replace=False)
        g = [real[i] for i in pick]
        sub.append((agg(C, g, base), agg(X, g, base)))
    S = pd.DataFrame(sub, columns=["auc_agg_centred", "auc_agg_uncentred"])
    S["difference"] = S.auc_agg_uncentred - S.auc_agg_centred
    S.round(4).to_csv(OUT / "subset163_draws.tsv", sep="\t", index=False)
    for col, obs in (("auc_agg_centred", full_c), ("auc_agg_uncentred", full_r),
                     ("difference", full_r - full_c)):
        rows.append({"analysis": "163 of the 422 Group A sets", "quantity": col,
                     "observed_all_422": obs, "median": S[col].median(),
                     "lo95": S[col].quantile(.025), "hi95": S[col].quantile(.975),
                     "min": S[col].min(), "max": S[col].max(), "n": len(S)})
    print(f"3-1 Group A 422 集合: centered {full_c:.4f} / uncentered {full_r:.4f}")
    print(f"    163 集合に間引き（{N_DRAW} 回）: "
          f"centered 中央 {S.auc_agg_centred.median():.4f} "
          f"[{S.auc_agg_centred.quantile(.025):.4f}, {S.auc_agg_centred.quantile(.975):.4f}] / "
          f"uncentered 中央 {S.auc_agg_uncentred.median():.4f} "
          f"[{S.auc_agg_uncentred.quantile(.025):.4f}, {S.auc_agg_uncentred.quantile(.975):.4f}]")
    print(f"    matched Group B の 0.996 に届いた draw: "
          f"{int((S.auc_agg_uncentred >= 0.996).sum())}/{len(S)}")

    # ---- 3-2 IRI 12 か月状態を外す
    drop = (P.a == "IRI_12mo") | (P.b == "IRI_12mo")
    m2 = base & ~drop.to_numpy()
    print(f"\n3-2 IRI 12 か月が関与するペア: 全 {int(drop.sum())} 組、"
          f"対照非共有のうち {int((base & drop.to_numpy()).sum())} 組")
    for lab, mask in (("87 pairs (as reported)", base),
                      ("without the IRI 12-month state", m2)):
        g_, c_, r_ = gene(X, mask), agg(C, real, mask), agg(X, real, mask)
        n_in = int((mask & wi).sum()); n_cr = int((mask & ~wi).sum())
        rows.append({"analysis": lab, "quantity": "auc_gene_level", "observed_all_422": g_,
                     "median": np.nan, "lo95": np.nan, "hi95": np.nan,
                     "min": np.nan, "max": np.nan, "n": n_in + n_cr})
        rows.append({"analysis": lab, "quantity": "auc_agg_centred", "observed_all_422": c_,
                     "median": np.nan, "lo95": np.nan, "hi95": np.nan,
                     "min": np.nan, "max": np.nan, "n": n_in + n_cr})
        rows.append({"analysis": lab, "quantity": "auc_agg_uncentred", "observed_all_422": r_,
                     "median": np.nan, "lo95": np.nan, "hi95": np.nan,
                     "min": np.nan, "max": np.nan, "n": n_in + n_cr})
        rows.append({"analysis": lab, "quantity": "difference", "observed_all_422": r_ - c_,
                     "median": np.nan, "lo95": np.nan, "hi95": np.nan,
                     "min": np.nan, "max": np.nan, "n": n_in + n_cr})
        print(f"    {lab:32s} n={n_in + n_cr} ({n_in} within / {n_cr} cross)  "
              f"gene {g_:.4f}  centered {c_:.4f}  uncentered {r_:.4f}  差 {r_ - c_:+.4f}")

    T = pd.DataFrame(rows)
    T.round(4).to_csv(OUT / "part3_controls.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'part3_controls.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
