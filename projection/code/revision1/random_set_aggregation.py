# -*- coding: utf-8 -*-
"""Part 2: サイズを合わせたランダム遺伝子集合で集約したときの centering の差。

本稿の主張は「pathway 平均という操作が centering の効果を増幅する」ことなので、
実際の 422 集合とサイズ分布だけを合わせたランダム集合で同じ差が出るかどうかが
直接の対照になる。4.13 の size-matched random set と同じ枠組み（300 draws）を、
経路単位ではなく集約段階に流用する。

  差が出る   → 効果は平均化という操作に由来する
  差が出ない → pathway の構造が効いているという積極的な所見

入力は Δ 行列・遺伝子リスト・ペア表・GMT だけ。既存の結果ファイルは読まない。
"""
from __future__ import annotations

import argparse
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=N_DRAW)
    a = ap.parse_args()

    missing = [fn for fn in PCS.COLL_FULL.values() if not (PCS.REF / fn).exists()]
    if missing:
        print(f"遺伝子セットの GMT が無い: {missing}", file=sys.stderr)
        return 2

    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    sizes = np.array([len(gi) for _c, _n, gi in paths])
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    n_genes = X.shape[0]

    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()

    def agg_auc(M, groups):
        A = np.vstack([M[g].mean(axis=0) for g in groups])
        An = A / np.linalg.norm(A, axis=0, keepdims=True)
        c = (An.T @ An)[ia, ib]
        return auc(c[m & wi], c[m & ~wi])

    real = [gi for _c, _n, gi in paths]
    obs_c, obs_r = agg_auc(C, real), agg_auc(X, real)
    print(f"実データ {len(paths)} 集合: centered {obs_c:.4f} / uncentered {obs_r:.4f} "
          f"/ 差 {obs_r - obs_c:+.4f}")

    rng = np.random.default_rng(SEED)
    rows = []
    for d in range(a.draws):
        groups = [rng.choice(n_genes, size=k, replace=False) for k in sizes]
        cc, rr = agg_auc(C, groups), agg_auc(X, groups)
        rows.append({"draw": d, "auc_agg_centred": cc, "auc_agg_uncentred": rr,
                     "difference": rr - cc})
    R = pd.DataFrame(rows)

    OUT.mkdir(parents=True, exist_ok=True)
    R.round(4).to_csv(OUT / "random_set_aggregation_draws.tsv", sep="\t", index=False)

    def q(col):
        v = R[col]
        return dict(median=v.median(), lo=v.quantile(.025), hi=v.quantile(.975),
                    mn=v.min(), mx=v.max())

    S = []
    for col in ("auc_agg_centred", "auc_agg_uncentred", "difference"):
        s = q(col)
        obs = {"auc_agg_centred": obs_c, "auc_agg_uncentred": obs_r,
               "difference": obs_r - obs_c}[col]
        S.append({"quantity": col, "observed_real_sets": obs, "random_median": s["median"],
                  "random_lo95": s["lo"], "random_hi95": s["hi"],
                  "random_min": s["mn"], "random_max": s["mx"],
                  "draws_at_or_above_observed": int((R[col] >= obs).sum()),
                  "n_draws": len(R)})
    T = pd.DataFrame(S)
    T.round(4).to_csv(OUT / "random_set_aggregation.tsv", sep="\t", index=False)

    print(f"\nサイズを合わせたランダム集合 {a.draws} 回:")
    for r in T.itertuples():
        print(f"  {r.quantity:20s} 実データ {r.observed_real_sets:.4f}   "
              f"ランダム 中央 {r.random_median:.4f} "
              f"[{r.random_lo95:.4f}, {r.random_hi95:.4f}] "
              f"（範囲 {r.random_min:.4f}–{r.random_max:.4f}）")
    d = R["difference"]
    print(f"\n  ランダム集合でも差は正: {int((d > 0).sum())}/{len(d)} draws")
    print(f"  実データの差 {obs_r - obs_c:.4f} 以上だった draw: "
          f"{int((d >= obs_r - obs_c).sum())}/{len(d)}")
    print(f"\n書き出し: {OUT / 'random_set_aggregation.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
