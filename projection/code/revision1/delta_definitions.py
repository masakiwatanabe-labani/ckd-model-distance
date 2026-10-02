# -*- coding: utf-8 -*-
"""AC: Δ の定義を変えて主解析（16状態・中心化の有無）を通す。

査読者 3 major 3 への対応。主解析は log2FC のまま固定し、他は感度解析として並置する。

  1. log2FC                現行。群平均 log2 発現の差
  2. per-gene z            各遺伝子を 16 状態にわたって標準化
  3. rank                  各状態の中で Δ を順位化し、[0,1] に写す
  4. quantile              Δ を状態間で分位正規化

各定義について次を出す。
  - 16 状態の state-wide mean Δ と、猫 4 状態とマウス 12 状態が分離するか
  - 遺伝子レベル AUC、中心化した遺伝子レベル AUC
  - 経路集約 AUC（中心化あり・なし）

注意: 順位化と分位正規化は、状態ごとの周辺分布を同じにする操作なので、
state-wide mean は構成上すべての状態で同じ値になる。分離が消えるのは当然で、
「定義依存かどうか」を見るには centering の対比のほうが情報を持つ。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402

OUT = HERE / "results" / "roundR"


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def as_log2fc(X):
    return X.copy()


def as_zscore(X):
    """遺伝子ごとに 16 状態にわたって標準化する。"""
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, ddof=1, keepdims=True)
    sd = np.where(sd > 0, sd, np.nan)
    return (X - mu) / sd


def as_rank(X):
    """状態ごとに遺伝子を順位化し、(0,1) に写す。"""
    out = np.empty_like(X)
    n = X.shape[0]
    for j in range(X.shape[1]):
        order = np.argsort(np.argsort(X[:, j], kind="mergesort"), kind="mergesort")
        out[:, j] = (order + 0.5) / n
    return out


def as_quantile(X):
    """Δ を状態間で分位正規化する（各状態の周辺分布を共通の平均分位に合わせる）。"""
    S = np.sort(X, axis=0)
    ref = S.mean(axis=1)
    out = np.empty_like(X)
    n = X.shape[0]
    for j in range(X.shape[1]):
        order = np.argsort(np.argsort(X[:, j], kind="mergesort"), kind="mergesort")
        out[:, j] = ref[order]
    assert n == X.shape[0]
    return out


DEFS = [("log2FC (main analysis)", as_log2fc),
        ("per-gene standardization", as_zscore),
        ("rank within state", as_rank),
        ("quantile normalization", as_quantile)]


def main() -> int:
    missing = [fn for fn in PCS.COLL_FULL.values() if not (PCS.REF / fn).exists()]
    if missing:
        print(f"GMT が無い: {missing}", file=sys.stderr)
        return 2

    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    states = list(D.columns)
    is_cat = np.array([s.startswith("cat") for s in states])
    idx = {s: i for i, s in enumerate(states)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()
    X0 = D.to_numpy(float)

    def cosmat(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        return (Mn.T @ Mn)[ia, ib]

    rows, means = [], []
    for label, fn in DEFS:
        X = fn(X0)
        ok = np.all(np.isfinite(X), axis=1)
        X = X[ok]
        C = X - X.mean(axis=0)
        gi_ok = {}
        keep = np.flatnonzero(ok)
        remap = {g: i for i, g in enumerate(keep)}
        for c, n_, gi in paths:
            gi2 = np.array([remap[g] for g in gi if g in remap])
            if len(gi2) >= 30:
                gi_ok[(c, n_)] = gi2
        agg_c = np.vstack([C[g].mean(axis=0) for g in gi_ok.values()])
        agg_r = np.vstack([X[g].mean(axis=0) for g in gi_ok.values()])

        mu = X.mean(axis=0)
        cat, mo = mu[is_cat], mu[~is_cat]
        sep = bool(cat.max() < mo.min() or mo.max() < cat.min())
        for s, v in zip(states, mu):
            means.append({"definition": label, "state": s,
                          "species": "feline" if s.startswith("cat") else "mouse",
                          "mean_delta": float(v)})

        def A(M):
            c = cosmat(M)
            return auc(c[m & wi], c[m & ~wi])

        rows.append({"definition": label, "n_genes": int(X.shape[0]),
                     "n_sets": len(gi_ok),
                     "cat_mean_min": float(cat.min()), "cat_mean_max": float(cat.max()),
                     "mouse_mean_min": float(mo.min()), "mouse_mean_max": float(mo.max()),
                     "mean_separates_species": sep,
                     "mean_gap": float(mo.min() - cat.max()),
                     "auc_gene_level": A(X), "auc_gene_centered": A(C),
                     "auc_agg_uncentered": A(agg_r), "auc_agg_centered": A(agg_c),
                     "centering_drop": A(agg_r) - A(agg_c)})

    # 条件ごとのコサイン分布も残す（本文が「0.99 付近に潰れる」と書く根拠）
    dist = []
    for label, fn in DEFS:
        X = fn(X0)
        X = X[np.all(np.isfinite(X), axis=1)]
        C = X - X.mean(axis=0)
        keep = np.flatnonzero(np.all(np.isfinite(fn(X0)), axis=1))
        remap = {g: i for i, g in enumerate(keep)}
        gg = [np.array([remap[g] for g in gi if g in remap]) for _c, _n, gi in paths]
        gg = [g for g in gg if len(g) >= 30]
        for cond, Mx in (("gene level", X),
                         ("aggregate, uncentered", np.vstack([X[g].mean(axis=0) for g in gg])),
                         ("aggregate, centered", np.vstack([C[g].mean(axis=0) for g in gg]))):
            c = cosmat(Mx)
            dist.append({"definition": label, "condition": cond,
                         "within_median": float(np.median(c[m & wi])),
                         "cross_median": float(np.median(c[m & ~wi])),
                         "min": float(c[m].min()), "max": float(c[m].max()),
                         "spread": float(c[m].max() - c[m].min())})
    pd.DataFrame(dist).round(4).to_csv(OUT / "delta_definitions_cosines.tsv", sep="\t",
                                       index=False)

    T, MU = pd.DataFrame(rows), pd.DataFrame(means)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "delta_definitions.tsv", sep="\t", index=False)
    MU.round(4).to_csv(OUT / "delta_definitions_state_means.tsv", sep="\t", index=False)

    print(f"{'definition':28s} {'genes':>6s} {'sets':>5s} {'gene':>7s} {'gene-c':>7s} "
          f"{'agg-unc':>8s} {'agg-cen':>8s} {'drop':>7s}  mean separates")
    for r in T.itertuples():
        print(f"{r.definition:28s} {r.n_genes:6d} {r.n_sets:5d} {r.auc_gene_level:7.4f} "
              f"{r.auc_gene_centered:7.4f} {r.auc_agg_uncentered:8.4f} "
              f"{r.auc_agg_centered:8.4f} {r.centering_drop:+7.4f}  "
              f"{'はい' if r.mean_separates_species else 'いいえ'}"
              f"（猫 {r.cat_mean_min:+.3f}〜{r.cat_mean_max:+.3f} / "
              f"マウス {r.mouse_mean_min:+.3f}〜{r.mouse_mean_max:+.3f}）")
    print(f"\n書き出し: {OUT / 'delta_definitions.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
