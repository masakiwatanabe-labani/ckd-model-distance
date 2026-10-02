# -*- coding: utf-8 -*-
"""AE2: ランダム集合の設計を 2 つ並べる（査読者 3 major 4）。

指摘は「サイズだけを合わせた乱択は、422 経路のあいだの遺伝子の重複・再利用を
保存していない」。そこで、重複構造を完全に保存する対照を追加する。

  design A: size-matched draws      各集合のサイズだけ合わせて独立に引く（現行）
  design B: gene-label permutation  422 集合の所属をそのまま保ち、Δ 行列の遺伝子
                                    ラベルだけを並べ替える。集合サイズ・全ペアの
                                    Jaccard・遺伝子の再利用回数が厳密に保存され、
                                    壊れるのは「どの遺伝子がどの経路に属するか」だけ。

入力は Δ 行列・ペア表・GMT のみ。既存の結果ファイルは読まない。
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

OUT = HERE / "results" / "roundR"
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
        print(f"GMT が無い: {missing}", file=sys.stderr)
        return 2

    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    real = [gi for _c, _n, gi in paths]
    sizes = np.array([len(gi) for gi in real])
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

    obs_c, obs_u = agg_auc(C, real), agg_auc(X, real)
    print(f"実データ {len(real)} 集合 / {int(m.sum())} ペア: "
          f"非中心化 {obs_u:.4f} / 中心化 {obs_c:.4f} / 差 {obs_u - obs_c:+.4f}")

    # 重複構造の記述（保存されていることを示す量）
    mem = np.zeros(n_genes, int)
    for gi in real:
        mem[gi] += 1
    print(f"  遺伝子 1 個が属する集合数: 中央 {int(np.median(mem))} / 最大 {mem.max()} / "
          f"どの集合にも属さない遺伝子 {int((mem == 0).sum())}")
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"quantity": "sets per gene, median", "value": float(np.median(mem))},
                  {"quantity": "sets per gene, maximum", "value": float(mem.max())},
                  {"quantity": "genes in no set", "value": float((mem == 0).sum())},
                  {"quantity": "genes in the universe", "value": float(n_genes)},
                  {"quantity": "pathway sets", "value": float(len(real))}]).to_csv(
        OUT / "random_set_designs_overlap.tsv", sep="\t", index=False)

    designs = {}
    for name in ("size-matched draws", "gene-label permutation"):
        rng = np.random.default_rng(SEED)
        rows = []
        for d in range(a.draws):
            if name == "size-matched draws":
                groups = [rng.choice(n_genes, size=k, replace=False) for k in sizes]
            else:
                perm = rng.permutation(n_genes)
                groups = [perm[gi] for gi in real]
            cc, uu = agg_auc(C, groups), agg_auc(X, groups)
            rows.append({"draw": d, "auc_agg_centred": cc, "auc_agg_uncentred": uu,
                         "difference": uu - cc})
        designs[name] = pd.DataFrame(rows)

    rec = []
    for name, R in designs.items():
        for col, obs in (("auc_agg_uncentred", obs_u), ("auc_agg_centred", obs_c),
                         ("difference", obs_u - obs_c)):
            v = R[col]
            rec.append({"design": name, "quantity": col, "n_draws": len(R),
                        "observed_real_sets": obs, "random_median": v.median(),
                        "random_min": v.min(), "random_max": v.max(),
                        "draws_at_or_above_observed": int((v >= obs).sum())})
    T = pd.DataFrame(rec)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "random_set_designs.tsv", sep="\t", index=False)
    for name, R in designs.items():
        tag = "sizematched" if name.startswith("size") else "labelperm"
        R.round(4).to_csv(OUT / f"random_set_designs_draws_{tag}.tsv", sep="\t", index=False)

    for name, R in designs.items():
        print(f"\n{name}（{a.draws} 回）")
        for col in ("auc_agg_uncentred", "auc_agg_centred", "difference"):
            v = R[col]
            print(f"  {col:20s} 中央 {v.median():.4f}（範囲 {v.min():.4f}–{v.max():.4f}）")
        d = R["difference"]
        print(f"  差が正の draw: {int((d > 0).sum())}/{len(d)}   "
              f"実データの差 {obs_u - obs_c:.4f} 以上: {int((d >= obs_u - obs_c).sum())}/{len(d)}")
    print(f"\n書き出し: {OUT / 'random_set_designs.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
