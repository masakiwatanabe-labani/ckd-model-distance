"""O: 経路粒度での cos θ。

Group A 空間（2,016 遺伝子）内で >=30 遺伝子を持つ経路に限り、
16状態 x 経路 の cos θ 行列を作る。基準軸は §2.2 と同じ 猫皮質 CKD3/4。

天井: 経路内の遺伝子だけを使った split-half から、セルごとに
      ceiling = sqrt(r_SB(基準軸) * r_SB(状態)) を出す。

制約への対応:
  - 遺伝子数が少ないと cos の分散が大きくなる。同じサイズのランダム遺伝子集合から
    cos の帰無分布を作り、経路の値がその幅の中か外かを必ず併記する。
  - 多重比較なので個別経路を「有意」とは呼ばない。分布とパターンのみ記述する。

Reactome はローカルの data/ref に無いため、GO BP / KEGG / Hallmark の3コレクションで実施。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from reliability import state_specs  # noqa: E402

OUT = HERE / "results" / "pathway"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260826
N_SPLIT = 200
MIN_GENES = 30
REF = "cat_CKD34"
COLLECTIONS = {"GO_BP": "geneset_gobp.gmt", "KEGG": "geneset_kegg.gmt",
               "Hallmark": "geneset_hallmark.gmt"}


def read_gmt(path):
    out = {}
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            out[f[0]] = {g.strip().upper() for g in f[2:] if g.strip()}
    return out


def cosv(A, Bm):
    """列ごとの cos。A,B は (genes, k)。"""
    an = A / np.linalg.norm(A, axis=0, keepdims=True)
    bn = Bm / np.linalg.norm(Bm, axis=0, keepdims=True)
    return (an * bn).sum(axis=0)


def main():
    rng = np.random.default_rng(SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes = D.index
    states = list(D.columns)
    X = D.to_numpy(float)
    gpos = {g: i for i, g in enumerate(genes)}
    print(f"Group A complete case: {len(genes)} 遺伝子 x {len(states)} 状態")

    # ---------------------------------------------- 経路
    paths = []
    for coll, fn in COLLECTIONS.items():
        gs = read_gmt(B.REF / fn)
        for name, members in gs.items():
            idx = np.array(sorted(gpos[g] for g in members & set(genes)))
            if len(idx) >= MIN_GENES:
                paths.append({"collection": coll, "pathway": name, "n_genes": len(idx),
                              "idx": idx})
        print(f"  {coll}: {len(gs)} sets → {sum(1 for p in paths if p['collection']==coll)} "
              f"sets with >= {MIN_GENES} Group A genes")
    print(f"合計 {len(paths)} 経路")

    # ---------------------------------------------- split-half Δ を状態ごとに1回だけ作る
    half = {}
    for state, mat, grp, label, ctrl, species in state_specs():
        if state not in states:
            continue
        A = mat[grp[grp == label].index]
        C = mat[grp[grp.isin(ctrl)].index]
        human = np.array([str(B.MAPS[species].get(g, g)).upper() for g in mat.index])
        seen, keep = set(), np.zeros(len(human), bool)
        for i, h in enumerate(human):
            if h not in seen:
                seen.add(h); keep[i] = True
        pos = pd.Series(np.arange(keep.sum()), index=human[keep]).reindex(genes)
        gi = pos.to_numpy()
        av, cv = A.to_numpy(float), C.to_numpy(float)
        na, nc = av.shape[1], cv.shape[1]
        d1 = np.empty((N_SPLIT, len(genes)), np.float32)
        d2 = np.empty((N_SPLIT, len(genes)), np.float32)
        for b in range(N_SPLIT):
            pa, pc = rng.permutation(na), rng.permutation(nc)
            h1, hk = na // 2, nc // 2
            e1 = (np.nanmean(av[:, pa[:h1]], 1) - np.nanmean(cv[:, pc[:hk]], 1))[keep]
            e2 = (np.nanmean(av[:, pa[h1:2*h1]], 1) - np.nanmean(cv[:, pc[hk:2*hk]], 1))[keep]
            d1[b] = e1[gi]; d2[b] = e2[gi]
        half[state] = (d1, d2)
    print("split-half Δ の構築完了")

    def reliability(state, idx):
        d1, d2 = half[state]
        a, b = d1[:, idx].T, d2[:, idx].T          # (genes, splits)
        r = np.median(cosv(a, b))
        return float(np.clip(2 * r / (1 + r), 0, 1)) if r > 0 else 0.0

    # ---------------------------------------------- ランダム集合の帰無分布（サイズ別）
    ref_i = states.index(REF)
    sizes = sorted({30, 50, 100, 200, 400, 800})
    null_rows = []
    for n in sizes:
        if n > len(genes):
            continue
        draws = np.array([rng.choice(len(genes), size=n, replace=False) for _ in range(300)])
        for s in states:
            si = states.index(s)
            v = np.array([cosv(X[d][:, [ref_i]], X[d][:, [si]])[0] for d in draws])
            null_rows.append({"n_genes": n, "state": s, "lo95": float(np.percentile(v, 2.5)),
                              "hi95": float(np.percentile(v, 97.5)),
                              "iqr": float(np.percentile(v, 75) - np.percentile(v, 25))})
    NULL = pd.DataFrame(null_rows)
    NULL.round(4).to_csv(OUT / "random_set_null.tsv", sep="\t", index=False)

    def null_band(n, state):
        d = NULL[NULL.state == state]
        i = (d.n_genes - n).abs().idxmin()
        return float(NULL.loc[i, "lo95"]), float(NULL.loc[i, "hi95"])

    # ---------------------------------------------- 本体
    rows = []
    for p in paths:
        idx = p["idx"]
        Xi = X[idx]
        c = cosv(Xi[:, [ref_i]] * np.ones((1, len(states))), Xi)
        r_ref = reliability(REF, idx)
        for si, s in enumerate(states):
            lo, hi = null_band(len(idx), s)
            ceil = float(np.sqrt(max(r_ref, 0) * reliability(s, idx)))
            rows.append({"collection": p["collection"], "pathway": p["pathway"],
                         "n_genes": p["n_genes"], "state": s, "cos": float(c[si]),
                         "ceiling": ceil, "gap_to_ceiling": ceil - float(c[si]),
                         "random_lo95": lo, "random_hi95": hi,
                         "above_random": bool(c[si] > hi)})
    R = pd.DataFrame(rows)
    R.round(4).to_csv(OUT / "pathway_cos.tsv", sep="\t", index=False)
    print(f"セル数: {len(R)}")

    piv = R.pivot_table(index=["collection", "pathway", "n_genes"],
                        columns="state", values="cos")
    piv.round(4).to_csv(OUT / "pathway_cos_matrix.tsv", sep="\t")
    print(f"書き出し: {OUT}")


if __name__ == "__main__":
    main()
