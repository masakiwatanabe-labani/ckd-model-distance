# -*- coding: utf-8 -*-
"""AB: 重複ヒトシンボルの「先頭を採用」ルールが入力行順に依存する件（査読者 3 minor）。

to_human() は `~index.duplicated(keep="first")` で重複を落とすので、
同じヒトシンボルに複数の元遺伝子が当たるとき、どれが残るかは入力の行順で決まる。
行順をシャッフルして主要結果が動くかを調べる。

やり方: 元空間の Δ（lfc）は行順に依存しないので、状態ごとに一度だけ計算し、
各シャッフルでは「重複シンボルの代表をどれにするか」だけを引き直す。
これは入力行をシャッフルして to_human を通すのと同じ結果になり、はるかに速い。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "code" / "revision1"))
import build_delta_matrix as B  # noqa: E402
import pathway_centering_steps as PCS  # noqa: E402

OUT = HERE / "results" / "roundR"
N_SHUFFLE = 100
SEED = 20260826


def source_deltas():
    """状態ごとの Δ を元空間（猫・マウスのシンボル）のまま返す。"""
    out = {}
    for state, s in B.cat_deltas().items():
        out[state] = ("cat", s[np.isfinite(s)])
    for state, s in B.podtreck_deltas().items():
        out[state] = ("mouse", s[np.isfinite(s)])
    for state, s in B.iri_deltas().items():
        out[state] = ("mouse", s[np.isfinite(s)])
    return out


def human_of(idx, species):
    mp = B.MAPS[species]
    return np.array([str(mp[g]).upper() if g in mp else str(g).upper() for g in idx])


def build(sd, rng=None):
    """代表の選び方を変えて Δ 行列を組む。rng が None なら現行（先頭を採用）。"""
    cols = {}
    for state, (sp, s) in sd.items():
        h = human_of(s.index, sp)
        df = pd.DataFrame({"h": h, "v": s.to_numpy()})
        if rng is None:
            pick = df.drop_duplicates("h", keep="first")
        else:
            order = rng.permutation(len(df))
            pick = df.iloc[order].drop_duplicates("h", keep="first")
        cols[state] = pd.Series(pick.v.to_numpy(), index=pick.h.to_numpy())
    D = pd.DataFrame(cols).sort_index()
    D.index.name = "gene"
    return D


def analyse(D, genes, paths_src, P):
    """主要 4 量を返す。paths_src は (collection, name, 遺伝子名の集合)。"""
    sub = D.loc[D.index.intersection(genes)].dropna(axis=0, how="any")
    pos = {g: i for i, g in enumerate(sub.index)}
    X = sub.to_numpy(float)
    C = X - X.mean(axis=0)
    groups = []
    for _c, _n, names in paths_src:
        gi = np.array(sorted(pos[g] for g in names if g in pos))
        if len(gi) >= 30:
            groups.append(gi)
    idx = {s: i for i, s in enumerate(sub.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()

    def a(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        c = (Mn.T @ Mn)[ia, ib]
        x, y = c[m & wi], c[m & ~wi]
        return float(np.mean((x[:, None] > y[None, :]) + 0.5 * (x[:, None] == y[None, :])))

    agg_r = np.vstack([X[g].mean(axis=0) for g in groups])
    agg_c = np.vstack([C[g].mean(axis=0) for g in groups])
    mu = X.mean(axis=0)
    is_cat = np.array([s.startswith("cat") for s in sub.columns])
    return {"n_genes": X.shape[0], "n_sets": len(groups),
            "auc_gene": a(X), "auc_agg_uncentered": a(agg_r), "auc_agg_centered": a(agg_c),
            "centering_drop": a(agg_r) - a(agg_c),
            "mean_separates": bool(mu[is_cat].max() < mu[~is_cat].min()),
            "mean_gap": float(mu[~is_cat].min() - mu[is_cat].max())}


def main() -> int:
    D0, P = PCS.load()
    genes = set(D0.index)
    # 経路は遺伝子名で持つ（シャッフルで行位置が変わるため）
    names = list(D0.index)
    paths_src = [(c, n, {names[i] for i in gi}) for c, n, gi in PCS.collect(D0, PCS.COLL_FULL)]

    print("元空間の Δ を計算中…")
    sd = source_deltas()
    ga = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}

    base = analyse(build(sd), ga, paths_src, P)
    print(f"現行（先頭を採用）: 遺伝子 {base['n_genes']} / 集合 {base['n_sets']} / "
          f"gene {base['auc_gene']:.4f} / 非中心化 {base['auc_agg_uncentered']:.4f} / "
          f"中心化 {base['auc_agg_centered']:.4f} / 差 {base['centering_drop']:+.4f}")

    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(N_SHUFFLE):
        r = analyse(build(sd, rng), ga, paths_src, P)
        r["shuffle"] = i
        rows.append(r)
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{N_SHUFFLE} 回")
    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(6).to_csv(OUT / "duplicate_symbol_shuffles.tsv", sep="\t", index=False)

    summ = []
    print(f"\n{'quantity':22s} {'current':>9s} {'median':>9s} {'min':>9s} {'max':>9s} "
          f"{'max |diff|':>10s}")
    for k in ("n_genes", "auc_gene", "auc_agg_uncentered", "auc_agg_centered",
              "centering_drop", "mean_gap"):
        v = T[k].astype(float)
        d = float((v - base[k]).abs().max())
        summ.append({"quantity": k, "current": base[k], "median": v.median(),
                     "min": v.min(), "max": v.max(), "max_abs_diff_from_current": d,
                     "n_shuffles": len(T)})
        print(f"{k:22s} {base[k]:9.4f} {v.median():9.4f} {v.min():9.4f} {v.max():9.4f} "
              f"{d:10.4f}")
    sep = T["mean_separates"].sum()
    summ.append({"quantity": "mean_separates_species", "current": base["mean_separates"],
                 "median": f"{sep}/{len(T)}", "min": "", "max": "",
                 "max_abs_diff_from_current": "", "n_shuffles": len(T)})
    print(f"{'mean_separates':22s} {'True':>9s} {sep}/{len(T)} のシャッフルで分離を保つ")
    # 代表選択の対象になる遺伝子数。写像は状態ごとの遺伝子リストに対して行われるので、
    # 数えるのは「いずれかの状態で同じヒトシンボルに 2 行以上が写る」遺伝子。
    D0 = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    ga = set((HERE / "groupA_intersection.txt").read_text().split())
    complete = set(D0.dropna(axis=0, how="any").index) & ga
    dup_any, per_sp = set(), {}
    for _state, (sp, s) in sd.items():
        h = human_of(s.index, sp)
        vc = pd.Series(h).value_counts()
        dup = set(vc[vc > 1].index)
        dup_any |= dup
        per_sp.setdefault(sp, set()).update(dup & complete)
    affected = sorted(dup_any & complete)
    for sp, g in sorted(per_sp.items()):
        summ.append({"quantity": f"analysed genes with duplicate source rows, {sp}",
                     "current": len(g), "median": "", "min": "", "max": "",
                     "max_abs_diff_from_current": "", "n_shuffles": ""})
    summ.append({"quantity": "analysed genes with duplicate source rows",
                 "current": len(affected), "median": "", "min": "", "max": "",
                 "max_abs_diff_from_current": "", "n_shuffles": ""})
    summ.append({"quantity": "share of the analysed genes",
                 "current": round(len(affected) / len(complete), 6), "median": "", "min": "",
                 "max": "", "max_abs_diff_from_current": "", "n_shuffles": ""})
    print(f"  代表選択の対象: {len(affected)} / {len(complete)} = "
          f"{len(affected)/len(complete):.2%}  内訳 "
          + " / ".join(f"{k} {len(v)}" for k, v in sorted(per_sp.items())))
    pd.DataFrame(summ).to_csv(OUT / "duplicate_symbol_sensitivity.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'duplicate_symbol_sensitivity.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
