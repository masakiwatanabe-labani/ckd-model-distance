# -*- coding: utf-8 -*-
"""AE2: 経路集合の設計を変えたときの中心化の対比（査読者 3 major 4）。

(a) コレクション別（Reactome / GO BP / KEGG / Hallmark）
(b) 経路の冗長性を減らしたあと（遺伝子重複 Jaccard で足切りした代表集合）
(c) 4 つの遺伝子空間に共通する経路 ID だけに限った比較
    これで「遺伝子空間が変わった影響」と「資格を満たす経路集合が変わった影響」を分ける。
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
SPACES = [("Group A (2,016)", "groupA_intersection.txt"),
          ("All 1:1 orthologues (7,897)", "aa_ortholog_all.txt"),
          ("Matched Group B (1,632)", "aa_groupB_matched2.txt"),
          ("Matched Group A (1,632)", "aa_groupA_matched.txt")]
JACCARD = 0.5


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def load_space(gene_file):
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / gene_file).read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    P = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
    return D, P


def contrast(D, P, groups):
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()

    def a(M):
        A = np.vstack([M[g].mean(axis=0) for g in groups])
        An = A / np.linalg.norm(A, axis=0, keepdims=True)
        c = (An.T @ An)[ia, ib]
        return auc(c[m & wi], c[m & ~wi])
    u, cn = a(X), a(C)
    return u, cn, u - cn


def reduce_redundancy(paths, thresh=JACCARD):
    """大きい集合から順に、既採用と Jaccard が閾値未満のものだけ残す。"""
    order = sorted(range(len(paths)), key=lambda i: -len(paths[i][2]))
    kept, sets = [], []
    for i in order:
        s = set(paths[i][2].tolist())
        if all(len(s & k) / len(s | k) < thresh for k in sets):
            kept.append(i)
            sets.append(s)
    return [paths[i] for i in sorted(kept)]


def main() -> int:
    missing = [fn for fn in PCS.COLL_FULL.values() if not (PCS.REF / fn).exists()]
    if missing:
        print(f"GMT が無い: {missing}", file=sys.stderr)
        return 2
    rows = []
    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)

    # ---- (a) コレクション別
    print("(a) コレクション別（Group A）")
    for coll in ["Reactome", "GO_BP", "KEGG", "Hallmark"]:
        sub = [p for p in paths if p[0] == coll]
        if len(sub) < 5:
            print(f"  {coll}: {len(sub)} セットのみ。飛ばす")
            continue
        u, c, d = contrast(D, P, [gi for _c, _n, gi in sub])
        rows.append({"analysis": "collection", "label": coll, "n_sets": len(sub),
                     "auc_uncentered": u, "auc_centered": c, "difference": d})
        print(f"  {coll:10s} {len(sub):4d} セット  非中心化 {u:.4f} / 中心化 {c:.4f} / 差 {d:+.4f}")
    u, c, d = contrast(D, P, [gi for _c, _n, gi in paths])
    rows.append({"analysis": "collection", "label": "all four (main)", "n_sets": len(paths),
                 "auc_uncentered": u, "auc_centered": c, "difference": d})
    print(f"  {'all four':10s} {len(paths):4d} セット  非中心化 {u:.4f} / 中心化 {c:.4f} / 差 {d:+.4f}")

    # ---- (b) 冗長性の削減
    print(f"\n(b) 冗長性の削減（Jaccard < {JACCARD} の代表集合）")
    red = reduce_redundancy(paths)
    u, c, d = contrast(D, P, [gi for _c, _n, gi in red])
    rows.append({"analysis": "redundancy", "label": f"Jaccard < {JACCARD}", "n_sets": len(red),
                 "auc_uncentered": u, "auc_centered": c, "difference": d})
    print(f"  {len(paths)} → {len(red)} セット  非中心化 {u:.4f} / 中心化 {c:.4f} / 差 {d:+.4f}")
    for th in (0.3, 0.2):
        r2 = reduce_redundancy(paths, th)
        u2, c2, d2 = contrast(D, P, [gi for _c, _n, gi in r2])
        rows.append({"analysis": "redundancy", "label": f"Jaccard < {th}", "n_sets": len(r2),
                     "auc_uncentered": u2, "auc_centered": c2, "difference": d2})
        print(f"  Jaccard < {th}: {len(r2)} セット  非中心化 {u2:.4f} / 中心化 {c2:.4f} / 差 {d2:+.4f}")

    # ---- (c) 4 空間に共通する経路 ID
    print("\n(c) 4 遺伝子空間に共通する経路 ID だけに限る")
    per_space, qualified = {}, []
    for label, gf in SPACES:
        Ds, Ps = load_space(gf)
        ps = PCS.collect(Ds, PCS.COLL_FULL)
        per_space[label] = (Ds, Ps, {(c_, n_): gi for c_, n_, gi in ps})
        qualified.append({(c_, n_) for c_, n_, _g in ps})
        print(f"  {label:30s} 資格を満たす経路 {len(ps)}")
    common = set.intersection(*qualified)
    print(f"  共通する経路 ID: {len(common)}")
    for label, _gf in SPACES:
        Ds, Ps, mp = per_space[label]
        gs_all = [mp[k] for k in mp]
        gs_com = [mp[k] for k in mp if k in common]
        ua, ca, da = contrast(Ds, Ps, gs_all)
        uc, cc, dc = contrast(Ds, Ps, gs_com)
        rows.append({"analysis": "common sets", "label": label, "n_sets": len(gs_com),
                     "auc_uncentered": uc, "auc_centered": cc, "difference": dc,
                     "n_sets_all": len(gs_all), "difference_all_sets": da})
        print(f"  {label:30s} 全 {len(gs_all):4d} セット 差 {da:+.4f}  →  "
              f"共通 {len(gs_com)} セット 差 {dc:+.4f}（非中心化 {uc:.4f} / 中心化 {cc:.4f}）")

    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "pathway_set_design.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'pathway_set_design.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
