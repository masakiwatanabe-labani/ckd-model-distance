# -*- coding: utf-8 -*-
"""C案の対照解析。centering なしの集約 AUC を、4つの遺伝子空間と 83 ペア部分集合で出す。

本文（Table 1）には centering ありの集約 AUC しか載っていないため、
「centering の差が Group A に固有かどうか」が読者に確かめられない。
同じ 4 空間と、コホートを共有する 4 ペアを除いた 83 ペアの部分集合について、
centering あり・なしの両方を計算する。

入力は Δ 行列・遺伝子空間のリスト・ペア表・GMT のみ。既存の結果ファイルは読まない。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402

OUT = HERE / "results" / "roundC"

# (表示名, 遺伝子リスト) — Table 1 の 4 列と同じ順
SPACES = [("Group A (2,016)", "groupA_intersection.txt"),
          ("All 1:1 orthologues (7,897)", "aa_ortholog_all.txt"),
          ("Matched Group B (1,632)", "aa_groupB_matched2.txt"),
          ("Matched Group A (1,632)", "aa_groupA_matched.txt")]

# コホート（同一のネコ個体）は共有するが対照は共有しないペア。83 ペア部分集合で外す。
COHORT_PAIRS = {("cat_CKD12", "cat_med_CKD12"), ("cat_CKD12", "cat_med_CKD34"),
                ("cat_CKD34", "cat_med_CKD12"), ("cat_CKD34", "cat_med_CKD34")}


def load_space(gene_file: str):
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / gene_file).read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    P = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
    return D, P


def conditions(D, P, paths, drop_cohort=False):
    """PCS と同じ四条件を、ペアの制限を変えて計算する。"""
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    agg_c = np.vstack([C[gi].mean(axis=0) for _c, _n, gi in paths])
    agg_r = np.vstack([X[gi].mean(axis=0) for _c, _n, gi in paths])
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    if drop_cohort:
        cohort = np.array([(a, b) in COHORT_PAIRS or (b, a) in COHORT_PAIRS
                           for a, b in zip(P.a, P.b)])
        m = m & ~cohort

    def cosmat(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        return (Mn.T @ Mn)[ia, ib]

    out = pd.DataFrame({"cos_raw": cosmat(X), "cos_centred": cosmat(C),
                        "cos_agg_centred": cosmat(agg_c), "cos_agg_raw": cosmat(agg_r),
                        "within": P.wi.astype(int).to_numpy()})
    return out.loc[m].reset_index(drop=True)


def auc(BA, col):
    a = BA[BA.within == 1][col].to_numpy(); b = BA[BA.within == 0][col].to_numpy()
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def pathway_median_auc(D, P, paths, drop_cohort=False):
    """経路ごとに AUC を出し、その中央値を返す（Table 1 の pathway-level AUC 行）。"""
    X = D.to_numpy(float)
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    wi = P.wi.to_numpy()
    m = (~P.shared_control).to_numpy()
    if drop_cohort:
        m = m & ~np.array([(a, b) in COHORT_PAIRS or (b, a) in COHORT_PAIRS
                           for a, b in zip(P.a, P.b)])

    def cosmat(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        return (Mn.T @ Mn)[ia, ib]

    out = []
    for _c, _n, gi in paths:
        c = cosmat(X[gi])
        a, b = c[m & wi], c[m & ~wi]
        out.append(float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :]))))
    return float(np.median(out))


def row(label, subset, BA, n_sets, n_genes, pw_median=float("nan"), collection="422 (+Reactome)"):
    return {"space": label, "pairs": subset, "collection": collection,
            "n_genes": n_genes, "n_sets": n_sets,
            "n_pairs": len(BA), "n_within": int(BA.within.sum()),
            "n_cross": int((1 - BA.within).sum()),
            "auc_gene_level": auc(BA, "cos_raw"),
            "auc_gene_centred": auc(BA, "cos_centred"),
            "auc_agg_centred": auc(BA, "cos_agg_centred"),
            "auc_agg_uncentred": auc(BA, "cos_agg_raw"),
            "pathway_auc_median": pw_median,
            "centring_drop": auc(BA, "cos_agg_raw") - auc(BA, "cos_agg_centred"),
            "cross_med_agg_centred": float(BA[BA.within == 0].cos_agg_centred.median()),
            "cross_med_agg_uncentred": float(BA[BA.within == 0].cos_agg_raw.median()),
            "within_med_agg_centred": float(BA[BA.within == 1].cos_agg_centred.median()),
            "within_med_agg_uncentred": float(BA[BA.within == 1].cos_agg_raw.median())}


def main() -> int:
    missing = [fn for fn in PCS.COLL_FULL.values() if not (PCS.REF / fn).exists()]
    if missing:
        print(f"遺伝子セットの GMT が無い: {missing}", file=sys.stderr)
        return 2

    rows = []
    for label, gf in SPACES:
        D, P = load_space(gf)
        paths = PCS.collect(D, PCS.COLL_FULL)
        BA = conditions(D, P, paths)
        rows.append(row(label, "87 (no shared controls)", BA, len(paths), D.shape[0],
                        pathway_median_auc(D, P, paths)))
        if gf == "groupA_intersection.txt":
            BA83 = conditions(D, P, paths, drop_cohort=True)
            rows.append(row(label, "83 (also no shared cohort)", BA83, len(paths), D.shape[0],
                            pathway_median_auc(D, P, paths, drop_cohort=True)))
            # 旧稿の 83 ペアの値は Reactome を含まない 193 セットで計算されていた。
            # どこから来た数字かが後で分かるように、その版も記録に残す。
            old = PCS.collect(D, PCS.COLL_BASE)
            BA83o = conditions(D, P, old, drop_cohort=True)
            rows.append(row(label, "83 (also no shared cohort)", BA83o, len(old), D.shape[0],
                            pathway_median_auc(D, P, old, drop_cohort=True),
                            collection="193 (GO+KEGG+Hallmark, superseded)"))

    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "centering_controls.tsv", sep="\t", index=False)

    print(f"{'space':28s} {'pairs':12s} {'sets':>5s} {'pathwayMed':>11s} "
          f"{'centred':>8s} {'uncentred':>10s} {'drop':>7s}")
    for r in T.itertuples():
        tag = " [193, 旧]" if r.collection.startswith("193") else ""
        print(f"{r.space:28s} {r.pairs[:12]:12s} {r.n_sets:5d} {r.pathway_auc_median:11.4f} "
              f"{r.auc_agg_centred:8.4f} {r.auc_agg_uncentred:10.4f} {r.centring_drop:+7.4f}{tag}")
    print(f"\n書き出し: {OUT / 'centering_controls.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
