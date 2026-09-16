"""Task 4. 猫側にフィルタを課していないことの影響を確認する。

(1) Group A 2,016 遺伝子の中で、猫とマウスの対照群発現の分位分布を比べる。
(2) 猫側にもマウスと同等の発現フィルタを課したとき、4つのネコ状態の
    split-half 信頼性がどう変わるかを計算する。

既存の結果は書き換えず results/revision1/ に出す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from reliability import cosine  # noqa: E402

OUT = HERE / "results" / "revision1"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260826
N_SPLIT = 500


def group_a_genes():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    return list(D.index)


def expression_quantiles(genes):
    """対照群平均発現の、データセット内パーセンタイル順位を種ごとに比べる。"""
    E = pd.read_csv(HERE / "control_log2cpm.tsv", sep="\t", index_col=0)
    rows = []
    for col, species in [("cat_ctx_control", "cat"), ("podtreck_control", "mouse"),
                         ("iri_young_sham", "mouse")]:
        s = E[col].dropna()
        pct = s.rank(pct=True)
        inA = pct.reindex(genes).dropna()
        rows.append({"column": col, "species": species,
                     "n_all": len(s), "n_groupA": len(inA),
                     "pct_median": float(inA.median()),
                     "pct_q25": float(inA.quantile(.25)), "pct_q75": float(inA.quantile(.75)),
                     "frac_below_0.25": float((inA < .25).mean()),
                     "frac_below_0.50": float((inA < .50).mean())})
    return pd.DataFrame(rows)


def feline_reliability(min_expr):
    """猫の4状態について、発現フィルタ min_expr（log2 値の群平均最大）を課したときの
    分割半信頼性。min_expr が None ならフィルタなし（現行）。"""
    genes = group_a_genes()
    rng = np.random.default_rng(SEED)
    out = []
    for state, label, tissue in B.CAT_STATES:
        mat, grp = B._cat_mat(tissue)
        case = mat[grp[grp == label].index]
        ctrl = mat[grp[grp == "Control"].index]
        if min_expr is not None:
            gm = pd.concat([case.mean(axis=1), ctrl.mean(axis=1)], axis=1).max(axis=1)
            keep = gm >= min_expr
            case, ctrl = case[keep], ctrl[keep]
        nc, nk = case.shape[1], ctrl.shape[1]
        r = np.empty(N_SPLIT)
        for i in range(N_SPLIT):
            a = rng.permutation(nc); b = rng.permutation(nk)
            d1 = case.iloc[:, a[: nc // 2]].mean(axis=1) - ctrl.iloc[:, b[: nk // 2]].mean(axis=1)
            d2 = case.iloc[:, a[nc // 2:]].mean(axis=1) - ctrl.iloc[:, b[nk // 2:]].mean(axis=1)
            h1, _ = B.to_human(d1, "cat"); h2, _ = B.to_human(d2, "cat")
            u, v = h1.reindex(genes), h2.reindex(genes)
            mm = u.notna() & v.notna()
            r[i] = cosine(u[mm].to_numpy(), v[mm].to_numpy()) if mm.sum() > 50 else np.nan
        rh = float(np.nanmedian(r))
        out.append({"state": state, "filter": "none" if min_expr is None else f">={min_expr}",
                    "n_genes_kept": int(len(case)),
                    "n_groupA_finite": int(mm.sum()),
                    "r_half": rh, "r_SB": (2 * rh / (1 + rh)) if rh > 0 else 0.0})
    return pd.DataFrame(out)


def main():
    genes = group_a_genes()
    E = expression_quantiles(genes)
    E.round(4).to_csv(OUT / "groupA_expression_by_species.tsv", sep="\t", index=False)
    print("=== Group A 遺伝子の対照群発現パーセンタイル（各データセット内）===")
    print(E.round(3).to_string(index=False))

    # 猫のシート値は log2 スケール。マウスの FPKM>=1 に相当する水準を分位で合わせる
    mat, grp = B._cat_mat("ctx")
    allmean = mat.mean(axis=1)
    thr = float(np.percentile(allmean, 25))
    print(f"\n猫側フィルタの閾値（皮質 全遺伝子平均の 25 パーセンタイル）= {thr:.3f}")

    frames = [feline_reliability(None), feline_reliability(thr)]
    T = pd.concat(frames, ignore_index=True)
    T.round(4).to_csv(OUT / "feline_filter_sensitivity.tsv", sep="\t", index=False)
    print("\n=== ネコ4状態の分割半信頼性 ===")
    print(T.round(3).to_string(index=False))
    piv = T.pivot(index="state", columns="filter", values="r_half")
    print("\nr_half の変化:")
    print((piv.assign(delta=lambda d: d.iloc[:, 0] - d.iloc[:, 1])).round(3).to_string())
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
