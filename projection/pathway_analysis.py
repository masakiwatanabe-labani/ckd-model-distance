"""O の解析: 経路粒度 cos θ から3つの問いに答える。

  1. 全体 cos が低いマウス状態でも、特定の経路では高いものがあるか
  2. 経路によってモデルの順位が入れ替わるか
  3. 旧稿の「モジュール集約で種の分離が消える」との関係
     — 集約という操作そのものを S1 空間で再現し、経路単位の cos と比べる

多重比較のため個別経路を「有意」とは呼ばない。分布とパターンのみ記述する。
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from pathway_cos import read_gmt, COLLECTIONS, MIN_GENES, REF  # noqa: E402

OUT = HERE / "results" / "pathway"
CAT = ["cat_CKD12", "cat_med_CKD12", "cat_med_CKD34"]      # 基準軸そのものは除く


def main():
    R = pd.read_csv(OUT / "pathway_cos.tsv", sep="\t")
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes, states = D.index, list(D.columns)
    glob = pd.read_csv(HERE / "results" / "alpha_groupA" / "projection.tsv",
                       sep="\t", index_col=0)["cos"]
    mouse = [s for s in states if not s.startswith("cat")]

    # -------------------------------------------------- 0. 遺伝子数と cos のばらつき
    NULL = pd.read_csv(OUT / "random_set_null.tsv", sep="\t")
    band = NULL.groupby("n_genes").iqr.median()
    print("=== 遺伝子数と cos のばらつき（ランダム集合 300 回の IQR 中央値） ===")
    for n, v in band.items():
        print(f"  n={n:4d}: IQR {v:.3f}")
    print(f"  経路の遺伝子数: 中央値 {R.n_genes.median():.0f}, "
          f"範囲 [{R.n_genes.min()}, {R.n_genes.max()}]")

    # -------------------------------------------------- 1. 全体が低い状態の高 cos 経路
    print("\n=== 1. 全体 cos が低いマウス状態で、経路単位では高いものがあるか ===")
    rows = []
    for s in mouse:
        d = R[R.state == s]
        g = float(glob[s])
        hi = d[(d.cos > d.random_hi95) & (d.n_genes >= 50)]
        rows.append({"state": s, "global_cos": g, "n_pathways": len(d),
                     "pathway_cos_median": d.cos.median(), "pathway_cos_max": d.cos.max(),
                     "n_above_random_band": len(hi),
                     "max_pathway": d.loc[d.cos.idxmax(), "pathway"][:44],
                     "max_n_genes": int(d.loc[d.cos.idxmax(), "n_genes"])})
    S1 = pd.DataFrame(rows).sort_values("global_cos")
    S1.round(3).to_csv(OUT / "state_pathway_summary.tsv", sep="\t", index=False)
    print(S1.round(3).to_string(index=False))

    # 全体最低の2状態で、猫の全体 cos (0.907) を超える経路があるか
    thr = float(glob["cat_CKD12"])
    exceed = R[(R.state.isin(mouse)) & (R.cos > thr) & (R.n_genes >= 50)]
    print(f"\n  マウス状態 x 経路のうち cos > {thr:.3f}（猫 CKD1/2 の全体値）: "
          f"{len(exceed)} / {len(R[R.state.isin(mouse)])} セル")
    if len(exceed):
        print(exceed.nlargest(8, "cos")[["state", "pathway", "n_genes", "cos",
                                         "ceiling", "random_hi95"]].round(3).to_string(index=False))

    # -------------------------------------------------- 2. 経路ごとのモデル順位の入れ替わり
    print("\n=== 2. 経路によってマウス状態の順位が入れ替わるか ===")
    piv = R[R.state.isin(mouse)].pivot_table(index="pathway", columns="state", values="cos")
    piv = piv.dropna()
    ranks = piv.rank(axis=1)
    # Kendall の W（一致係数）
    n_p, k = ranks.shape
    Rj = ranks.sum(axis=0).to_numpy()
    W = 12 * ((Rj - Rj.mean()) ** 2).sum() / (n_p ** 2 * (k ** 3 - k))
    print(f"  Kendall の W（{n_p} 経路 x {k} マウス状態）= {W:.3f}  "
          f"（1 = 全経路で順位が同一、0 = 無関係）")
    top = piv.idxmax(axis=1).value_counts()
    print("  各経路で cos が最大になったマウス状態:")
    for s, c in top.items():
        print(f"    {s:12s} {c:4d} 経路 ({c/n_p:.0%})")
    pair_rho = [stats.spearmanr(ranks.iloc[i], ranks.iloc[j])[0]
                for i, j in combinations(range(min(n_p, 60)), 2)]
    print(f"  経路ペア間の順位相関: 中央 {np.median(pair_rho):+.3f} "
          f"[{np.percentile(pair_rho,5):+.3f}, {np.percentile(pair_rho,95):+.3f}]")

    # -------------------------------------------------- 3. 集約 vs 経路単位
    print("\n=== 3. 旧稿の『集約で種の分離が消える』との関係 ===")
    # (a) 経路単位: 各経路で 猫3状態 vs マウス12状態 の cos を AUC で比較
    aucs = []
    for pth, d in R.groupby("pathway"):
        d = d.set_index("state")
        c = d.loc[[x for x in CAT if x in d.index], "cos"].to_numpy()
        m = d.loc[[x for x in mouse if x in d.index], "cos"].to_numpy()
        if len(c) < 2 or len(m) < 5:
            continue
        a = ((c[:, None] > m[None, :]).sum() + 0.5 * (c[:, None] == m[None, :]).sum()) / (c.size * m.size)
        aucs.append({"pathway": pth, "n_genes": int(d.n_genes.iloc[0]), "auc_cat_over_mouse": float(a),
                     "cat_median": float(np.median(c)), "mouse_median": float(np.median(m))})
    A = pd.DataFrame(aucs)
    A.round(4).to_csv(OUT / "pathway_species_separation.tsv", sep="\t", index=False)
    print(f"  経路単位の分離 AUC（猫 > マウス）: 中央 {A.auc_cat_over_mouse.median():.3f}, "
          f"範囲 [{A.auc_cat_over_mouse.min():.3f}, {A.auc_cat_over_mouse.max():.3f}]")
    print(f"    AUC = 1.00（完全分離）の経路: {int((A.auc_cat_over_mouse >= 0.999).sum())} / {len(A)} "
          f"({(A.auc_cat_over_mouse >= 0.999).mean():.0%})")
    print(f"    AUC >= 0.90 の経路          : {int((A.auc_cat_over_mouse >= 0.90).sum())} / {len(A)} "
          f"({(A.auc_cat_over_mouse >= 0.90).mean():.0%})")
    print(f"    AUC <= 0.60 の経路          : {int((A.auc_cat_over_mouse <= 0.60).sum())} / {len(A)} "
          f"({(A.auc_cat_over_mouse <= 0.60).mean():.0%})")

    # (b) 集約: 各経路を1つのスコアに潰し、193次元の経路スコア空間で状態間 cos
    paths = []
    gpos = {g: i for i, g in enumerate(genes)}
    for coll, fn in COLLECTIONS.items():
        for name, members in read_gmt(B.REF / fn).items():
            idx = np.array(sorted(gpos[g] for g in members & set(genes)))
            if len(idx) >= MIN_GENES:
                paths.append((name, idx))
    X = D.to_numpy(float)
    Z = (X - X.mean(axis=0)) / X.std(axis=0)          # 状態内で標準化してから集約
    M = np.vstack([Z[idx].mean(axis=0) for _, idx in paths])   # (経路, 状態)
    agg = {}
    for i, j in combinations(range(len(states)), 2):
        u, v = M[:, i], M[:, j]
        agg[(states[i], states[j])] = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
    ai = pd.DataFrame([{"a": a, "b": b, "cos_aggregated": c,
                        "class": "within species" if (a.startswith("cat") == b.startswith("cat"))
                        else "cross species"} for (a, b), c in agg.items()])
    ai.round(4).to_csv(OUT / "aggregated_module_cos.tsv", sep="\t", index=False)
    w = ai[ai["class"] == "within species"].cos_aggregated
    x = ai[ai["class"] == "cross species"].cos_aggregated
    ov = ((w.to_numpy()[:, None] > x.to_numpy()[None, :]).mean())
    print(f"\n  集約（{len(paths)} 経路スコア空間）での状態間 cos:")
    print(f"    種内 (n={len(w)}): 中央 {w.median():.3f} [{w.min():.3f}, {w.max():.3f}]")
    print(f"    種間 (n={len(x)}): 中央 {x.median():.3f} [{x.min():.3f}, {x.max():.3f}]")
    print(f"    AUC(種内 > 種間) = {ov:.3f}")
    print(f"  （参考）遺伝子レベルの同じ量: 種内 0.279 / 種間 0.552 の非類似度 → AUC は L の §3 参照")
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
