"""R: 経路単位の種分離を §2.4 と同一設計に揃えて計算し直す。

旧版の問題:
    経路単位は「猫3状態 vs マウス12状態を基準軸への cos θ で分離」する設計だった。
    猫3状態は基準軸とコホート・対照を共有するので、§2.4 が定量した共有対照バイアス
    （対照共有ペア中央 0.717 対 非共有 0.494）をそのまま受ける。
    しかも集約後の AUC 0.541 はペア単位・ペア間 cos θ で、単位も量も違っていた。
    比較になっていない。

新設計（経路単位と集約で完全に揃える）:
    1. 各経路内で全120ペアの cos θ を計算する（基準軸を固定しない）
    2. §2.4 と同じ制限をかける
         (a) 対照検体を共有するペアを除外
         (b) さらに厳しく、同一コホート由来のペアも除外
    3. 種内ペア vs 種間ペアの AUC を経路ごとに出す
    4. 集約（各経路を1スコアに潰した193次元空間）でも**同じ制限**で AUC を出す
"""
from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]

# 遺伝子空間の差し替え（AA. 感度解析）。既定は Group A intersection で、
# 環境変数を与えない限り従来と完全に同じ動作になる。
TAG = os.environ.get("XSP_TAG", "")
GENES = Path(os.environ.get("XSP_GENES", str(HERE / "groupA_intersection.txt")))
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathway_cos import read_gmt, COLLECTIONS, MIN_GENES  # noqa: E402

OUT = HERE / "results" / "revision1" / "pathway_reactome"
OUT.mkdir(parents=True, exist_ok=True)


def auc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if x.size < 3 or y.size < 3:
        return np.nan
    return float(((x[:, None] > y[None, :]).sum()
                  + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size))


def main():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes, states = D.index, list(D.columns)
    X = D.to_numpy(float)
    gpos = {g: i for i, g in enumerate(genes)}

    # ---- ペアの分類（§2.4 と同一の定義を control_analysis から引く）
    P = pd.read_csv(HERE / "results" / ("control_analysis" + TAG) / "pairs.tsv", sep="\t")
    P["species_a"] = np.where(P.a.str.startswith("cat"), "cat", "mouse")
    P["species_b"] = np.where(P.b.str.startswith("cat"), "cat", "mouse")
    P["within_species"] = P.species_a == P.species_b
    P["shared_cohort"] = ((P.dataset_a == "cat") & (P.dataset_b == "cat")) | P.shared_control
    idx = {s: i for i, s in enumerate(states)}
    ia = P.a.map(idx).to_numpy()
    ib = P.b.map(idx).to_numpy()
    wi = P.within_species.to_numpy()

    masks = {
        "all pairs (旧版と同じ制限なし)": np.ones(len(P), bool),
        "exclude shared controls": ~P.shared_control.to_numpy(),
        "exclude shared cohort": ~P.shared_cohort.to_numpy(),
    }
    for lab, m in masks.items():
        print(f"{lab}: 全 {int(m.sum())} ペア（種内 {int((m & wi).sum())} / 種間 {int((m & ~wi).sum())}）")

    paths = []
    for coll, fn in COLLECTIONS.items():
        for name, members in read_gmt(B.REF / fn).items():
            gi = np.array(sorted(gpos[g] for g in members & set(genes)))
            if len(gi) >= MIN_GENES:
                paths.append((coll, name, gi))
    print(f"経路 {len(paths)}")

    # ---- 1-3. 経路単位: 全120ペアの cos → 制限 → AUC
    rows = []
    for coll, name, gi in paths:
        Z = X[gi]
        Zn = Z / np.linalg.norm(Z, axis=0, keepdims=True)
        C = Zn.T @ Zn
        cos = C[ia, ib]
        rec = {"collection": coll, "pathway": name, "n_genes": len(gi)}
        for lab, m in masks.items():
            rec[f"auc[{lab}]"] = auc(cos[m & wi], cos[m & ~wi])
            rec[f"within_med[{lab}]"] = float(np.median(cos[m & wi]))
            rec[f"cross_med[{lab}]"] = float(np.median(cos[m & ~wi]))
        rows.append(rec)
    R = pd.DataFrame(rows)
    R.round(4).to_csv(OUT / "pathway_separation_paired.tsv", sep="\t", index=False)

    # ---- 4. 集約: 同じ制限で
    Zs = (X - X.mean(axis=0)) / X.std(axis=0)
    M = np.vstack([Zs[gi].mean(axis=0) for _, _, gi in paths])      # (経路, 状態)
    Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
    Cagg = Mn.T @ Mn
    cos_agg = Cagg[ia, ib]

    print("\n=== 経路単位 vs 集約（同一のペア集合・同一の量）===")
    print(f"{'制限':28s} {'経路単位 AUC 中央':>18s} {'>=0.90 の割合':>14s} {'集約 AUC':>10s}")
    summ = []
    for lab, m in masks.items():
        a_path = R[f"auc[{lab}]"].dropna()
        a_agg = auc(cos_agg[m & wi], cos_agg[m & ~wi])
        summ.append({"restriction": lab, "n_pairs": int(m.sum()),
                     "n_within": int((m & wi).sum()), "n_cross": int((m & ~wi).sum()),
                     "pathway_auc_median": float(a_path.median()),
                     "pathway_auc_q25": float(a_path.quantile(.25)),
                     "pathway_auc_q75": float(a_path.quantile(.75)),
                     "pathway_frac_ge_090": float((a_path >= 0.90).mean()),
                     "aggregated_auc": a_agg,
                     "within_med_agg": float(np.median(cos_agg[m & wi])),
                     "cross_med_agg": float(np.median(cos_agg[m & ~wi]))})
        print(f"{lab:28s} {a_path.median():18.3f} {(a_path>=0.90).mean():14.0%} {a_agg:10.3f}")
    S = pd.DataFrame(summ)
    S.round(4).to_csv(OUT / "separation_summary.tsv", sep="\t", index=False)

    print("\n=== 判定に使う対比（対照共有を除外した設計）===")
    lab = "exclude shared controls"
    a_path = R[f"auc[{lab}]"].dropna()
    a_agg = float(S[S.restriction == lab].aggregated_auc.iloc[0])
    print(f"  経路単位 AUC: 中央 {a_path.median():.3f}  四分位 [{a_path.quantile(.25):.3f}, {a_path.quantile(.75):.3f}]")
    print(f"                 >=0.90 {int((a_path>=0.90).sum())}/{len(a_path)} ({(a_path>=0.90).mean():.0%})"
          f"、<=0.60 {int((a_path<=0.60).sum())}/{len(a_path)}")
    print(f"  集約 AUC     : {a_agg:.3f}")
    print(f"  差           : {a_path.median()-a_agg:+.3f}")
    lab2 = "exclude shared cohort"
    a2 = R[f"auc[{lab2}]"].dropna(); a2g = float(S[S.restriction == lab2].aggregated_auc.iloc[0])
    print(f"\n  （より厳しい制限）経路単位 中央 {a2.median():.3f} / 集約 {a2g:.3f} / 差 {a2.median()-a2g:+.3f}")
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
