"""対照解析: 基準軸を固定せず、全状態ペアの cos 分布の中で 0.907 を読む。

問い:
    cat_CKD12 x cat_CKD34 の cos = 0.907 は「猫CKD1/2は方向的に猫CKD3/4に近い」
    ことを示すのか、それとも同一データセット・同一対照群から計算した Δ どうしなら
    自動的に高く出るだけなのか。

ペアの分類:
    within_dataset            同一データセット内（cat / podtreck / iri それぞれ）
    same_species_diff_dataset 同種・異データセット（Pod-TRECK x IRI）
    cross_species             異種

対照:
    上振れ  GSE98622 内の時点ペア全て（同一データセット・同一種・異なる時点）
    下振れ  猫 cortex x 猫 medulla（同一データセット・同一個体・異区画）

            旧稿の値の対応に注意:
              CKD3/4 vs Control 同士の 皮質x髄質 rho = 0.763（06_compartment 晩期）
              進行軸(CKD3/4 vs CKD1/2)同士の 皮質x髄質 rho = 0.204（06_compartment）
            「進行軸が直交」は後者の話。Δ行列の状態（対照群比）で作れるのは前者なので、
            進行軸は delta_matrix_aux.tsv から別途読み、両方を下振れ対照として出す。

共有対照バイアス:
    2つの状態を同じ対照群サンプルに対する差として計算すると、対照側の
    サンプリング誤差が両方の Δ に同じ符号で乗り、cos が押し上げられる。
    cat_CKD12 と cat_CKD34 はまさにこの形（同じ6頭の Control）。
    分割対照解析で、対照を重ならない半分ずつに分けて計算し直し、
    上振れ対照にも同じ処理を施して like-for-like で比較する。
"""
from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

# 遺伝子空間の差し替え（AA. 感度解析）。既定は Group A intersection で、
# 環境変数を与えない限り従来と完全に同じ動作になる。
TAG = os.environ.get("XSP_TAG", "")
GENES = Path(os.environ.get("XSP_GENES", str(HERE / "groupA_intersection.txt")))
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402

OUT = HERE / "results" / ("control_analysis" + TAG)
OUT.mkdir(parents=True, exist_ok=True)

DATASET = {"cat_CKD12": "cat", "cat_CKD34": "cat",
           "cat_med_CKD12": "cat", "cat_med_CKD34": "cat",
           "mouse_5D": "podtreck", "mouse_2W": "podtreck", "mouse_3W": "podtreck"}
DATASET.update({f"IRI_{t}": "iri" for t in
                ["2h", "4h", "24h", "48h", "72h", "7d", "14d", "28d", "12mo"]})
SPECIES = {k: ("cat" if v == "cat" else "mouse") for k, v in DATASET.items()}
# 同じ対照群サンプルに対する差として計算されている状態のグループ
CTRL_GROUP = {"cat_CKD12": "cat_ctx_ctrl", "cat_CKD34": "cat_ctx_ctrl",
              "cat_med_CKD12": "cat_med_ctrl", "cat_med_CKD34": "cat_med_ctrl",
              "mouse_5D": "podtreck_ctrl", "mouse_2W": "podtreck_ctrl",
              "mouse_3W": "podtreck_ctrl", "IRI_12mo": "sham12m"}
CTRL_GROUP.update({f"IRI_{t}": "young_sham" for t in
                   ["2h", "4h", "24h", "48h", "72h", "7d", "14d", "28d"]})
COMPARTMENT = {"cat_CKD12": "ctx", "cat_CKD34": "ctx",
               "cat_med_CKD12": "med", "cat_med_CKD34": "med"}
STAGE = {"cat_CKD12": "CKD12", "cat_CKD34": "CKD34",
         "cat_med_CKD12": "CKD12", "cat_med_CKD34": "CKD34"}


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def cat_pair_subclass(a, b):
    if COMPARTMENT[a] == COMPARTMENT[b]:
        return "same_compartment_diff_stage"
    return ("same_stage_diff_compartment" if STAGE[a] == STAGE[b]
            else "diff_compartment_diff_stage")


# ------------------------------------------------------------- 分割対照
def split_control_deltas(mat, grp, ctrl_label, case_labels, genes, species):
    """対照サンプルを重ならない半分に分けた Δ を、全ての分け方について作る。

    返り値: {case_label: {subset_index: Δ(Series, genes に整列)}}, 補集合の対応表
    """
    ctrl_cols = list(grp[grp.isin(ctrl_label)].index if isinstance(ctrl_label, list)
                     else grp[grp == ctrl_label].index)
    n = len(ctrl_cols)
    half = n // 2
    subsets = [frozenset(c) for c in combinations(range(n), half)]
    comp = {s: frozenset(set(range(n)) - s) for s in subsets}
    sub_idx = {s: i for i, s in enumerate(subsets)}

    ctrl_means = {}
    for s in subsets:
        cols = [ctrl_cols[i] for i in sorted(s)]
        ctrl_means[s] = mat[cols].mean(axis=1)

    out = {}
    for lab in case_labels:
        case_mean = mat[grp[grp == lab].index].mean(axis=1)
        d = {}
        for s in subsets:
            v = case_mean - ctrl_means[s]
            v = v[np.isfinite(v)]
            hv, _ = B.to_human(v, species)
            d[sub_idx[s]] = hv.reindex(genes)
        out[lab] = d
    return out, {sub_idx[s]: sub_idx[comp[s]] for s in subsets}


def split_control_cos(deltas, comp_map, a_lab, b_lab):
    """状態 a と b を、重ならない対照半分に対して計算した cos の分布。"""
    vals = []
    for si, ci in comp_map.items():
        u, v = deltas[a_lab][si], deltas[b_lab][ci]
        m = u.notna() & v.notna()
        if m.sum() > 50:
            vals.append(cosine(u[m].to_numpy(), v[m].to_numpy()))
    return np.array(vals, float)


def main():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes = D.index
    states = list(D.columns)
    X = D.to_numpy(float)
    print(f"cos 行列: {len(states)} 状態 x {len(genes)} 遺伝子（Group A intersection, complete case）")

    # ---------------------------------------------------------- 行列とペア表
    C = pd.DataFrame(np.eye(len(states)), index=states, columns=states)
    R = pd.DataFrame(np.eye(len(states)), index=states, columns=states)
    rows = []
    for i, j in combinations(range(len(states)), 2):
        a, b = states[i], states[j]
        c = cosine(X[:, i], X[:, j])
        r = float(stats.spearmanr(X[:, i], X[:, j])[0])
        C.loc[a, b] = C.loc[b, a] = c
        R.loc[a, b] = R.loc[b, a] = r
        same_ds = DATASET[a] == DATASET[b]
        cls = ("within_dataset" if same_ds else
               "same_species_diff_dataset" if SPECIES[a] == SPECIES[b] else "cross_species")
        rows.append({"a": a, "b": b, "cos": c, "spearman_rho": r, "class": cls,
                     "dataset_a": DATASET[a], "dataset_b": DATASET[b],
                     "shared_control": CTRL_GROUP[a] == CTRL_GROUP[b],
                     "cat_subclass": cat_pair_subclass(a, b) if same_ds and DATASET[a] == "cat" else ""})
    P = pd.DataFrame(rows)
    C.round(6).to_csv(OUT / "cos_matrix.tsv", sep="\t")
    R.round(6).to_csv(OUT / "spearman_matrix.tsv", sep="\t")
    P.round(6).to_csv(OUT / "pairs.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 分布
    def describe(d, label):
        v = d["cos"].to_numpy(float)
        v = v[np.isfinite(v)]
        return {"set": label, "n_pairs": len(v), "min": v.min(), "q25": np.percentile(v, 25),
                "median": np.median(v), "q75": np.percentile(v, 75), "max": v.max()}

    dist = [describe(P[P["class"] == c], c) for c in
            ["within_dataset", "same_species_diff_dataset", "cross_species"]]
    dist += [describe(P[(P["class"] == "within_dataset") & (P.dataset_a == ds)], f"within_{ds}")
             for ds in ["cat", "podtreck", "iri"]]
    dist += [describe(P[P.shared_control], "shared_control_pairs"),
             describe(P[~P.shared_control], "no_shared_control_pairs")]
    DIST = pd.DataFrame(dist)
    DIST.round(4).to_csv(OUT / "cos_distributions.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 対照と検定対象
    def pick(a, b, col="cos"):
        m = ((P.a == a) & (P.b == b)) | ((P.a == b) & (P.b == a))
        return float(P.loc[m, col].iloc[0])

    test_cos = pick("cat_CKD12", "cat_CKD34")
    iri = P[(P["class"] == "within_dataset") & (P.dataset_a == "iri")]["cos"].to_numpy(float)
    within = P[P["class"] == "within_dataset"]["cos"].to_numpy(float)

    # 下振れ対照（進行軸）は補助行列から
    A = pd.read_csv(HERE / "delta_matrix_aux.tsv", sep="\t", index_col=0).reindex(genes)
    am = A.notna().all(axis=1)
    prog_cos = cosine(A.loc[am, "cat_ctx_prog"].to_numpy(), A.loc[am, "cat_med_prog"].to_numpy())
    prog_rho = float(stats.spearmanr(A.loc[am, "cat_ctx_prog"], A.loc[am, "cat_med_prog"])[0])

    def pct(v, arr):
        return float((arr < v).mean() * 100)

    ctrls = pd.DataFrame([
        {"role": "検定対象", "pair": "cat_CKD12 x cat_CKD34",
         "cos": test_cos, "spearman_rho": pick("cat_CKD12", "cat_CKD34", "spearman_rho")},
        {"role": "上振れ対照(中央値)", "pair": "IRI 時点ペア 36組",
         "cos": float(np.median(iri)), "spearman_rho": np.nan},
        {"role": "上振れ対照(最大)", "pair": "IRI 時点ペア 36組",
         "cos": float(iri.max()), "spearman_rho": np.nan},
        {"role": "下振れ対照(指示どおり)", "pair": "cat_CKD34 x cat_med_CKD34",
         "cos": pick("cat_CKD34", "cat_med_CKD34"),
         "spearman_rho": pick("cat_CKD34", "cat_med_CKD34", "spearman_rho")},
        {"role": "下振れ対照(進行軸)", "pair": "cat_ctx_prog x cat_med_prog",
         "cos": prog_cos, "spearman_rho": prog_rho},
    ])
    ctrls["percentile_in_within_dataset"] = [pct(c, within) for c in ctrls["cos"]]
    ctrls["percentile_in_IRI_pairs"] = [pct(c, iri) for c in ctrls["cos"]]
    ctrls.round(4).to_csv(OUT / "controls.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 分割対照
    cat_mat, cat_grp = B._cat_mat("ctx")
    cat_d, cat_comp = split_control_deltas(cat_mat, cat_grp, "Control",
                                           ["CKD1/2", "CKD3/4"], genes, "cat")
    cat_split = split_control_cos(cat_d, cat_comp, "CKD1/2", "CKD3/4")

    imat = B.read_int("mouse_iri_matrix")
    igrp = B.read_int("mouse_iri_grp").squeeze()
    ikeep = np.log2(imat + 1)
    young = [s for s, _l, _h, c in B.IRI_STATES if c == B.YOUNG_SHAM]
    ilabels = [l for _s, l, _h, c in B.IRI_STATES if c == B.YOUNG_SHAM]
    iri_d, iri_comp = split_control_deltas(ikeep, igrp, B.YOUNG_SHAM, ilabels, genes, "mouse")
    iri_split = {}
    for x, y in combinations(range(len(ilabels)), 2):
        iri_split[(young[x], young[y])] = split_control_cos(
            iri_d, iri_comp, ilabels[x], ilabels[y])
    iri_split_med = np.array([np.median(v) for v in iri_split.values()])

    sc = pd.DataFrame([
        {"pair": "cat_CKD12 x cat_CKD34", "n_ctrl_full": 6, "n_ctrl_split": 3,
         "cos_shared_control": test_cos,
         "cos_split_median": float(np.median(cat_split)),
         "cos_split_lo": float(np.percentile(cat_split, 2.5)),
         "cos_split_hi": float(np.percentile(cat_split, 97.5)),
         "n_splits": len(cat_split)},
    ] + [
        {"pair": f"{a} x {b}", "n_ctrl_full": 6, "n_ctrl_split": 3,
         "cos_shared_control": pick(a, b),
         "cos_split_median": float(np.median(v)),
         "cos_split_lo": float(np.percentile(v, 2.5)),
         "cos_split_hi": float(np.percentile(v, 97.5)),
         "n_splits": len(v)}
        for (a, b), v in iri_split.items()])
    sc["inflation"] = sc["cos_shared_control"] - sc["cos_split_median"]
    sc.round(4).to_csv(OUT / "split_control.tsv", sep="\t", index=False)

    cat_split_med = float(np.median(cat_split))
    pct_split = float((iri_split_med < cat_split_med).mean() * 100)

    # ---------------------------------------------------------- 報告
    L = ["# 対照解析: cos 分布の中での cat_CKD12 x cat_CKD34", "",
         f"Group A intersection, 全16状態で complete case の {len(genes)} 遺伝子。",
         "基準軸は固定していない（全 120 ペアの cos を対称に計算）。", "",
         "## クラス別 cos 分布", "",
         "| 集合 | ペア数 | 最小 | Q25 | 中央値 | Q75 | 最大 |", "|---|---|---|---|---|---|---|"]
    for _, r in DIST.iterrows():
        L.append(f"| {r['set']} | {r['n_pairs']} | {r['min']:.3f} | {r['q25']:.3f} "
                 f"| {r['median']:.3f} | {r['q75']:.3f} | {r['max']:.3f} |")
    L += ["", "## 対照と検定対象", "",
          "| 役割 | ペア | cos | Spearman rho | within_dataset 内の百分位 | IRIペア内の百分位 |",
          "|---|---|---|---|---|---|"]
    for _, r in ctrls.iterrows():
        rr = "-" if not np.isfinite(r["spearman_rho"]) else f"{r['spearman_rho']:.3f}"
        L.append(f"| {r['role']} | {r['pair']} | {r['cos']:.3f} | {rr} "
                 f"| {r['percentile_in_within_dataset']:.0f} | {r['percentile_in_IRI_pairs']:.0f} |")
    L += ["", "## ネコ内ペアの内訳", "",
          "| ペア | 種別 | cos | Spearman rho |", "|---|---|---|---|"]
    for _, r in P[(P["class"] == "within_dataset") & (P.dataset_a == "cat")].iterrows():
        L.append(f"| {r['a']} x {r['b']} | {r['cat_subclass']} | {r['cos']:.3f} | {r['spearman_rho']:.3f} |")
    L += ["", "## 分割対照（共有対照バイアスの除去）", "",
          "対照サンプルを重ならない半分に分け、2つの状態を別々の半分に対する差として",
          "計算し直した cos。上振れ対照（IRI 若齢sham を対照とする8時点）にも",
          "同じ処理を施してある（どちらも対照 n=6 -> 3）。", "",
          f"- cat_CKD12 x cat_CKD34: 共有対照 {test_cos:.3f} -> 分割対照 中央値 "
          f"{cat_split_med:.3f} [95% {np.percentile(cat_split, 2.5):.3f}, "
          f"{np.percentile(cat_split, 97.5):.3f}]（{len(cat_split)} 通り）",
          f"- IRI 時点ペア 28組の分割対照 cos 中央値: 最小 {iri_split_med.min():.3f} / "
          f"中央値 {np.median(iri_split_med):.3f} / 最大 {iri_split_med.max():.3f}",
          f"- 分割対照どうしで比べたとき、cat_CKD12 x cat_CKD34 は IRI 28 ペア中の "
          f"第 {pct_split:.0f} 百分位", "",
          "分割対照では各 Δ の対照 n が半分になるので、独立ノイズが増えて cos は",
          "全体に下がる。両側に同じ処理をしているので比較は成立するが、",
          "分割後の絶対値を共有対照の値と直接比べてはいけない。"]
    (OUT / "report.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
