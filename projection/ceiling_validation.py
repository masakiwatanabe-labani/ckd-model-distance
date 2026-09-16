"""G: 天井推定の妥当性検査。

補正後 cos が 1 を超えるペアが 120 中 5 組ある。真の cos は 1 を超えられないので、
天井 sqrt(r_XX * r_YY) が過小か、観測 cos が上振れしているかのどちらか。

検査:
  1. split-half の r を、分割に使った検体数 n の関数として出す（各 500 分割）。
     ただし disjoint な2半分を作るには 2n <= min(n_case, n_ctrl) が必要で、
     n=3 の群（Pod-TRECK 全状態、IRI 全状態の症例群）は n=1 しか取れない。
     症例側を n=1 に固定して対照側だけ変える副曲線も併せて出す。
  2. Spearman-Brown は r(kn) = k r(n) / (1 + (k-1) r(n)) を予言する。
     曲線が取れる状態で、r(1) からの予言と実測 r(2), r(3) を突き合わせる。
     予言より実測が低ければ SB は上振れ（= 天井を過大）、
     予言より高ければ SB は下振れ（= 天井を過小 → 補正後 >1 の説明になる）。
  3. 3通りの天井で「異種の観測最大 0.548 は天井から離れている」が成り立つか。
       ceiling_raw = sqrt(r_half_X * r_half_Y)   補正なし。最も保守的（低い天井）
       ceiling_SB  = sqrt(r_SB_X * r_SB_Y)       現行
       ceiling_sat = 飽和曲線の漸近値から          曲線が取れる状態のみ

補正後 >1 のもう一つの説明: 天井の式は2状態の誤差が独立であることを仮定する。
ネコの4状態は同一コホート（17頭）由来で、対照検体も共有する。誤差が相関していれば
観測 cos は「独立な誤差から計算した天井」を超えうる。異種ペアは動物も対照も
バッチも共有しないので、この仮定が最もよく成り立つ。
"""
from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

# 遺伝子空間の差し替え（AA. 感度解析）。既定は Group A intersection。
TAG = os.environ.get("XSP_TAG", "")
GENES = Path(os.environ.get("XSP_GENES", str(HERE / "groupA_intersection.txt")))
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from build_precision import human_projector  # noqa: E402
from reliability import state_specs, cosine  # noqa: E402

OUT = HERE / "results" / ("reliability" + TAG)
OUT.mkdir(parents=True, exist_ok=True)
N_SPLIT = 500
SEED = 20260826


def r_curve(mat, grp, case, ctrl, genes, species, rng, n_case=None, n_ctrl=None):
    """検体数を指定した disjoint 2分割の cos 分布。"""
    A = mat[grp[grp == case].index].to_numpy(float)
    C = mat[grp[grp.isin(ctrl)].index].to_numpy(float)
    human, keep = human_projector(mat.index, species)
    pos = pd.Series(np.arange(len(human)), index=human).reindex(genes).dropna()
    gi = pos.to_numpy(int)
    na, nc = A.shape[1], C.shape[1]
    if 2 * n_case > na or 2 * n_ctrl > nc:
        return None
    vals = np.empty(N_SPLIT)
    for i in range(N_SPLIT):
        pa, pc = rng.permutation(na), rng.permutation(nc)
        a1, a2 = A[:, pa[:n_case]], A[:, pa[n_case:2 * n_case]]
        c1, c2 = C[:, pc[:n_ctrl]], C[:, pc[n_ctrl:2 * n_ctrl]]
        d1 = (np.nanmean(a1, 1) - np.nanmean(c1, 1))[keep][gi]
        d2 = (np.nanmean(a2, 1) - np.nanmean(c2, 1))[keep][gi]
        m = np.isfinite(d1) & np.isfinite(d2)
        vals[i] = cosine(d1[m], d2[m]) if m.sum() > 50 else np.nan
    v = vals[np.isfinite(vals)]
    return float(np.median(v)) if v.size else np.nan


def sb(r, k=2):
    """Spearman-Brown: 長さ k 倍にしたときの信頼性。"""
    r = np.clip(r, 0, 0.999999)
    return k * r / (1 + (k - 1) * r)


def main():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes = D.index
    rng = np.random.default_rng(SEED)

    # ------------------------------------------------ 1. r(n) 曲線
    rows = []
    for state, mat, grp, case, ctrl, species in state_specs():
        na = int((grp == case).sum())
        nc = int(grp.isin(ctrl).sum())
        nmax = min(na, nc) // 2
        for n in range(1, nmax + 1):
            r = r_curve(mat, grp, case, ctrl, genes, species, rng, n_case=n, n_ctrl=n)
            rows.append({"state": state, "curve": "both", "n_case": n, "n_ctrl": n,
                         "r": r, "n_case_avail": na, "n_ctrl_avail": nc, "n_max": nmax})
        # 副曲線: 症例を 1 に固定して対照側だけ変える（n=3 群でも点が取れる）
        for k in range(1, nc // 2 + 1):
            r = r_curve(mat, grp, case, ctrl, genes, species, rng, n_case=1, n_ctrl=k)
            rows.append({"state": state, "curve": "ctrl_only", "n_case": 1, "n_ctrl": k,
                         "r": r, "n_case_avail": na, "n_ctrl_avail": nc, "n_max": nmax})
        print(f"{state:14s} case={na} ctrl={nc} n_max={nmax}")
    CUR = pd.DataFrame(rows)
    CUR.round(4).to_csv(OUT / "r_curve.tsv", sep="\t", index=False)

    # ------------------------------------------------ 2. SB の内部検証
    val = []
    both = CUR[CUR.curve == "both"]
    for state, d in both.groupby("state"):
        d = d.set_index("n_case")["r"]
        if 1 not in d.index:
            continue
        for k in [2, 3]:
            if k in d.index and np.isfinite(d[1]) and np.isfinite(d[k]):
                val.append({"state": state, "k": k, "r_1": d[1], "r_k_observed": d[k],
                            "r_k_SB_predicted": float(sb(d[1], k)),
                            "observed_minus_predicted": d[k] - float(sb(d[1], k))})
    VAL = pd.DataFrame(val)
    if len(VAL):
        VAL.round(4).to_csv(OUT / "sb_validation.tsv", sep="\t", index=False)

    # 対照側だけの曲線でも同じ検証（マウス状態を含められる）
    val2 = []
    co = CUR[CUR.curve == "ctrl_only"]
    for state, d in co.groupby("state"):
        d = d.set_index("n_ctrl")["r"]
        for k in [2, 3]:
            if 1 in d.index and k in d.index and np.isfinite(d[1]) and np.isfinite(d[k]):
                val2.append({"state": state, "k": k, "r_1": d[1], "r_k_observed": d[k],
                             "r_k_SB_predicted": float(sb(d[1], k)),
                             "observed_minus_predicted": d[k] - float(sb(d[1], k))})
    VAL2 = pd.DataFrame(val2)
    if len(VAL2):
        VAL2.round(4).to_csv(OUT / "sb_validation_ctrl_only.tsv", sep="\t", index=False)

    # ------------------------------------------------ 3. 3通りの天井
    REL = pd.read_csv(OUT / "reliability.tsv", sep="\t").set_index("state")
    P = pd.read_csv(OUT / "pairs_corrected.tsv", sep="\t")
    r_raw = REL["r_half_median"].clip(lower=0)
    r_sb = REL["reliability_SB_median"]
    # 飽和推定: n_max の実測 r を、全検体まで SB で外挿する（曲線が取れる状態のみ）
    r_sat = {}
    for state, d in both.groupby("state"):
        d = d.sort_values("n_case")
        nmax = int(d.n_max.iloc[0])
        if nmax >= 2 and np.isfinite(d.r.iloc[-1]):
            k = min(REL.loc[state, "n_case"], REL.loc[state, "n_ctrl"]) / nmax
            r_sat[state] = float(sb(d.r.iloc[-1], k))
    r_sat = pd.Series(r_sat)

    def ceil_col(rv, name):
        return P.apply(lambda x: float(np.sqrt(rv.get(x["a"], np.nan) * rv.get(x["b"], np.nan)))
                       if np.isfinite(rv.get(x["a"], np.nan)) and np.isfinite(rv.get(x["b"], np.nan))
                       else np.nan, axis=1).rename(name)

    P["ceiling_raw"] = ceil_col(r_raw, "ceiling_raw")
    P["ceiling_SB"] = ceil_col(r_sb, "ceiling_SB")
    P["ceiling_sat"] = ceil_col(r_sat, "ceiling_sat")
    # 同一コホート/同一対照検体の共有（天井の独立誤差仮定が破れる条件）
    P["shared_cohort"] = ((P.dataset_a == "cat") & (P.dataset_b == "cat")) | P.shared_control
    P.round(4).to_csv(OUT / "pairs_ceilings.tsv", sep="\t", index=False)

    S = []
    for cls in ["within_dataset", "same_species_diff_dataset", "cross_species"]:
        d = P[P["class"] == cls]
        row = {"class": cls, "n": len(d), "cos_median": float(d.cos.median()),
               "cos_max": float(d.cos.max())}
        for c in ["ceiling_raw", "ceiling_SB", "ceiling_sat"]:
            v = d[c].dropna()
            row[f"{c}_median"] = float(v.median()) if len(v) else np.nan
            row[f"{c}_min"] = float(v.min()) if len(v) else np.nan
            row[f"n_exceed_{c}"] = int((d.cos > d[c]).sum())
        S.append(row)
    SS = pd.DataFrame(S)
    SS.round(4).to_csv(OUT / "ceiling_comparison.tsv", sep="\t", index=False)

    # ------------------------------------------------ 報告
    cs = P[P["class"] == "cross_species"]
    L = ["# G: 天井推定の妥当性検査", "",
         "## 1. r(n) 曲線が取れる状態", "",
         "disjoint な2半分を作るには 2n <= min(n_case, n_ctrl) が必要。", "",
         "| 状態 | 症例 n | 対照 n | 取れる n の上限 |", "|---|---|---|---|"]
    for state, d in both.groupby("state", sort=False):
        r = d.iloc[0]
        L.append(f"| {state} | {int(r.n_case_avail)} | {int(r.n_ctrl_avail)} | {int(r.n_max)} |")
    L += ["", "**曲線（2点以上）が取れるのはネコ4状態のみ。マウス12状態は症例 n=3 のため "
          "n=1 しか取れず、飽和の検証ができません。**", "",
          "| 状態 | n | r(n) |", "|---|---|---|"]
    for _, r in both.iterrows():
        if np.isfinite(r.r):
            L.append(f"| {r.state} | {int(r.n_case)} | {r.r:.3f} |")
    L += ["", "副曲線（症例を1に固定し対照側のみ変化）:", "",
          "| 状態 | 対照 n | r |", "|---|---|---|"]
    for _, r in co.iterrows():
        if np.isfinite(r.r):
            L.append(f"| {r.state} | {int(r.n_ctrl)} | {r.r:.3f} |")
    if len(VAL):
        L += ["", "## 2. Spearman-Brown の内部検証（両側を変える曲線）", "",
              "| 状態 | k | r(1) | 実測 r(k) | SB 予言 | 実測 − 予言 |", "|---|---|---|---|---|---|"]
        for _, r in VAL.iterrows():
            L.append(f"| {r.state} | {int(r.k)} | {r.r_1:.3f} | {r.r_k_observed:.3f} "
                     f"| {r.r_k_SB_predicted:.3f} | {r.observed_minus_predicted:+.3f} |")
    if len(VAL2):
        L += ["", "対照側のみの曲線での検証（マウス状態を含む）:", "",
              "| 状態 | k | r(1) | 実測 r(k) | SB 予言 | 実測 − 予言 |", "|---|---|---|---|---|---|"]
        for _, r in VAL2.iterrows():
            L.append(f"| {r.state} | {int(r.k)} | {r.r_1:.3f} | {r.r_k_observed:.3f} "
                     f"| {r.r_k_SB_predicted:.3f} | {r.observed_minus_predicted:+.3f} |")
    L += ["", "## 3. 3通りの天井", "",
          "| クラス | n | cos 中央 | cos 最大 | 天井(raw) 中央 | 天井(SB) 中央 | 天井(sat) 中央 |"
          " | 超過数 raw / SB |", "|---|---|---|---|---|---|---|---|"]
    for _, r in SS.iterrows():
        L.append(f"| {r['class']} | {r['n']} | {r['cos_median']:.3f} | {r['cos_max']:.3f} "
                 f"| {r['ceiling_raw_median']:.3f} | {r['ceiling_SB_median']:.3f} "
                 f"| {r['ceiling_sat_median'] if np.isfinite(r['ceiling_sat_median']) else float('nan'):.3f} "
                 f"| {int(r['n_exceed_ceiling_raw'])} / {int(r['n_exceed_ceiling_SB'])} |")
    L += ["", "### 異種ペアでの判定", "",
          f"- 観測 cos 最大: **{cs.cos.max():.3f}**",
          f"- 天井(raw, 補正なし・最も保守的): 中央 {cs.ceiling_raw.median():.3f}、"
          f"最小 {cs.ceiling_raw.min():.3f}",
          f"- 天井(SB): 中央 {cs.ceiling_SB.median():.3f}、最小 {cs.ceiling_SB.min():.3f}",
          f"- 異種ペアで天井を超えたもの: raw {int((cs.cos > cs.ceiling_raw).sum())} / "
          f"SB {int((cs.cos > cs.ceiling_SB).sum())} （全 {len(cs)} 組）", "",
          "## 4. 天井を超えたペアの内訳（誤差の独立性が破れているか）", "",
          "| ペア | クラス | 同一コホート/対照共有 | cos | 天井(SB) | 天井(raw) |",
          "|---|---|---|---|---|---|"]
    for _, r in P[P.cos > P.ceiling_SB].iterrows():
        L.append(f"| {r.a} × {r.b} | {r['class']} | {'あり' if r.shared_cohort else 'なし'} "
                 f"| {r.cos:.3f} | {r.ceiling_SB:.3f} | {r.ceiling_raw:.3f} |")
    (OUT / "ceiling_validation.md").write_text("\n".join(L) + "\n")
    print()
    print("\n".join(L[L.index("## 2. Spearman-Brown の内部検証（両側を変える曲線）"):]
                    if "## 2. Spearman-Brown の内部検証（両側を変える曲線）" in L else L))


if __name__ == "__main__":
    main()
