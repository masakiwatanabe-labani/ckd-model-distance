"""F: Δ の測定精度をマッチング共変量として作る。

「量的バイアス」の本体は平均発現量ではなく Δlog2FC の測定精度である、
という定式化に対応する。2種類の精度指標を出す。

  1. moderated t の標準誤差（limma 相当の経験ベイズ縮小後）
         se_g = sqrt( s2_post_g * (1/n_case + 1/n_ctrl) )
  2. split-half による経験的な Δ のノイズ
         症例群・対照群をそれぞれ半分に割って独立に Δ を2本作り、
         その差 (Δ1 - Δ2) の標準偏差を分割を繰り返して求める。
         分布の仮定を置かない代わりに n=3 の群では粗くなる。

出力: delta_se_moderated.tsv / delta_se_splithalf.tsv
      列 = 対象の2状態（cat_CKD34, mouse_2W）、行 = ヒトシンボル
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from lib_stats import moderated_ttest  # noqa: E402

N_SPLIT = 400
SEED = 20260826


def human_projector(index, species):
    """行列 index -> ヒトシンボルへの写像を一度だけ作る（分割ごとの再計算を避ける）。"""
    m = B.MAPS[species]
    human = np.array([str(m[g]).upper() if g in m else str(g).upper() for g in index])
    seen, keep = set(), np.zeros(len(human), bool)
    for i, h in enumerate(human):
        if h not in seen:
            seen.add(h)
            keep[i] = True
    return human[keep], keep


def targets():
    cat_mat, cat_grp = B._cat_mat("ctx")
    mr = B.read_int("mouse_rna")
    mgrp = B.read_int("mouse_rna_grp").squeeze()
    gm = mr.T.groupby(mgrp).mean().T
    m = np.log2(mr[gm.max(axis=1) >= B.TH["fpkm_min"]] + 1)
    return [("cat_CKD34", cat_mat, cat_grp, "CKD3/4", ["Control"], "cat"),
            ("mouse_2W", m, mgrp, "DT100ng-2W", ["Ctrl"], "mouse")]


def main():
    rng = np.random.default_rng(SEED)
    mod, sh = {}, {}
    for state, mat, grp, case, ctrl, species in targets():
        A = mat[grp[grp == case].index]
        C = mat[grp[grp.isin(ctrl)].index]
        res = moderated_ttest(A, C)
        se = np.sqrt(res["s2_post"] * (1.0 / res["n_a"] + 1.0 / res["n_b"]))
        se = se[np.isfinite(se)]
        hv, _ = B.to_human(se, species)
        mod[state] = hv

        human, keep = human_projector(mat.index, species)
        av, cv = A.to_numpy(float), C.to_numpy(float)
        na, nc = av.shape[1], cv.shape[1]
        diffs = np.empty((N_SPLIT, keep.sum()))
        for i in range(N_SPLIT):
            pa, pc = rng.permutation(na), rng.permutation(nc)
            a1, a2 = av[:, pa[: na // 2]], av[:, pa[na // 2:]]
            c1, c2 = cv[:, pc[: nc // 2]], cv[:, pc[nc // 2:]]
            d = (np.nanmean(a1, 1) - np.nanmean(c1, 1)) - (np.nanmean(a2, 1) - np.nanmean(c2, 1))
            diffs[i] = d[keep]
        s = pd.Series(np.nanstd(diffs, axis=0, ddof=1), index=human)
        sh[state] = s[~s.index.duplicated(keep="first")]
        print(f"{state}: moderated SE 中央値 {np.nanmedian(hv):.4f} / "
              f"split-half ノイズSD 中央値 {np.nanmedian(sh[state]):.4f} "
              f"(症例 n={na} 対照 n={nc})")

    pd.DataFrame(mod).rename_axis("gene").sort_index().to_csv(
        HERE / "delta_se_moderated.tsv", sep="\t", na_rep="NA")
    pd.DataFrame(sh).rename_axis("gene").sort_index().to_csv(
        HERE / "delta_se_splithalf.tsv", sep="\t", na_rep="NA")

    # 2つの精度指標がどれくらい一致するか（どちらを使っても結論が変わらないかの目安）
    a = pd.DataFrame(mod).rank(pct=True).mean(axis=1)
    b = pd.DataFrame(sh).rank(pct=True).mean(axis=1)
    i = a.index.intersection(b.index)
    from scipy import stats as st
    print(f"2指標の順位一致: Spearman rho = {st.spearmanr(a[i], b[i])[0]:.3f} (n={len(i)})")


if __name__ == "__main__":
    main()
