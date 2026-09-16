"""D: split-half 信頼性と減衰の天井。

各状態について、症例群と対照群をそれぞれランダムに半分に割り、
半分ごとに独立に Δ を作って半分同士の cos を取る。Spearman-Brown 補正

    r = 2 r_half / (1 + r_half)

で全検体版 Δ の信頼性を推定する。ペア (X, Y) の観測 cos の天井は

    ceiling = sqrt(r_XX * r_YY)

減衰補正 cos = cos_obs / ceiling。

判定:
    異種ペアの天井が観測最大 0.548 に近ければ、「マウスは 0.55 を超えない」は
    測定限界の記述であって種差の主張にならない。天井が十分高ければ主張は立つ。

注意（結果の読み方を左右する）:
    Pod-TRECK は Ctrl n=3・各症例群 n=3、IRI も症例 n=3。半分に割ると 1 対 2 に
    なり、半分ごとの Δ は極めてノイジーで r_half は下振れする。Spearman-Brown は
    それを引き上げるが、n=3 の分割から得た信頼性は不確かさが大きい。
    ネコ（対照 6、症例 5〜8）はまだしも安定する。
    したがって「天井」は点推定ではなく CI で読むこと。
"""
from __future__ import annotations

import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

# 遺伝子空間の差し替え（AA. 感度解析）。既定は Group A intersection で、
# 環境変数を与えない限り従来と完全に同じ動作になる。
TAG = os.environ.get("XSP_TAG", "")
GENES = Path(os.environ.get("XSP_GENES", str(HERE / "groupA_intersection.txt")))
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
import cos_matrix as CM  # noqa: E402

OUT = HERE / "results" / ("reliability" + TAG)
OUT.mkdir(parents=True, exist_ok=True)
N_SPLIT = 500
SEED = 20260826


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def state_specs():
    """(状態名, 発現行列, 群ラベル, 症例ラベル, 対照ラベル群, 種) を返す。"""
    out = []
    for state, label, tissue in B.CAT_STATES:
        mat, grp = B._cat_mat(tissue)
        out.append((state, mat, grp, label, ["Control"], "cat"))

    mr = B.read_int("mouse_rna")
    mgrp = B.read_int("mouse_rna_grp").squeeze()
    gm = mr.T.groupby(mgrp).mean().T
    m = np.log2(mr[gm.max(axis=1) >= B.TH["fpkm_min"]] + 1)
    for state, label, _h in B.PODTRECK_STATES:
        out.append((state, m, mgrp, label, ["Ctrl"], "mouse"))

    imat = B.read_int("mouse_iri_matrix")
    igrp = B.read_int("mouse_iri_grp").squeeze()
    il = np.log2(imat + 1)
    for state, label, _h, ctrl in B.IRI_STATES:
        out.append((state, il, igrp, label, ctrl, "mouse"))
    return out


def split_half_reliability(mat, grp, case_label, ctrl_labels, genes, species, rng):
    """半分同士の cos と Spearman-Brown 補正後の信頼性の分布。"""
    case = list(grp[grp == case_label].index)
    ctrl = list(grp[grp.isin(ctrl_labels)].index)
    means = {}

    def dmean(cols):
        key = tuple(sorted(cols))
        if key not in means:
            means[key] = mat[list(cols)].mean(axis=1)
        return means[key]

    r_half = np.empty(N_SPLIT)
    for i in range(N_SPLIT):
        c = rng.permutation(case)
        k = rng.permutation(ctrl)
        c1, c2 = c[: len(c) // 2], c[len(c) // 2:]
        k1, k2 = k[: len(k) // 2], k[len(k) // 2:]
        d1 = dmean(c1) - dmean(k1)
        d2 = dmean(c2) - dmean(k2)
        h1, _ = B.to_human(d1[np.isfinite(d1)], species)
        h2, _ = B.to_human(d2[np.isfinite(d2)], species)
        u, v = h1.reindex(genes), h2.reindex(genes)
        m = u.notna() & v.notna()
        r_half[i] = cosine(u[m].to_numpy(), v[m].to_numpy()) if m.sum() > 50 else np.nan

    rh = r_half[np.isfinite(r_half)]
    with np.errstate(invalid="ignore", divide="ignore"):
        sb = 2 * rh / (1 + rh)
    sb = np.clip(np.where(rh > 0, sb, 0.0), 0.0, 1.0)
    return rh, sb, {"n_case": len(case), "n_ctrl": len(ctrl),
                    "half_case": len(case) // 2, "half_ctrl": len(ctrl) // 2,
                    "frac_r_half_nonpositive": float((rh <= 0).mean())}


def orthology_strata(genes):
    """ヒト遺伝子ごとに、猫側とマウス側のオルソログ写像の品質を引く。"""
    def load(fname, key):
        o = pd.read_csv(B.REF / fname, sep="\t")
        o = o[o.orthology_type.isin(B.cfg["orthology"]["keep_types"])]
        o = o.dropna(subset=["src_symbol", "tgt_symbol"]).copy()
        o["rank"] = (-o.confidence.fillna(0)) * 1000 - o.perc_id.fillna(0)
        o = o.sort_values("rank").drop_duplicates("src_symbol")
        o["h"] = o.tgt_symbol.astype(str).str.upper()
        o = o.drop_duplicates("h")
        return o.set_index("h")[["confidence", "perc_id"]].rename(
            columns={"confidence": f"conf_{key}", "perc_id": f"pid_{key}"})

    df = load("orthologs_cat2human.tsv", "cat").join(
        load("orthologs_mouse2human.tsv", "mouse"), how="outer")
    return df.reindex(genes)


def main():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    genes = D.index
    rng = np.random.default_rng(SEED)

    # ---------------------------------------------------------- 1. 信頼性
    rows, sb_draws = [], {}
    for state, mat, grp, label, ctrl, species in state_specs():
        rh, sb, meta = split_half_reliability(mat, grp, label, ctrl, genes, species, rng)
        sb_draws[state] = sb
        rows.append({"state": state, "species": species, **meta,
                     "r_half_median": float(np.median(rh)),
                     "reliability_SB_median": float(np.median(sb)),
                     "reliability_SB_lo95": float(np.percentile(sb, 2.5)),
                     "reliability_SB_hi95": float(np.percentile(sb, 97.5)),
                     "n_splits": len(rh)})
        print(f"{state:14s} n_case={meta['n_case']:2d}/{meta['half_case']} "
              f"n_ctrl={meta['n_ctrl']:2d}/{meta['half_ctrl']} "
              f"r_half={np.median(rh):+.3f} -> r_SB={np.median(sb):.3f} "
              f"[{np.percentile(sb, 2.5):.3f}, {np.percentile(sb, 97.5):.3f}]")
    REL = pd.DataFrame(rows)
    REL.round(4).to_csv(OUT / "reliability.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 2-3. 天井と補正
    P = pd.read_csv(HERE / "results" / ("control_analysis" + TAG) / "pairs.tsv", sep="\t")
    out = []
    for _, r in P.iterrows():
        a, b = r["a"], r["b"]
        sa, sb_ = sb_draws[a], sb_draws[b]
        ceil = np.sqrt(sa * sb_)
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = np.where(ceil > 0, r["cos"] / ceil, np.nan)
        corr = corr[np.isfinite(corr)]
        out.append({**r.to_dict(),
                    "ceiling_median": float(np.median(ceil)),
                    "ceiling_lo95": float(np.percentile(ceil, 2.5)),
                    "ceiling_hi95": float(np.percentile(ceil, 97.5)),
                    "cos_corrected": float(np.median(corr)) if corr.size else np.nan,
                    "cos_corrected_lo95": float(np.percentile(corr, 2.5)) if corr.size else np.nan,
                    "cos_corrected_hi95": float(np.percentile(corr, 97.5)) if corr.size else np.nan,
                    "cos_exceeds_ceiling": bool(r["cos"] > np.median(ceil))})
    PC = pd.DataFrame(out)
    PC.round(4).to_csv(OUT / "pairs_corrected.tsv", sep="\t", index=False)

    def summ(d, lab, col):
        v = d[col].to_numpy(float)
        v = v[np.isfinite(v)]
        return {"set": lab, "metric": col, "n": len(v), "min": v.min(),
                "median": float(np.median(v)), "max": v.max()}

    S = []
    for cls in ["within_dataset", "same_species_diff_dataset", "cross_species"]:
        d = PC[PC["class"] == cls]
        for col in ["cos", "ceiling_median", "cos_corrected"]:
            S.append(summ(d, cls, col))
    SS = pd.DataFrame(S)
    SS.round(4).to_csv(OUT / "corrected_distributions.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 4. 写像品質の層別
    strat = orthology_strata(genes)
    both = strat.dropna(subset=["pid_cat", "pid_mouse"])
    pid = both[["pid_cat", "pid_mouse"]].min(axis=1)
    q = pd.qcut(pid.rank(method="first"), 4, labels=[1, 2, 3, 4])
    cs = PC[PC["class"] == "cross_species"]
    ss = PC[PC["class"] == "same_species_diff_dataset"]
    X = D
    strat_rows = []
    for qq in [1, 2, 3, 4]:
        g = both.index[(q == qq).to_numpy()]
        g = X.index.intersection(g)
        for lab, dd in [("cross_species", cs), ("same_species_diff_dataset", ss)]:
            vals = [cosine(X.loc[g, r["a"]].to_numpy(), X.loc[g, r["b"]].to_numpy())
                    for _, r in dd.iterrows()]
            strat_rows.append({"pid_quartile": qq, "n_genes": len(g), "class": lab,
                               "median_cos": float(np.median(vals)),
                               "min_pid": float(pid[q == qq].min()),
                               "max_pid": float(pid[q == qq].max())})
    ST = pd.DataFrame(strat_rows)
    ST.round(4).to_csv(OUT / "orthology_strata.tsv", sep="\t", index=False)

    # 高信頼度 1:1 のみ
    hi = both.index[(both.conf_cat == 1) & (both.conf_mouse == 1)]
    hi = X.index.intersection(hi)
    hi_rows = []
    for lab, dd in [("cross_species", cs), ("same_species_diff_dataset", ss)]:
        vals = [cosine(X.loc[hi, r["a"]].to_numpy(), X.loc[hi, r["b"]].to_numpy())
                for _, r in dd.iterrows()]
        hi_rows.append({"set": "high_confidence_1to1", "n_genes": len(hi), "class": lab,
                        "median_cos": float(np.median(vals)), "max_cos": float(np.max(vals))})
    HI = pd.DataFrame(hi_rows)
    HI.round(4).to_csv(OUT / "orthology_high_confidence.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- 報告
    cs_ceil = PC[PC["class"] == "cross_species"]["ceiling_median"]
    cs_obs_max = PC[PC["class"] == "cross_species"]["cos"].max()
    L = ["# D: split-half 信頼性と減衰の天井", "",
         f"Group A intersection {len(genes)} 遺伝子、ランダム分割 {N_SPLIT} 回。", "",
         "## 1. 各状態の信頼性", "",
         "| 状態 | 症例 n (半分) | 対照 n (半分) | r_half 中央値 | 信頼性 r (SB) | 95% CI |",
         "|---|---|---|---|---|---|"]
    for _, r in REL.iterrows():
        L.append(f"| {r['state']} | {r['n_case']} ({r['half_case']}) | {r['n_ctrl']} ({r['half_ctrl']}) "
                 f"| {r['r_half_median']:+.3f} | {r['reliability_SB_median']:.3f} "
                 f"| [{r['reliability_SB_lo95']:.3f}, {r['reliability_SB_hi95']:.3f}] |")
    L += ["", "## 2-3. 天井と減衰補正後の cos", "",
          "| クラス | 指標 | n | 最小 | 中央値 | 最大 |", "|---|---|---|---|---|---|"]
    for _, r in SS.iterrows():
        L.append(f"| {r['set']} | {r['metric']} | {r['n']} | {r['min']:.3f} "
                 f"| {r['median']:.3f} | {r['max']:.3f} |")
    L += ["", "## 判定材料", "",
          f"- 異種ペアの天井: 中央値 {cs_ceil.median():.3f}、"
          f"範囲 [{cs_ceil.min():.3f}, {cs_ceil.max():.3f}]",
          f"- 異種ペアの観測 cos 最大: {cs_obs_max:.3f}",
          f"- 天井が観測最大を下回るペア数: {int(PC['cos_exceeds_ceiling'].sum())} / {len(PC)}",
          "", "## 4. オルソログ写像品質の層別（異種 cos が写像由来に減衰しているか）", "",
          "| perc_id 四分位 | pid 範囲 | 遺伝子数 | クラス | cos 中央値 |",
          "|---|---|---|---|---|"]
    for _, r in ST.iterrows():
        L.append(f"| Q{r['pid_quartile']} | {r['min_pid']:.1f}–{r['max_pid']:.1f} | {r['n_genes']} "
                 f"| {r['class']} | {r['median_cos']:.3f} |")
    L += ["", "高信頼度 1:1 オルソログのみ:", ""]
    for _, r in HI.iterrows():
        L.append(f"- {r['class']}: n_genes={r['n_genes']}, cos 中央値 {r['median_cos']:.3f}, "
                 f"最大 {r['max_cos']:.3f}")
    (OUT / "report.md").write_text("\n".join(L) + "\n")
    print()
    print("\n".join(L[L.index("## 2-3. 天井と減衰補正後の cos"):]))


if __name__ == "__main__":
    main()
