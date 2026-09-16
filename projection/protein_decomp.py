"""Q: 蛋白層での分解（傍証。主解析は RNA のまま）。

猫 proteome（cortex/medulla、CKD1/2・CKD3/4）と Pod-TRECK proteome（Day14/Day21）で
RNA と同じ分解を回す。IRI に proteome が無いため 6 状態しかない。

猫プロテオームの中心化について:
    値は蛋白ごとに行中心化されている（行平均 ≈ 0）。ただし行 SD は 0.27–0.49 と
    ばらついており、一律のスケーリング（z化）は施されていない。
    中心化は行ごとの定数シフトなので Δ = mean(症例) − mean(対照) では相殺され、
    Δ ベクトルの向きも長さも影響を受けない。したがって cos と nrm は計算できる。
    旧稿の「log2FC の絶対値をマウスと直接比較してはいけない」は絶対存在量の話で、
    群間差の可否とは別。ただし寄託前に蛋白ごとのスケーリングが施されていた可能性は
    完全には否定できないため、傍証の位置づけに留める。

制約: n が少ないので順位の記述に留め、検定はしない。
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
from alpha_decomp import decompose, check_identities  # noqa: E402
from lib_stats import moderated_ttest  # noqa: E402

OUT = HERE / "results" / "protein"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260826
N_SPLIT = 500
DET = 0.5
REF = "cat_prot_CKD34"

# (状態名, 行列, 群ラベル列, 症例, 対照, 種)
SPECS = [("cat_prot_CKD12", "cat_prot_ctx", "CKD1/2", "Control", "cat"),
         ("cat_prot_CKD34", "cat_prot_ctx", "CKD3/4", "Control", "cat"),
         ("cat_med_prot_CKD12", "cat_prot_med", "CKD1/2", "Control", "cat"),
         ("cat_med_prot_CKD34", "cat_prot_med", "CKD3/4", "Control", "cat"),
         ("mouse_prot_D14", "mouse_prot", "Day14", "Ctrl", "mouse"),
         ("mouse_prot_D21", "mouse_prot", "Day21", "Ctrl", "mouse")]


def load(name):
    d = pd.read_parquet(B.INT / f"{name}.parquet")
    g = pd.read_parquet(B.INT / f"{name}_grp.parquet").squeeze()
    return d[d.notna().mean(axis=1) >= DET], g


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def main():
    rng = np.random.default_rng(SEED)
    series, mats = {}, {}
    for state, mname, case, ctrl, species in SPECS:
        mat, grp = load(mname)
        mats[state] = (mat, grp, case, ctrl, species)
        res = moderated_ttest(mat[grp[grp == case].index], mat[grp[grp == ctrl].index])
        v = res["lfc"]
        hv, _ = B.to_human(v[np.isfinite(v)], species)
        series[state] = hv
        print(f"{state:20s} case n={int((grp==case).sum())} ctrl n={int((grp==ctrl).sum())} "
              f"蛋白 {len(res)} → ヒト空間 {len(hv)}")

    D = pd.DataFrame(series)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D_all = D.dropna(axis=0, how="any")
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    print(f"\n全6状態 complete case: {len(D_all)} 蛋白 / うち Group A: {len(D)} 蛋白")
    D.round(6).to_csv(OUT / "delta_matrix_protein.tsv", sep="\t")

    # ------------------------------------------------ 分解
    states = list(D.columns)
    a = D[REF].to_numpy(float)
    P = decompose(a, D.to_numpy(float))
    P.index = states
    check_identities(P)

    acc = {k: np.empty((2000, len(states))) for k in ("alpha", "cos", "nrm")}
    X = D.to_numpy(float)
    for b in range(2000):
        i = rng.integers(0, len(D), len(D))
        d = decompose(a[i], X[i])
        for k in acc:
            acc[k][b] = d[k].to_numpy()
    for k, M in acc.items():
        P[f"{k}_lo95"] = np.percentile(M, 2.5, axis=0)
        P[f"{k}_hi95"] = np.percentile(M, 97.5, axis=0)
    P.insert(0, "n_proteins", len(D))

    # ------------------------------------------------ 天井
    rel = {}
    for state, (mat, grp, case, ctrl, species) in mats.items():
        cs = list(grp[grp == case].index)
        ks = list(grp[grp == ctrl].index)
        vals = []
        for _ in range(N_SPLIT):
            c, k = rng.permutation(cs), rng.permutation(ks)
            h, hk = len(c) // 2, len(k) // 2
            d1 = mat[list(c[:h])].mean(axis=1) - mat[list(k[:hk])].mean(axis=1)
            d2 = mat[list(c[h:2*h])].mean(axis=1) - mat[list(k[hk:2*hk])].mean(axis=1)
            h1, _ = B.to_human(d1[np.isfinite(d1)], species)
            h2, _ = B.to_human(d2[np.isfinite(d2)], species)
            u, v = h1.reindex(D.index), h2.reindex(D.index)
            m = u.notna() & v.notna()
            if m.sum() > 50:
                vals.append(cosine(u[m].to_numpy(), v[m].to_numpy()))
        r = np.median(vals)
        rel[state] = float(np.clip(2 * r / (1 + r), 0, 1)) if r > 0 else 0.0
        P.loc[state, "r_half"] = r
        P.loc[state, "reliability_SB"] = rel[state]
    P["ceiling_vs_ref"] = [float(np.sqrt(rel[REF] * rel[s])) for s in states]
    P.round(4).to_csv(OUT / "projection_protein.tsv", sep="\t")

    print("\n=== 蛋白層の分解（基準軸 = 猫皮質 CKD3/4 の protein Δ） ===")
    print(P[["n_proteins", "alpha", "cos", "nrm", "R_perp",
             "reliability_SB", "ceiling_vs_ref"]].round(3).to_string())

    # ------------------------------------------------ RNA との比較
    R = pd.read_csv(HERE / "results" / "alpha_groupA" / "projection.tsv", sep="\t", index_col=0)
    pairs = [("cat_prot_CKD12", "cat_CKD12"), ("cat_prot_CKD34", "cat_CKD34"),
             ("cat_med_prot_CKD12", "cat_med_CKD12"), ("cat_med_prot_CKD34", "cat_med_CKD34"),
             ("mouse_prot_D14", "mouse_2W"), ("mouse_prot_D21", "mouse_3W")]
    cmp = pd.DataFrame([{"state_protein": p, "state_rna": r,
                         "cos_protein": P.loc[p, "cos"], "cos_rna": R.loc[r, "cos"],
                         "alpha_protein": P.loc[p, "alpha"], "alpha_rna": R.loc[r, "alpha"],
                         "nrm_protein": P.loc[p, "nrm"], "nrm_rna": R.loc[r, "nrm"]}
                        for p, r in pairs])
    cmp.round(4).to_csv(OUT / "protein_vs_rna.tsv", sep="\t", index=False)
    print("\n=== RNA との対応 ===")
    print(cmp.round(3).to_string(index=False))
    nr = cmp.dropna()
    print(f"\n  cos の順位一致（6状態、基準軸自身を含む）: Spearman "
          f"{stats.spearmanr(nr.cos_protein, nr.cos_rna)[0]:+.3f}")
    m5 = nr[nr.state_protein != REF]
    print(f"  基準軸を除く5状態: Spearman {stats.spearmanr(m5.cos_protein, m5.cos_rna)[0]:+.3f}")
    print(f"  cos: protein 中央 {nr.cos_protein.median():.3f} / RNA 中央 {nr.cos_rna.median():.3f}")

    # ------------------------------------------------ 旧稿 rho の S1 再計算
    print("\n=== 旧稿 protein x protein rho の S1 条件での再計算 ===")
    x, y = series["cat_prot_CKD34"], series["mouse_prot_D21"]
    for lab, idx in [("全共通蛋白（ortholog 写像後）", x.index.intersection(y.index)),
                     ("Group A に限定", D.index)]:
        u, v = x.reindex(idx), y.reindex(idx)
        m = u.notna() & v.notna()
        print(f"  {lab:26s} n={int(m.sum()):5d}  Spearman rho={stats.spearmanr(u[m], v[m])[0]:+.4f}  "
              f"cos={cosine(u[m].to_numpy(), v[m].to_numpy()):+.4f}")
    print("  （旧稿: cat_prot_ctx_late x m_prot_d21, ortholog, n=2184, rho=0.5996）")
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
