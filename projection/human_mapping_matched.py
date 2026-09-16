"""§2.5 再検定: オルソログ写像品質を遺伝子ごとに揃えた設計。

問題:
    cat x human はネコ側だけが、mouse x human はマウス側だけが写像を経る。
    写像品質は非対称（cat→human perc_id 中央 89.0、mouse→human 92.4）で、
    層別すると種差が写像品質の関数として動き、両種とも良い層では逆転した。

設計:
    層別ではなく **遺伝子ごとの caliper マッチング**を使う。
    |perc_id(cat→human) − perc_id(mouse→human)| <= delta を満たす遺伝子だけを残す。
    こうすると集合内の全遺伝子について両種の写像品質が同程度になり、
    「高cat/低mouse」と「低cat/高mouse」が層内で相殺される問題が起きない。

判定:
    マッチ後も cat x human > mouse x human が残れば §2.5 は成立し、
    KPMP snRNA-seq を独立検証コホートに使う段階へ進む。
    消えるなら、この主張は現行データでは立たない。

対照:
    遺伝子集合を絞ると検出力が落ちるだけかもしれないので、
    同じ大きさのランダム部分集合（同一プールから抽出）の AUC 分布を必ず併記する。
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import human_extrapolation as K  # noqa: E402
import human_checks as HC  # noqa: E402
from lib_stats import moderated_ttest  # noqa: E402

OUT = HERE / "results" / "human_exploratory"
SEED = 20260826
CALIPERS = [2.0, 5.0, 10.0]


def main():
    rng = np.random.default_rng(SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    plats = K.load_ercb()

    STRONG = ["hs_ERCB_allCKD_A", "hs_ERCB_DN_A", "hs_ERCB_RPGN_A"]
    human = {}
    for name, stem, case, ctrl in K.ERCB_STATES:
        if name not in STRONG:
            continue
        mat, grp = plats[stem]
        v = moderated_ttest(mat[grp[grp.isin(case)].index],
                            mat[grp[grp.isin(ctrl)].index])["lfc"]
        human[name] = v[np.isfinite(v)]

    G = pd.Index(sorted(set(D.index) & set(plats["GSE104954-GPL22945"][0].index)))
    H = pd.DataFrame({k: v.reindex(G) for k, v in human.items()})
    G = G[D.reindex(G).notna().all(axis=1) & H.notna().all(axis=1)]
    Dh, H = D.reindex(G), H.reindex(G)

    pid = HC.perc_id_tables(G)
    pool = G.intersection(pid.dropna().index)          # 両種の perc_id が引ける遺伝子
    pidp = pid.loc[pool]
    print(f"全共通遺伝子 {len(G)} / 両種の perc_id が引ける {len(pool)}")

    DV = Dh.loc[pool].to_numpy(float)
    HV = H.loc[pool].to_numpy(float)
    is_cat = np.array([c.startswith("cat") for c in Dh.columns])
    states = list(Dh.columns)
    n_cat = int(is_cat.sum())

    def cos_matrix(pos):
        d, h = DV[pos], HV[pos]
        dn = d / np.linalg.norm(d, axis=0, keepdims=True)
        hn = h / np.linalg.norm(h, axis=0, keepdims=True)
        return dn.T @ hn, dn

    def auc(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return float(((x[:, None] > y[None, :]).sum()
                      + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size))

    def obs_auc(pos):
        C, _ = cos_matrix(pos)
        return auc(C[is_cat].ravel(), C[~is_cat].ravel())

    def exact_p(pos):
        """種ラベルの状態レベル網羅並べ替え C(16,4)=1820。"""
        C, _ = cos_matrix(pos)
        obs = auc(C[is_cat].ravel(), C[~is_cat].ravel())
        null = []
        for c in combinations(range(len(states)), n_cat):
            m = np.zeros(len(states), bool)
            m[list(c)] = True
            null.append(auc(C[m].ravel(), C[~m].ravel()))
        null = np.array(null)
        return obs, float((null >= obs).sum() / null.size), float(np.median(null))

    def summary(pos, label, extra=None):
        C, dn = cos_matrix(pos)
        cm = dn.T @ dn
        obs, p, nullmed = exact_p(pos)
        # 同サイズのランダム部分集合（同一プールから）
        rnd = np.array([obs_auc(rng.choice(len(pool), size=len(pos), replace=False))
                        for _ in range(500)])
        row = {"subset": label, "n_genes": len(pos),
               "pid_cat_median": float(pidp.pid_cat.to_numpy()[pos].mean()),
               "pid_mouse_median": float(pidp.pid_mouse.to_numpy()[pos].mean()),
               "cat_x_human": float(np.median(C[is_cat])),
               "mouse_x_human": float(np.median(C[~is_cat])),
               "diff": float(np.median(C[is_cat]) - np.median(C[~is_cat])),
               "cat_x_mouse": float(np.median(cm[np.ix_(is_cat, ~is_cat)])),
               "auc": obs, "p_exact": p, "null_median": nullmed,
               "size_null_median": float(np.median(rnd)),
               "size_null_p5": float(np.percentile(rnd, 5)),
               "size_null_p95": float(np.percentile(rnd, 95))}
        row.update(extra or {})
        return row

    rows = [summary(np.arange(len(pool)), "全遺伝子（マッチなし）")]
    dpid = (pidp.pid_cat - pidp.pid_mouse).to_numpy()
    for delta in CALIPERS:
        pos = np.where(np.abs(dpid) <= delta)[0]
        if len(pos) < 200:
            print(f"  caliper {delta}: n={len(pos)} — 少なすぎるため評価しない")
            continue
        rows.append(summary(pos, f"caliper |Δperc_id| <= {delta:g}"))
    # マッチしたうえで、なお写像品質が高い遺伝子だけに絞る
    pos = np.where((np.abs(dpid) <= 5.0)
                   & (pidp.pid_cat.to_numpy() >= np.median(pidp.pid_cat)))[0]
    if len(pos) >= 200:
        rows.append(summary(pos, "caliper 5 かつ perc_id 上位半分"))

    R = pd.DataFrame(rows)
    R.round(4).to_csv(OUT / "mapping_matched_test.tsv", sep="\t", index=False)

    print()
    cols = ["subset", "n_genes", "pid_cat_median", "pid_mouse_median",
            "cat_x_human", "mouse_x_human", "diff", "cat_x_mouse"]
    print(R[cols].round(3).to_string(index=False))
    print()
    print(R[["subset", "n_genes", "auc", "p_exact",
             "size_null_median", "size_null_p5", "size_null_p95"]].round(4).to_string(index=False))
    print()
    for _, r in R.iterrows():
        if r.subset.startswith("全"):
            continue
        inside = r.size_null_p5 <= r.auc <= r.size_null_p95
        print(f"  {r.subset}: AUC {r.auc:.3f} は同サイズ乱択の 90% 区間 "
              f"[{r.size_null_p5:.3f}, {r.size_null_p95:.3f}] の"
              f"{'内側 → サイズ低下では説明できない差はない' if inside else '外側'}")
    print(f"\n書き出し: {OUT / 'mapping_matched_test.tsv'}")


if __name__ == "__main__":
    main()
