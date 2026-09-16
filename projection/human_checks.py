"""K の本採用に伴う検査。§2.3 と同じ一巡をヒト軸に対して行う。

この解析に固有の交絡が2つある。

  (i)  写像の非対称性。cat x human はネコ側だけが、mouse x human はマウス側だけが
       オルソログ写像を経る。mouse-human のオルソロジーは cat-human より整備されて
       いるので、写像品質が効くならマウス側が有利になるはず（＝観測と逆向き）。
       これを実測で確認し、写像品質を揃えても差が残るかを見る。
  (ii) 状態の信頼性。信頼性が高い状態ほどヒトに近く見えるなら、それは種ではなく
       測定の話になる。マウス側のほうが信頼性が高いので、これも観測と逆向きのはず。

平均発現量は、両クラスが同一の 1,690 遺伝子で計算されるため、クラス間で差が
つきようがない（遺伝子集合が共通）。念のためヒト側対照の発現量で層別もする。

共有対照については、ネコ／マウスの状態とヒトの状態は検体を一切共有しないため
cross-species ペアには効かない。ヒト3状態どうしは同じ LD 18 検体を共有し、さらに
allCKD_A は DN_A / RPGN_A を包含するので、ヒト内 cos は上振れしている。
ヒト内の値は参照値としてのみ扱い、入れ子でないヒト1状態だけを使った感度分析も出す。
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
import human_extrapolation as K  # noqa: E402

OUT = HERE / "results" / "human_exploratory"


def auc(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return float(((x[:, None] > y[None, :]).sum()
                  + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size))


def perm(hp, drop=(), human_keep=None):
    if len(hp) == 0 or "nh" not in hp.columns:
        return {"n_states": 0, "n_pairs": 0, "auc": np.nan,
                "null_median": np.nan, "n_perm": 0, "p_exact": np.nan}
    d = hp[~hp.nh.isin(drop)]
    if human_keep is not None:
        d = d[d.hs.isin(human_keep)]
    st = sorted(d.nh.unique())
    nc = sum(x.startswith("cat") for x in st)

    def stat(cs):
        m = d.nh.isin(cs)
        return auc(d.cos[m], d.cos[~m]) if m.sum() >= 3 and (~m).sum() >= 3 else np.nan

    obs = stat([x for x in st if x.startswith("cat")])
    null = np.array([v for v in (stat([st[i] for i in c])
                                 for c in combinations(range(len(st)), nc))
                     if np.isfinite(v)])
    return {"n_states": len(st), "n_pairs": len(d), "auc": obs,
            "null_median": float(np.median(null)), "n_perm": int(null.size),
            "p_exact": float((null >= obs).sum() / null.size)}


def perc_id_tables(genes):
    def load(fname, key):
        o = pd.read_csv(B.REF / fname, sep="\t")
        o = o[o.orthology_type.isin(B.cfg["orthology"]["keep_types"])]
        o = o.dropna(subset=["src_symbol", "tgt_symbol"]).copy()
        o["rank"] = (-o.confidence.fillna(0)) * 1000 - o.perc_id.fillna(0)
        o = o.sort_values("rank").drop_duplicates("src_symbol")
        o["h"] = o.tgt_symbol.astype(str).str.upper()
        return o.drop_duplicates("h").set_index("h")["perc_id"].rename(key)
    return pd.concat([load("orthologs_cat2human.tsv", "pid_cat"),
                      load("orthologs_mouse2human.tsv", "pid_mouse")],
                     axis=1).reindex(genes)


def main():
    rng = np.random.default_rng(K.SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    plats = K.load_ercb()

    STRONG = ["hs_ERCB_allCKD_A", "hs_ERCB_DN_A", "hs_ERCB_RPGN_A"]
    human, hmat = {}, {}
    for name, stem, case, ctrl in K.ERCB_STATES:
        if name not in STRONG:
            continue
        mat, grp = plats[stem]
        from lib_stats import moderated_ttest
        A = mat[grp[grp.isin(case)].index]
        C = mat[grp[grp.isin(ctrl)].index]
        v = moderated_ttest(A, C)["lfc"]
        human[name] = v[np.isfinite(v)]
        hmat[name] = (mat, grp, case, ctrl)

    G = pd.Index(sorted(set(D.index) & set(plats["GSE104954-GPL22945"][0].index)))
    H = pd.DataFrame({k: v.reindex(G) for k, v in human.items()})
    G = G[D.reindex(G).notna().all(axis=1) & H.notna().all(axis=1)]
    Dh, H = D.reindex(G), H.reindex(G)
    print(f"共通遺伝子 {len(G)}")

    def pairs(gene_idx):
        rows = []
        for nh in Dh.columns:
            for hs in H.columns:
                u, v = Dh.loc[gene_idx, nh], H.loc[gene_idx, hs]
                m = u.notna() & v.notna()
                if m.sum() < 200:
                    continue
                rows.append({"nh": nh, "hs": hs,
                             "class": "cat x human" if nh.startswith("cat") else "mouse x human",
                             "cos": K.cosine(u[m].to_numpy(), v[m].to_numpy())})
        return pd.DataFrame(rows)

    base = pairs(G)
    print("\n=== 基準（全 %d 遺伝子） ===" % len(G))
    print(pd.DataFrame([perm(base)]).round(4).to_string(index=False))

    # ---------------------------------------- (i) 写像の非対称性
    pid = perc_id_tables(G)
    ok = pid.dropna()
    print("\n=== (i) オルソログ写像の非対称性 ===")
    print(f"  cat→human   perc_id: 中央 {ok.pid_cat.median():.1f} "
          f"[Q1 {ok.pid_cat.quantile(.25):.1f}, Q3 {ok.pid_cat.quantile(.75):.1f}]")
    print(f"  mouse→human perc_id: 中央 {ok.pid_mouse.median():.1f} "
          f"[Q1 {ok.pid_mouse.quantile(.25):.1f}, Q3 {ok.pid_mouse.quantile(.75):.1f}]")
    print(f"  mouse のほうが高い遺伝子の割合: {(ok.pid_mouse > ok.pid_cat).mean():.1%}")

    # 遺伝子集合を絞ると AUC が下がるのは、写像品質のせいか単にサイズのせいか。
    # 同じサイズのランダム部分集合を引いて帰無分布を作り、両者を区別する。
    # 観測 AUC だけを計算する軽い版（並べ替えは回さない）。
    # cos は遺伝子集合を絞るだけで再計算できるので、行列演算で一括処理する。
    DV = Dh.to_numpy(float)
    HV = H.to_numpy(float)
    is_cat = np.array([c.startswith("cat") for c in Dh.columns])

    def obs_auc(pos):
        d = DV[pos]
        h = HV[pos]
        dn = d / np.linalg.norm(d, axis=0, keepdims=True)
        hn = h / np.linalg.norm(h, axis=0, keepdims=True)
        C = dn.T @ hn                       # (非ヒト状態 x ヒト状態) の cos
        x, y = C[is_cat].ravel(), C[~is_cat].ravel()
        return auc(x, y)

    def size_matched_null(k, n_draw=500):
        return np.array([obs_auc(rng.choice(len(G), size=k, replace=False))
                         for _ in range(n_draw)])

    rows = []
    subsets = []
    for q, lab in [(0.50, "両種とも perc_id 上位50%"), (0.25, "両種とも perc_id 上位25%")]:
        thr_c, thr_m = ok.pid_cat.quantile(1 - q), ok.pid_mouse.quantile(1 - q)
        subsets.append((lab, G.intersection(ok.index[(ok.pid_cat >= thr_c)
                                                     & (ok.pid_mouse >= thr_m)])))
    subsets.append(("両種とも perc_id 下位50%",
                    G.intersection(ok.index[(ok.pid_cat <= ok.pid_cat.median())
                                            & (ok.pid_mouse <= ok.pid_mouse.median())])))
    subsets.append(("cat の perc_id >= mouse", G.intersection(ok.index[ok.pid_cat >= ok.pid_mouse])))
    subsets.append(("cat の perc_id <  mouse", G.intersection(ok.index[ok.pid_cat < ok.pid_mouse])))

    for lab, gi in subsets:
        r = perm(pairs(gi))
        if not np.isfinite(r["auc"]):
            print(f"  {lab}: n={len(gi)} — 遺伝子が少なくペアが作れず評価不能")
            rows.append({"subset": lab, "n_genes": len(gi), **r,
                         "size_null_median": np.nan, "size_null_p5": np.nan,
                         "explained_by_size": None})
            continue
        nullsz = size_matched_null(len(gi))  # 同サイズのランダム遺伝子集合
        p5 = float(np.percentile(nullsz, 5)) if nullsz.size else np.nan
        rows.append({"subset": lab, "n_genes": len(gi), **r,
                     "size_null_median": float(np.median(nullsz)) if nullsz.size else np.nan,
                     "size_null_p5": p5,
                     "explained_by_size": bool(r["auc"] >= p5) if np.isfinite(p5) else None})
        print(f"  {lab}: n={len(gi)}, AUC {r['auc']:.4f}, p={r['p_exact']:.4f} | "
              f"同サイズ乱択 AUC 中央 {np.median(nullsz):.4f} "
              f"[5%点 {p5:.4f}] → サイズで説明"
              f"{'できる' if r['auc'] >= p5 else 'できない'}")

    # ---------------------------------------- (ii) 状態の信頼性
    rel = pd.read_csv(HERE / "results" / "reliability" / "reliability.tsv",
                      sep="\t").set_index("state")["reliability_SB_median"]
    per = base.groupby("nh").cos.median()
    df = pd.DataFrame({"cos_to_human": per, "reliability": rel.reindex(per.index),
                       "species": ["cat" if s.startswith("cat") else "mouse" for s in per.index]})
    from scipy import stats as st
    rho = st.spearmanr(df.reliability, df.cos_to_human)[0]
    print("\n=== (ii) 状態の信頼性は説明するか ===")
    print(f"  信頼性 x 対ヒト cos の順位相関: rho = {rho:+.3f} (n={len(df)})")
    print(f"  猫4状態の信頼性 中央 {df[df.species=='cat'].reliability.median():.3f} / "
          f"マウス12状態 {df[df.species=='mouse'].reliability.median():.3f}")
    df.round(4).to_csv(OUT / "state_reliability_vs_human.tsv", sep="\t")

    # ---------------------------------------- (iii) ヒト状態の入れ子・共有対照
    print("\n=== (iii) ヒト状態の入れ子と共有対照 ===")
    hs_rows = []
    for hs in STRONG:
        r = perm(base, human_keep=[hs])
        hs_rows.append({"human_state": hs, **r})
        print(f"  {hs:20s} のみ: AUC {r['auc']:.4f}, p={r['p_exact']:.4f} "
              f"({r['n_pairs']} ペア)")
    pd.DataFrame(hs_rows).round(4).to_csv(OUT / "per_human_state_test.tsv", sep="\t", index=False)

    # ヒト内ペアの共有対照バイアス（分割対照）
    from lib_stats import moderated_ttest  # noqa: F811
    mat, grp = plats["GSE104954-GPL22945"]
    ctrl_cols = list(grp[grp == "control_LD"].index)
    sub = []
    for a, b in combinations(STRONG, 2):
        _, _, ca, _ = hmat[a]
        _, _, cb, _ = hmat[b]
        nested = set(ca) < set(cb) or set(cb) < set(ca)
        vals = []
        for _ in range(200):
            k = rng.permutation(ctrl_cols)
            h = len(k) // 2
            d1 = mat[grp[grp.isin(ca)].index].mean(axis=1) - mat[list(k[:h])].mean(axis=1)
            d2 = mat[grp[grp.isin(cb)].index].mean(axis=1) - mat[list(k[h:2 * h])].mean(axis=1)
            d1, d2 = d1.reindex(G), d2.reindex(G)
            m = d1.notna() & d2.notna()
            vals.append(K.cosine(d1[m].to_numpy(), d2[m].to_numpy()))
        u = H[a]; v = H[b]; m = u.notna() & v.notna()
        sub.append({"pair": f"{a} x {b}", "case_nested": nested,
                    "cos_shared_control": K.cosine(u[m].to_numpy(), v[m].to_numpy()),
                    "cos_split_control": float(np.median(vals))})
    SUB = pd.DataFrame(sub)
    SUB.round(4).to_csv(OUT / "human_within_shared_control.tsv", sep="\t", index=False)
    print(SUB.round(3).to_string(index=False))

    pd.DataFrame(rows).round(4).to_csv(OUT / "human_mapping_checks.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
