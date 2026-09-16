"""K: ヒト外挿（探索的）。主解析には組み込まない。

ERCB (GSE104954, 尿細管間質マイクロアレイ) と KPMP regional proteomics を
ヒト状態として加え、cat x human / mouse x human / cat x mouse の cos を出す。

本採用の可否を決めるために、先に3点を報告する:
  1. 各ヒトデータセットで split-half の天井が推定可能か（群サイズの制約）
  2. プラットフォーム差（ERCB=マイクロアレイ、他=RNA-seq）が cos に効く程度
     ERCB は2つのアレイ版(GPL22945/GPL24120)に同じ疾患(DN)を持つので、
     種・疾患・組織を固定したままプラットフォームだけ変えたペアが作れる。
  3. Group A の3種共通空間に遺伝子が何個残るか

制約:
  - 天井が推定できないヒトペアの cos は主張に使わない。
  - KPMP は per-sample データが無く（事前計算済みコントラストのみ）、
    さらにプロテオーム層なので RNA 層の Δ 行列とは測定量が違う。
    参考値としてのみ出し、cat x human の主張には使わない。
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
import importlib.util as _iu

_sp = _iu.spec_from_file_location("ercb", ROOT / "src" / "11_ercb.py")
ercb = _iu.module_from_spec(_sp); _sp.loader.exec_module(ercb)
from lib_stats import moderated_ttest  # noqa: E402
import build_delta_matrix as B  # noqa: E402

OUT = HERE / "results" / "human_exploratory"
OUT.mkdir(parents=True, exist_ok=True)
N_SPLIT = 300
SEED = 20260826

# (状態名, プラットフォーム, 症例ラベル群, 対照ラベル群)
ERCB_STATES = [
    ("hs_ERCB_allCKD_A", "GSE104954-GPL22945",
     ["DN", "RPGN", "FSGS", "FSGS_MCD", "MCD"], ["control_LD"]),
    ("hs_ERCB_DN_A", "GSE104954-GPL22945", ["DN"], ["control_LD"]),
    ("hs_ERCB_RPGN_A", "GSE104954-GPL22945", ["RPGN"], ["control_LD"]),
    ("hs_ERCB_allCKD_B", "GSE104954-GPL24120",
     ["SLE", "IgAN", "HT", "MGN", "DN", "FSGS", "MCD", "TMD"], ["control_LD", "control_TN"]),
    ("hs_ERCB_DN_B", "GSE104954-GPL24120", ["DN"], ["control_LD", "control_TN"]),
    ("hs_ERCB_HT_B", "GSE104954-GPL24120", ["HT"], ["control_LD", "control_TN"]),
]
SPECIES = {"cat": "cat", "mouse": "mouse", "human": "human"}


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def load_ercb():
    e2s = ercb.entrez_to_symbol()
    plats = {}
    for stem in ["GSE104954-GPL22945", "GSE104954-GPL24120"]:
        mat, grp = ercb.load_platform(stem)
        sym = [e2s.get(str(i).replace("_at", ""), f"ENTREZ:{str(i).replace('_at','')}")
               for i in mat.index]
        mat = mat.set_axis(sym, axis=0).groupby(level=0).mean()
        mat.index = [str(i).upper() for i in mat.index]
        plats[stem] = (mat[~mat.index.duplicated()], grp)
    return plats


def kpmp_lfc():
    p = ROOT / "data" / "external" / "KPMP" / "DataLake_DEPs.txt"
    if not p.exists():
        return {}
    df = pd.read_csv(p, sep="\t")
    out = {}
    for comp, name in [("CKD.vs.HRT.in.TI", "hs_KPMP_TI"), ("CKD.vs.HRT.in.G", "hs_KPMP_G")]:
        d = df[df.Comparison == comp].copy()
        gn = d.Gene_name.astype(str).str.strip()
        d = d[~(d.Gene_name.isna() | (gn == "") | (gn.str.lower() == "nan") | d.LogFC.isna())]
        d = d.assign(_sym=d.Gene_name.astype(str).str.strip().str.upper(), _abs=d.LogFC.abs())
        d = d.sort_values(["_sym", "Adj_pvalue", "_abs"], ascending=[True, True, False])
        out[name] = d.drop_duplicates("_sym", keep="first").set_index("_sym")["LogFC"]
    return out


def split_half_reliability(mat, grp, case_labels, ctrl_labels, genes, rng):
    case = list(grp[grp.isin(case_labels)].index)
    ctrl = list(grp[grp.isin(ctrl_labels)].index)
    nmax = min(len(case), len(ctrl)) // 2
    if nmax < 1:
        return None, {"n_case": len(case), "n_ctrl": len(ctrl), "n_max": nmax}
    M = mat.reindex(genes)
    vals = np.empty(N_SPLIT)
    for i in range(N_SPLIT):
        c = rng.permutation(case); k = rng.permutation(ctrl)
        h = len(c) // 2; hk = len(k) // 2
        d1 = M[list(c[:h])].mean(axis=1) - M[list(k[:hk])].mean(axis=1)
        d2 = M[list(c[h:2 * h])].mean(axis=1) - M[list(k[hk:2 * hk])].mean(axis=1)
        m = d1.notna() & d2.notna()
        vals[i] = cosine(d1[m].to_numpy(), d2[m].to_numpy()) if m.sum() > 50 else np.nan
    rh = vals[np.isfinite(vals)]
    sb = np.clip(np.where(rh > 0, 2 * rh / (1 + rh), 0.0), 0, 1)
    return (rh, sb), {"n_case": len(case), "n_ctrl": len(ctrl), "n_max": nmax,
                      "half_case": len(case) // 2, "half_ctrl": len(ctrl) // 2}


def main():
    rng = np.random.default_rng(SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    plats = load_ercb()

    # ---------------------------------------------- ヒト状態の Δ
    human, meta = {}, []
    for name, stem, case, ctrl in ERCB_STATES:
        mat, grp = plats[stem]
        A = mat[grp[grp.isin(case)].index]
        C = mat[grp[grp.isin(ctrl)].index]
        res = moderated_ttest(A, C)
        human[name] = res["lfc"][np.isfinite(res["lfc"])]
        meta.append({"state": name, "platform": stem.split("-")[1], "layer": "microarray",
                     "n_case": A.shape[1], "n_ctrl": C.shape[1]})
    kp = kpmp_lfc()
    for name, v in kp.items():
        human[name] = v
        meta.append({"state": name, "platform": "KPMP-MS", "layer": "proteome",
                     "n_case": np.nan, "n_ctrl": np.nan})
    META = pd.DataFrame(meta)

    # ---------------------------------------------- 3. 3種共通の遺伝子空間
    ercb_genes = set(plats["GSE104954-GPL22945"][0].index)
    common_ercb = sorted(set(D.index) & ercb_genes)
    kp_genes = set(kp["hs_KPMP_TI"].index) if "hs_KPMP_TI" in kp else set()
    common_kpmp = sorted(set(D.index) & kp_genes)
    gene_report = {
        "Group A intersection (cat+mouse, complete case)": len(D.index),
        "ERCB genes (Entrez CDF)": len(ercb_genes),
        "3-species common (Group A ∩ ERCB)": len(common_ercb),
        "KPMP proteins": len(kp_genes),
        "3-species common (Group A ∩ KPMP)": len(common_kpmp),
    }
    print("=== 3. 遺伝子空間 ===")
    for k, v in gene_report.items():
        print(f"  {k}: {v}")

    G = pd.Index(common_ercb)
    Dh = D.reindex(G)
    H = pd.DataFrame({k: v.reindex(G) for k, v in human.items()})
    keep = Dh.notna().all(axis=1) & H[[c for c in H.columns if c.startswith("hs_ERCB")]].notna().all(axis=1)
    G = G[keep]
    print(f"  ERCB 状態すべてで有限: {len(G)}")

    # ---------------------------------------------- 1. 天井の推定可能性
    rel_rows, sb_draws = [], {}
    for name, stem, case, ctrl in ERCB_STATES:
        mat, grp = plats[stem]
        r, info = split_half_reliability(mat, grp, case, ctrl, G, rng)
        row = {"state": name, **info}
        if r is None:
            row.update({"r_half_median": np.nan, "reliability_SB_median": np.nan,
                        "estimable": False})
        else:
            rh, sb = r
            sb_draws[name] = sb
            row.update({"r_half_median": float(np.median(rh)),
                        "reliability_SB_median": float(np.median(sb)),
                        "reliability_SB_lo95": float(np.percentile(sb, 2.5)),
                        "reliability_SB_hi95": float(np.percentile(sb, 97.5)),
                        "estimable": True})
        rel_rows.append(row)
    for name in kp:
        rel_rows.append({"state": name, "n_case": np.nan, "n_ctrl": np.nan, "n_max": 0,
                         "r_half_median": np.nan, "reliability_SB_median": np.nan,
                         "estimable": False})
    REL = pd.DataFrame(rel_rows)
    REL.round(4).to_csv(OUT / "human_reliability.tsv", sep="\t", index=False)
    print("\n=== 1. 天井の推定可能性 ===")
    print(REL[["state", "n_case", "n_ctrl", "n_max", "r_half_median",
               "reliability_SB_median", "estimable"]].round(3).to_string(index=False))

    # 既存16状態の信頼性（cat/mouse 側）
    main_rel = pd.read_csv(HERE / "results" / "reliability" / "reliability.tsv",
                           sep="\t").set_index("state")

    # ---------------------------------------------- 4. cos の3クラス
    allv = {}
    for c in D.columns:
        allv[c] = Dh.loc[G, c]
    for c in H.columns:
        allv[c] = H.loc[G, c]
    sp = {c: ("cat" if c.startswith("cat") else "human" if c.startswith("hs_") else "mouse")
          for c in allv}

    rows = []
    for a, b in combinations(list(allv), 2):
        u, v = allv[a], allv[b]
        m = u.notna() & v.notna()
        if m.sum() < 200:
            continue
        cls = "".join(sorted([sp[a], sp[b]]))
        label = {"catmouse": "cat x mouse", "cathuman": "cat x human",
                 "humanmouse": "mouse x human", "catcat": "within cat",
                 "mousemouse": "within mouse", "humanhuman": "within human"}[cls]
        ra = main_rel["reliability_SB_median"].get(a, sb_draws.get(a, [np.nan]))
        rb = main_rel["reliability_SB_median"].get(b, sb_draws.get(b, [np.nan]))
        ra = float(np.median(ra)) if not np.isscalar(ra) else float(ra)
        rb = float(np.median(rb)) if not np.isscalar(rb) else float(rb)
        ceil = float(np.sqrt(ra * rb)) if np.isfinite(ra) and np.isfinite(rb) else np.nan
        rows.append({"a": a, "b": b, "class": label, "n_genes": int(m.sum()),
                     "cos": cosine(u[m].to_numpy(), v[m].to_numpy()),
                     "ceiling_SB": ceil,
                     "ceiling_estimable": bool(np.isfinite(ceil)),
                     "layer_mismatch": bool(a.startswith("hs_KPMP") or b.startswith("hs_KPMP"))})
    P = pd.DataFrame(rows)
    P.round(4).to_csv(OUT / "human_pairs.tsv", sep="\t", index=False)

    # 信頼性が低すぎるヒト状態は除外する。GPL24120 は対照が LD3+TN2=5 しかなく、
    # 半分に割ると 2 検体になるため Δ が対照のサンプリング誤差に埋もれる。
    REL_MIN = 0.70
    weak = sorted(REL.loc[REL.estimable & (REL.reliability_SB_median < REL_MIN), "state"])
    strong = sorted(REL.loc[REL.estimable & (REL.reliability_SB_median >= REL_MIN), "state"])
    print(f"\n信頼性 {REL_MIN} 未満で除外するヒト状態: {weak}")
    print(f"残るヒト状態: {strong}")

    def summarise(d_all, tag):
        summ = []
        for cls in ["cat x human", "mouse x human", "cat x mouse"]:
            d = d_all[d_all["class"] == cls]
            if len(d):
                summ.append({"subset": tag, "class": cls, "n_pairs": len(d),
                             "cos_median": d.cos.median(),
                             "cos_min": d.cos.min(), "cos_max": d.cos.max(),
                             "ceiling_median": d.ceiling_SB.median(),
                             "gap_min": (d.ceiling_SB - d.cos).min(),
                             "n_above_ceiling": int((d.cos > d.ceiling_SB).sum())})
        return pd.DataFrame(summ)

    ok = P[P.ceiling_estimable & ~P.layer_mismatch]
    bad_set = set(weak)
    ok_strong = ok[~ok.a.isin(bad_set) & ~ok.b.isin(bad_set)]
    print("\n=== 4a. 全 ERCB 状態（弱い状態を含む・参考） ===")
    print(summarise(ok, "all_ERCB").round(3).to_string(index=False))
    print("\n=== 4b. 信頼性 >= 0.70 のヒト状態のみ（判断はこちらで行う） ===")
    SUM = pd.concat([summarise(ok, "all_ERCB"), summarise(ok_strong, "reliable_only")],
                    ignore_index=True)
    print(summarise(ok_strong, "reliable_only").round(3).to_string(index=False))
    SUM.round(4).to_csv(OUT / "human_cos_summary.tsv", sep="\t", index=False)
    ok_strong.round(4).to_csv(OUT / "human_pairs_reliable.tsv", sep="\t", index=False)

    print("\n--- 状態別の内訳（信頼性 >= 0.70） ---")
    for cls in ["cat x human", "mouse x human"]:
        d = ok_strong[ok_strong["class"] == cls]
        for h in strong:
            dd = d[(d.a == h) | (d.b == h)]
            if len(dd):
                print(f"  {cls:15s} vs {h:18s}: n={len(dd):2d} "
                      f"cos 中央 {dd.cos.median():.3f} [{dd.cos.min():.3f}, {dd.cos.max():.3f}]")

    print("\n--- KPMP（層が違い天井も推定不可。参考値、主張には使わない） ---")
    kpp = P[P.layer_mismatch]
    for cls in ["cat x human", "mouse x human"]:
        d = kpp[kpp["class"] == cls]
        if len(d):
            print(f"  {cls}: n={len(d)}, cos 中央 {d.cos.median():.3f} "
                  f"[{d.cos.min():.3f}, {d.cos.max():.3f}]")

    # ---------------------------------------------- 2. プラットフォーム差
    print("\n=== 2. プラットフォーム差の見積もり ===")
    plat = []
    for a, b, note in [("hs_ERCB_DN_A", "hs_ERCB_DN_B", "同一疾患(DN)・同一組織・別アレイ版"),
                       ("hs_ERCB_allCKD_A", "hs_ERCB_allCKD_B", "全CKD・別アレイ版")]:
        r = P[((P.a == a) & (P.b == b)) | ((P.a == b) & (P.b == a))]
        if len(r):
            plat.append({"pair": f"{a} x {b}", "note": note,
                         "cos": float(r.cos.iloc[0]), "ceiling_SB": float(r.ceiling_SB.iloc[0])})
    within_A = P[(P["class"] == "within human") & P.a.str.endswith("_A") & P.b.str.endswith("_A")]
    if len(within_A):
        plat.append({"pair": "GPL22945 内のヒト状態ペア", "note": "同一アレイ版・別疾患",
                     "cos": float(within_A.cos.median()), "ceiling_SB": float(within_A.ceiling_SB.median())})
    PLAT = pd.DataFrame(plat)
    PLAT.round(4).to_csv(OUT / "platform_effect.tsv", sep="\t", index=False)
    print(PLAT.round(3).to_string(index=False))

    # ---------------------------------------------- 5. cat x human vs mouse x human
    # ペアは 16 の非ヒト状態から作られ独立でないので、p は状態レベルの
    # 網羅並べ替え（C(16,4)=1820）から取る（§2.3 と同じ手続き）。
    hp = ok_strong[ok_strong["class"].isin(["cat x human", "mouse x human"])].copy()
    hp["nh"] = hp.apply(lambda r: r["a"] if r["b"].startswith("hs_") else r["b"], axis=1)

    def auc2(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return float(((x[:, None] > y[None, :]).sum()
                      + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size))

    def perm_test(drop, label):
        d = hp[~hp.nh.isin(drop)]
        st = sorted(d.nh.unique())
        nc = sum(x.startswith("cat") for x in st)

        def stat(cs):
            m = d.nh.isin(cs)
            return auc2(d.cos[m], d.cos[~m]) if m.sum() >= 3 and (~m).sum() >= 3 else np.nan

        obs = stat([x for x in st if x.startswith("cat")])
        null = np.array([v for v in (stat([st[i] for i in c])
                                     for c in combinations(range(len(st)), nc))
                         if np.isfinite(v)])
        return {"subset": label, "n_states": len(st), "n_pairs": len(d), "auc": obs,
                "null_median": float(np.median(null)), "null_p95": float(np.percentile(null, 95)),
                "n_perm": int(null.size), "p_exact": float((null >= obs).sum() / null.size)}

    AT = pd.DataFrame([
        perm_test([], "all mouse states"),
        perm_test(["IRI_2h"], "drop IRI 2h"),
        perm_test(["IRI_2h", "IRI_4h"], "drop acute IRI (2h, 4h)"),
        perm_test(["IRI_2h", "IRI_4h", "IRI_24h", "IRI_48h", "IRI_72h"], "only >= 7 d"),
        perm_test(["IRI_2h", "IRI_4h", "mouse_5D"], "drop acute IRI + Pod-TRECK D5"),
    ])
    AT.round(4).to_csv(OUT / "cat_vs_mouse_to_human.tsv", sep="\t", index=False)
    print("\n=== 5. cat x human vs mouse x human（状態レベル網羅並べ替え） ===")
    print(AT.round(4).to_string(index=False))
    per = hp.groupby("nh").cos.median().sort_values(ascending=False)
    per.round(4).to_frame("cos_median_vs_human").to_csv(OUT / "per_state_to_human.tsv", sep="\t")
    catmed = per[[i for i in per.index if i.startswith("cat")]]
    momed = per[[i for i in per.index if not i.startswith("cat")]]
    print(f"  猫状態の cos 中央: [{catmed.min():.3f}, {catmed.max():.3f}]")
    print(f"  マウス最良状態: {momed.index[0]} = {momed.iloc[0]:.3f}（猫の最低 {catmed.min():.3f} と重なる）")

    pd.Series(gene_report).to_frame("n").to_csv(OUT / "gene_space.tsv", sep="\t")
    META.to_csv(OUT / "human_states.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
