"""ペアの cos θ と天井の差を、個体レベルの不確実性つきで評価する。

遺伝子ブートストラップの区間は遺伝子の入れ替わりに対する不確実性であって、
個体の入れ替わりに対するものではない。マウスの多くは n=3 なので、個体レベルの
誤差はこれより大きい。ここでは個体 ID を単位に復元抽出する。

設計上の注意（2026-09-15 の修正で作り直した点）:

1. 「異なる個体が2頭未満」の判定は**個体 ID** で行う。発現行列の値の一意性で
   判定してはならない（旧実装の誤り。未定義反復が常に0になっていた）。
2. 分割半は復元抽出後の**異なる個体**に対して行う。多重度をそのまま分割すると
   同じ個体が両半分に入り、r_half が上振れして天井を押し上げる。
3. 個体の対応を保つ。対照群を共有する状態は同じ抽出を使う。ネコの皮質と髄質は
   同一個体なので、ID 集合が交わる群はひとつのプールとしてまとめて抽き、
   各群はそのうち自分に属する ID を取る。
4. 未定義になった反復は再抽出せず、割合をそのまま報告する。区間は定義できた
   反復に条件つけた値になる。

併せて、退化の起きない相互確認として leave-one-animal-out ジャックナイフも出す。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402

OUT = HERE / "results" / "pair_uncertainty"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260826
N_BOOT = int(os.environ.get("PU_BOOT", 1000))
N_INNER = int(os.environ.get("PU_INNER", 50))
MIN_DISTINCT = 2


def human_matrix(mat, cols, species, genes):
    out = np.full((len(genes), len(cols)), np.nan)
    idx = pd.Index(genes)
    for j, c in enumerate(cols):
        v, _ = B.to_human(mat[c], species)
        out[:, j] = v.reindex(idx).to_numpy(dtype=float)
    return out


def build_states(genes):
    """状態ごとに (個体 ID -> 遺伝子ベクトル) を case / ctrl で持つ。"""
    S = {}

    def put(state, case_df, ctrl_df, species, ctrl_key):
        S[state] = dict(
            case_ids=list(case_df.columns), ctrl_ids=list(ctrl_df.columns),
            case=human_matrix(case_df, list(case_df.columns), species, genes),
            ctrl=human_matrix(ctrl_df, list(ctrl_df.columns), species, genes),
            species=species, ctrl_key=ctrl_key)

    for state, label, tissue in B.CAT_STATES:
        mat, grp = B._cat_mat(tissue)
        put(state, mat[grp[grp == label].index], mat[grp[grp == "Control"].index],
            "cat", ("cat_control",))          # 皮質と髄質の対照は同一6個体
    mr = B.read_int("mouse_rna"); mgrp = B.read_int("mouse_rna_grp").squeeze()
    gm = mr.T.groupby(mgrp).mean().T
    m = np.log2(mr[gm.max(axis=1) >= B.TH["fpkm_min"]] + 1)
    for state, label, _h in B.PODTRECK_STATES:
        put(state, m[mgrp[mgrp == label].index], m[mgrp[mgrp == "Ctrl"].index],
            "mouse", ("podtreck",))
    imat = B.read_int("mouse_iri_matrix"); igrp = B.read_int("mouse_iri_grp").squeeze()
    for state, label, _h, ctrl_labels in B.IRI_STATES:
        A = np.log2(imat[igrp[igrp == label].index] + 1)
        C = np.log2(imat[igrp[igrp.isin(ctrl_labels)].index] + 1)
        keep = pd.concat([A, C], axis=1).max(axis=1) > 1
        put(state, A[keep], C[keep], "mouse", ("iri", tuple(sorted(ctrl_labels))))
    return S


def pools(groups):
    """ID 集合が交わる群を同じプールにまとめる（union-find）。"""
    keys = list(groups)
    parent = {k: k for k in keys}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if set(groups[a]) & set(groups[b]):
                parent[find(a)] = find(b)
    out = {}
    for k in keys:
        out.setdefault(find(k), []).append(k)
    return list(out.values())


def delta(cols, mats, ids):
    """個体 ID の多重集合から Δ を作る。"""
    cidx = [mats["case_pos"][i] for i in ids["case"]]
    kidx = [mats["ctrl_pos"][i] for i in ids["ctrl"]]
    return mats["case"][:, cidx].mean(1) - mats["ctrl"][:, kidx].mean(1)


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def r_half_distinct(mats, ids, rng, n_inner):
    """異なる個体だけで分割半 cos の中央値。2頭未満なら NaN。"""
    uc = sorted(set(ids["case"])); uk = sorted(set(ids["ctrl"]))
    if len(uc) < MIN_DISTINCT or len(uk) < MIN_DISTINCT:
        return np.nan
    cpos = np.array([mats["case_pos"][i] for i in uc])
    kpos = np.array([mats["ctrl_pos"][i] for i in uk])
    vals = np.empty(n_inner)
    for i in range(n_inner):
        a = rng.permutation(len(cpos)); b = rng.permutation(len(kpos))
        c1, c2 = cpos[a[: len(a) // 2]], cpos[a[len(a) // 2:]]
        k1, k2 = kpos[b[: len(b) // 2]], kpos[b[len(b) // 2:]]
        d1 = mats["case"][:, c1].mean(1) - mats["ctrl"][:, k1].mean(1)
        d2 = mats["case"][:, c2].mean(1) - mats["ctrl"][:, k2].mean(1)
        vals[i] = cosine(d1, d2)
    return float(np.nanmedian(vals))


def main():
    gl = [g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()]
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    genes = sorted(set(gl) & set(D.index))
    genes = [g for g in genes if D.loc[g].notna().all()]
    print(f"遺伝子 {len(genes)}")

    S = build_states(genes)
    for st, d in S.items():
        d["case_pos"] = {i: j for j, i in enumerate(d["case_ids"])}
        d["ctrl_pos"] = {i: j for j, i in enumerate(d["ctrl_ids"])}
        obs = delta(None, d, {"case": d["case_ids"], "ctrl": d["ctrl_ids"]})
        assert np.nanmax(np.abs(obs - D.loc[genes, st].to_numpy())) < 1e-6, f"{st} の Δ が再現しない"
    print("Δ の再現を確認（全16状態）")

    # 退化率の理論値との照合（n=3 群の 3 個体復元抽出で全部同一 = 3(1/3)^3 = 1/9）
    rng0 = np.random.default_rng(12345)
    hit = sum(len(set(rng0.integers(0, 3, 3))) < 2 for _ in range(200000))
    print(f"自己検査: n=3 で異なる個体が2頭未満になる割合 実測 {hit/200000:.4f} / 理論 {1/9:.4f}")

    P = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    CE = pd.read_csv(HERE / "results" / "reliability" / "pairs_ceilings.tsv",
                     sep="\t").set_index(["a", "b"])
    rng = np.random.default_rng(SEED)
    rows = []
    for _, pr in P.iterrows():
        a, b = pr.a, pr.b
        A, Bs = S[a], S[b]
        groups = {("case", a): A["case_ids"], ("ctrl", a): A["ctrl_ids"],
                  ("case", b): Bs["case_ids"], ("ctrl", b): Bs["ctrl_ids"]}
        # 対照を共有する状態は ID 集合が一致するので pools が同じプールにまとめる
        grouping = pools(groups)
        cosb = np.full(N_BOOT, np.nan); ceil = np.full(N_BOOT, np.nan)
        n_undef = 0
        for i in range(N_BOOT):
            drawn = {}
            for pool in grouping:
                ids = sorted({x for k in pool for x in groups[k]})
                take = [ids[j] for j in rng.integers(0, len(ids), len(ids))]
                for k in pool:
                    member = set(groups[k])
                    drawn[k] = [x for x in take if x in member]
            if any(len(v) == 0 for v in drawn.values()):
                n_undef += 1
                continue
            ia = {"case": drawn[("case", a)], "ctrl": drawn[("ctrl", a)]}
            ib = {"case": drawn[("case", b)], "ctrl": drawn[("ctrl", b)]}
            cosb[i] = cosine(delta(None, A, ia), delta(None, Bs, ib))
            ra = r_half_distinct(A, ia, rng, N_INNER)
            rb = r_half_distinct(Bs, ib, rng, N_INNER)
            if np.isfinite(ra) and np.isfinite(rb):
                ceil[i] = np.sqrt(max(ra, 0.0) * max(rb, 0.0))
            else:
                n_undef += 1
        gap = ceil - cosb
        ok = np.isfinite(gap)

        # --- ジャックナイフ（個体を1頭ずつ抜く。退化しない）---
        jack = []
        allids = sorted({(k, x) for k, v in groups.items() for x in v})
        for k_drop, x_drop in allids:
            ids = {k: [x for x in v if not (k == k_drop and x == x_drop)]
                   for k, v in groups.items()}
            if any(len(v) < 2 for v in ids.values()):
                continue
            ja = {"case": ids[("case", a)], "ctrl": ids[("ctrl", a)]}
            jb = {"case": ids[("case", b)], "ctrl": ids[("ctrl", b)]}
            c = cosine(delta(None, A, ja), delta(None, Bs, jb))
            ra = r_half_distinct(A, ja, rng, N_INNER)
            rb = r_half_distinct(Bs, jb, rng, N_INNER)
            if np.isfinite(ra) and np.isfinite(rb):
                jack.append(np.sqrt(max(ra, 0.0) * max(rb, 0.0)) - c)
        jack = np.array(jack)

        rows.append(dict(
            a=a, b=b, cls=pr["class"], shared_control=bool(pr.shared_control),
            cos_obs=float(pr.cos), ceiling_obs=float(CE.loc[(a, b), "ceiling_raw"]),
            n_ok=int(ok.sum()), frac_undefined=float(n_undef / N_BOOT),
            cos_lo=float(np.nanpercentile(cosb, 2.5)) if ok.any() else np.nan,
            cos_hi=float(np.nanpercentile(cosb, 97.5)) if ok.any() else np.nan,
            ceil_lo=float(np.nanpercentile(ceil, 2.5)) if ok.any() else np.nan,
            ceil_hi=float(np.nanpercentile(ceil, 97.5)) if ok.any() else np.nan,
            gap_median=float(np.nanmedian(gap)) if ok.any() else np.nan,
            gap_lo=float(np.nanpercentile(gap, 2.5)) if ok.any() else np.nan,
            gap_hi=float(np.nanpercentile(gap, 97.5)) if ok.any() else np.nan,
            p_gap_le0=float(np.mean(gap[ok] <= 0)) if ok.any() else np.nan,
            jack_n=len(jack), jack_min=float(jack.min()) if len(jack) else np.nan,
            jack_median=float(np.median(jack)) if len(jack) else np.nan))
    T = pd.DataFrame(rows)
    T.round(4).to_csv(OUT / "pair_uncertainty.tsv", sep="\t", index=False)

    cs = T[T.cls == "cross_species"]
    ws = T[(T.cls != "cross_species") & (~T.shared_control)]
    print(f"\n異種ペア {len(cs)} 組（個体ブートストラップ {N_BOOT} 回、内側分割 {N_INNER} 回）")
    print(f"  未定義（分割不能）反復の割合: 中央 {cs.frac_undefined.median():.3f} "
          f"[{cs.frac_undefined.min():.3f}, {cs.frac_undefined.max():.3f}]")
    print(f"  gap = 天井 − cos の中央値: {cs.gap_median.median():.3f}")
    print(f"  gap の 2.5%点の最小: {cs.gap_lo.min():.3f}")
    print(f"  gap の区間が 0 を含まないペア: {int((cs.gap_lo > 0).sum())} / {len(cs)}")
    print(f"  P(gap <= 0) の最大: {cs.p_gap_le0.max():.3f}")
    print(f"  観測 cos の 97.5%点の最大: {cs.cos_hi.max():.3f}")
    print(f"  天井の 2.5%点の最小: {cs.ceil_lo.min():.3f}")
    print(f"  ジャックナイフ gap の最小: {cs.jack_min.min():.3f} "
          f"（0 以下のペア {int((cs.jack_min <= 0).sum())} / {len(cs)}）")
    print(f"\n対照非共有の同種ペア {len(ws)} 組: gap 中央 {ws.gap_median.median():.3f}, "
          f"区間が 0 を含まない {int((ws.gap_lo > 0).sum())} / {len(ws)}")
    print(f"  同 未定義割合: 中央 {ws.frac_undefined.median():.3f}")
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
