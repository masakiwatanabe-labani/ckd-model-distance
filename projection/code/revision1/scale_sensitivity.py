"""Part 3-2. 入力尺度（offset と低発現フィルタ）を変えたときに、主要なコサイン量と
経路集約後の AUC がどう動くかを調べる。

Table S2 は Spearman 相関の感度分析にとどまっていた。log2(FPKM+1) と
log2(count+1) では同じ +1 でも遺伝子ごとの縮小効果が違い、それは振幅だけでなく
方向と平均中心化の効果にも及ぶ。ここでは Δ 行列を offset とフィルタを変えて作り直し、
中心結論を直接支える 4 つの量を再計算する。

  - 猫皮質 CKD1/2 と髄質 CKD1/2 の alpha / cos / nrm（順位逆転が保たれるか）
  - 全 120 ペアのうち対照を共有しない 87 ペアでの種内 vs 種間 AUC（遺伝子単位）
  - 同じ 87 ペアを経路平均に集約したときの AUC、中心化あり／なし
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_delta_matrix as B  # noqa: E402
from pathway_cos import read_gmt, COLLECTIONS, MIN_GENES  # noqa: E402

OUT = HERE / "results" / "round7"
OUT.mkdir(parents=True, exist_ok=True)
REF = "cat_CKD34"
OFFSETS = [1.0, 0.5, 0.1]
FILTERS = [1.0, 0.0, 5.0]        # 群平均 FPKM の下限。既定は 1.0


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def build_delta(offset, fpkm_min):
    """既定の手続きと同じ順序で Δ 行列を作る。offset とフィルタだけを差し替える。"""
    cols = {}
    for state, label, tissue in B.CAT_STATES:                  # 猫は寄託 log2 値、offset 非該当
        mat, grp = B._cat_mat(tissue)
        case = mat[grp[grp == label].index]
        ctrl = mat[grp[grp == "Control"].index]
        cols[state] = case.mean(axis=1) - ctrl.mean(axis=1)
    mr = B.read_int("mouse_rna"); mgrp = B.read_int("mouse_rna_grp").squeeze()
    gm = mr.T.groupby(mgrp).mean().T
    keep = gm.max(axis=1) >= fpkm_min
    m = np.log2(mr[keep] + offset)
    for state, label, _h in B.PODTRECK_STATES:
        cols[state] = (m[mgrp[mgrp == label].index].mean(axis=1)
                       - m[mgrp[mgrp == "Ctrl"].index].mean(axis=1))
    imat = B.read_int("mouse_iri_matrix"); igrp = B.read_int("mouse_iri_grp").squeeze()
    for state, label, _h, ctrl_labels in B.IRI_STATES:
        A = np.log2(imat[igrp[igrp == label].index] + offset)
        Cc = np.log2(imat[igrp[igrp.isin(ctrl_labels)].index] + offset)
        k = pd.concat([A, Cc], axis=1).max(axis=1) > np.log2(1 + offset)
        cols[state] = A[k].mean(axis=1) - Cc[k].mean(axis=1)
    out = {}
    for state, v in cols.items():
        species = "cat" if state.startswith("cat") else "mouse"
        hv, _ = B.to_human(v, species)
        out[state] = hv
    return pd.DataFrame(out)


def main() -> int:
    genes0 = [g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()]
    pairs = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    keep = pairs[~pairs.shared_control]
    sets = {}
    for coll, fn in COLLECTIONS.items():
        for name, members in read_gmt(B.REF / fn).items():
            sets[(coll, name)] = members
    rows = []
    for offset in OFFSETS:
        for fmin in FILTERS:
            D = build_delta(offset, fmin)
            D = D.loc[D.index.intersection(sorted(set(genes0)))].dropna(axis=0, how="any")
            n = len(D)
            a = D[REF].to_numpy()
            def st(s):
                v = D[s].to_numpy()
                c = cosine(v, a)
                nr = float(np.linalg.norm(v) / np.linalg.norm(a))
                return float(v @ a / (a @ a)), c, nr
            ac, cc, nc = st("cat_CKD12")
            am, cm, nm = st("cat_med_CKD12")

            X = D.to_numpy(float)
            idx = {s: i for i, s in enumerate(D.columns)}
            gene_cos, within = [], []
            for _, r in keep.iterrows():
                gene_cos.append(cosine(X[:, idx[r.a]], X[:, idx[r.b]]))
                within.append(r["class"] != "cross_species")
            within = np.array(within); gene_cos = np.array(gene_cos)
            auc_gene = auc(gene_cos[within], gene_cos[~within])

            Z = (X - X.mean(axis=0)) / X.std(axis=0)          # 状態ごとに遺伝子方向で標準化
            have = set(D.index)
            mats = {}
            for lab, mat in [("centred", Z), ("raw", X)]:
                rowsP = []
                for (coll, name), members in sets.items():
                    ix = [D.index.get_loc(g) for g in (members & have)]
                    if len(ix) >= MIN_GENES:
                        rowsP.append(mat[ix, :].mean(axis=0))
                mats[lab] = np.array(rowsP)
            out = {}
            for lab, Mt in mats.items():
                vals = [cosine(Mt[:, idx[r.a]], Mt[:, idx[r.b]]) for _, r in keep.iterrows()]
                vals = np.array(vals)
                out[lab] = auc(vals[within], vals[~within])
            rows.append(dict(
                offset=offset, fpkm_min=fmin, n_genes=n, n_sets=mats["centred"].shape[0],
                alpha_cortex=round(ac, 4), cos_cortex=round(cc, 4), nrm_cortex=round(nc, 4),
                alpha_medulla=round(am, 4), cos_medulla=round(cm, 4), nrm_medulla=round(nm, 4),
                inversion_holds=bool(am > ac and cm < cc),
                alpha_medulla_above_1=bool(am > 1.0),
                auc_gene_level=round(auc_gene, 4),
                auc_pathway_centred=round(out["centred"], 4),
                auc_pathway_uncentred=round(out["raw"], 4),
                centring_drop=round(out["raw"] - out["centred"], 4)))
            print(f"offset {offset} filter {fmin}: genes {n} sets {rows[-1]['n_sets']} "
                  f"gene AUC {auc_gene:.3f} centred {out['centred']:.3f} "
                  f"uncentred {out['raw']:.3f} inversion {rows[-1]['inversion_holds']}",
                  flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "scale_sensitivity.tsv", sep="\t", index=False)
    base = T[(T.offset == 1.0) & (T.fpkm_min == 1.0)].iloc[0]
    print("\n既定 (offset 1, filter 1) に対する振れ幅:")
    for c in ("auc_gene_level", "auc_pathway_centred", "auc_pathway_uncentred", "centring_drop"):
        print(f"  {c:24s} 既定 {base[c]:.3f}  範囲 {T[c].min():.3f} - {T[c].max():.3f}")
    print(f"  順位逆転が保たれた条件: {int(T.inversion_holds.sum())} / {len(T)}")
    print(f"  alpha > 1 が保たれた条件: {int(T.alpha_medulla_above_1.sum())} / {len(T)}")
    print(f"\n書き出し: {OUT / 'scale_sensitivity.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
