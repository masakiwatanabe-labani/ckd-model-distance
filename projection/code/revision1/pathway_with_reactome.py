"""Task 3. Reactome を加えて §2.3 の主要数値を再計算する。

既存の結果は書き換えず、results/revision1/ に新しいファイルとして出す。
セットの取り方・サイズ基準・遺伝子空間・ペアの制限はすべて既存と同一。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402
from pathway_cos import read_gmt, MIN_GENES  # noqa: E402

OUT = HERE / "results" / "revision1"
OUT.mkdir(parents=True, exist_ok=True)
COLL4 = {"GO_BP": "geneset_gobp.gmt", "KEGG": "geneset_kegg.gmt",
         "Hallmark": "geneset_hallmark.gmt", "Reactome": "geneset_reactome.gmt"}
SEED = 20260826
N_RAND = 300


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def load():
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    P = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
    return D, P


def collect(D, colls):
    genes = set(D.index)
    gpos = {g: i for i, g in enumerate(D.index)}
    out = []
    for coll, fn in colls.items():
        for name, members in read_gmt(B.REF / fn).items():
            gi = np.array(sorted(gpos[g] for g in members & genes))
            if len(gi) >= MIN_GENES:
                out.append((coll, name, gi))
    return out


def analyse(D, P, paths, label):
    X = D.to_numpy(float)
    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()
    wi = P.wi.to_numpy()

    def cosmat(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        return (Mn.T @ Mn)[ia, ib]

    rows = []
    for coll, name, gi in paths:
        c = cosmat(X[gi])
        rows.append({"collection": coll, "pathway": name, "n_genes": len(gi),
                     "auc": auc(c[m & wi], c[m & ~wi])})
    R = pd.DataFrame(rows)

    C = X - X.mean(axis=0)
    cond = {"cos_raw": cosmat(X), "cos_centred": cosmat(C),
            "cos_agg_centred": cosmat(np.vstack([C[gi].mean(axis=0) for _c, _n, gi in paths])),
            "cos_agg_raw": cosmat(np.vstack([X[gi].mean(axis=0) for _c, _n, gi in paths]))}
    summ = {"label": label, "n_sets": len(paths),
            "n_within": int((m & wi).sum()), "n_cross": int((m & ~wi).sum()),
            "pathway_auc_median": R.auc.median(),
            "pathway_auc_q25": R.auc.quantile(.25), "pathway_auc_q75": R.auc.quantile(.75),
            "frac_ge_090": float((R.auc >= 0.90).mean()), "n_le_060": int((R.auc <= 0.60).sum()),
            "size_auc_spearman": float(pd.Series(R.n_genes).corr(R.auc, method="spearman"))}
    for k, v in cond.items():
        summ[f"auc[{k}]"] = auc(v[m & wi], v[m & ~wi])
        summ[f"within_med[{k}]"] = float(np.median(v[m & wi]))
        summ[f"cross_med[{k}]"] = float(np.median(v[m & ~wi]))
    return R, summ


def main():
    D, P = load()
    base = collect(D, {k: v for k, v in COLL4.items() if k != "Reactome"})
    full = collect(D, COLL4)
    print(f"既存 {len(base)} セット / Reactome 追加後 {len(full)} セット")
    per_coll = pd.Series([c for c, _n, _g in full]).value_counts()
    print(per_coll.to_string())

    # KEGG との重複（Jaccard）
    kegg = {n: set(g) for c, n, g in full if c == "KEGG"}
    reac = {n: set(g) for c, n, g in full if c == "Reactome"}
    best = []
    for rn, rg in reac.items():
        j = max((len(rg & kg) / len(rg | kg) for kg in kegg.values()), default=0.0)
        best.append(j)
    print(f"Reactome 各セットの KEGG との最大 Jaccard: 中央 {np.median(best):.3f}, "
          f">=0.5 の割合 {np.mean(np.array(best) >= 0.5):.3f}")

    Rb, sb = analyse(D, P, base, "GO+KEGG+Hallmark (193)")
    Rf, sf = analyse(D, P, full, f"+ Reactome ({len(full)})")
    S = pd.DataFrame([sb, sf])
    S.round(4).to_csv(OUT / "pathway_with_reactome_summary.tsv", sep="\t", index=False)
    Rf.round(4).to_csv(OUT / "pathway_with_reactome_per_set.tsv", sep="\t", index=False)
    pd.DataFrame({"reactome_set": list(reac), "max_jaccard_vs_kegg": best}).round(4) \
        .to_csv(OUT / "reactome_kegg_overlap.tsv", sep="\t", index=False)
    cols = ["label", "n_sets", "pathway_auc_median", "pathway_auc_q25", "pathway_auc_q75",
            "n_le_060", "size_auc_spearman", "auc[cos_raw]", "auc[cos_centred]",
            "auc[cos_agg_centred]", "auc[cos_agg_raw]"]
    print()
    print(S[cols].round(4).to_string(index=False))
    print(f"\n書き出し: {OUT}")


if __name__ == "__main__":
    main()
