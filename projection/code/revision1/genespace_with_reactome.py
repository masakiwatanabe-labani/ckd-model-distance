"""Task 3 の続き。3つの遺伝子空間での経路単位 / 集約 AUC を Reactome 込みで再計算する。"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_delta_matrix as B  # noqa: E402
from pathway_cos import read_gmt, COLLECTIONS, MIN_GENES  # noqa: E402

OUT = HERE / "results" / "revision1"
SPACES = [("Group A", "groupA_intersection.txt", "control_analysis"),
          ("all 1:1 orthologues", "aa_ortholog_all.txt", "control_analysis_ortholog_all"),
          ("matched Group B", "aa_groupB_matched.txt", "control_analysis_groupB_matched")]


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def main():
    rows = []
    for label, genefile, pairdir in SPACES:
        D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
        want = {g.strip() for g in (HERE / genefile).read_text().split() if g.strip()}
        D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
        gpos = {g: i for i, g in enumerate(D.index)}
        genes = set(D.index)
        paths = []
        for coll, fn in COLLECTIONS.items():
            for name, members in read_gmt(B.REF / fn).items():
                gi = np.array(sorted(gpos[g] for g in members & genes))
                if len(gi) >= MIN_GENES:
                    paths.append(gi)
        P = pd.read_csv(HERE / "results" / pairdir / "pairs.tsv", sep="\t")
        P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
        idx = {s: i for i, s in enumerate(D.columns)}
        ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
        m = (~P.shared_control).to_numpy()
        wi = P.wi.to_numpy()
        X = D.to_numpy(float)

        def cosmat(M):
            Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
            return (Mn.T @ Mn)[ia, ib]

        per = [auc(c[m & wi], c[m & ~wi]) for c in (cosmat(X[gi]) for gi in paths)]
        C = X - X.mean(axis=0)
        agg = cosmat(np.vstack([C[gi].mean(axis=0) for gi in paths]))
        whole = cosmat(X)
        rows.append({"space": label, "n_genes": len(D), "n_sets": len(paths),
                     "whole_panel_auc": auc(whole[m & wi], whole[m & ~wi]),
                     "pathway_auc_median": float(np.median(per)),
                     "aggregated_auc_centred": auc(agg[m & wi], agg[m & ~wi])})
    T = pd.DataFrame(rows)
    T.round(4).to_csv(OUT / "genespace_with_reactome.tsv", sep="\t", index=False)
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
