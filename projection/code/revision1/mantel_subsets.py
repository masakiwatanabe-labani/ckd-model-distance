"""Part E-1. 経過時間との関連が、部分集合でも保たれるかを確かめる。

本文は12のモデル状態で「経過時間が方向を約2倍強く追う」と書いているが、
12状態には Pod-TRECK（糸球体起点）と IRI（尿細管起点）が混ざっており、
起点とモデル種別が完全に交絡している。2つの部分集合で確認する。

  - IRI 単独（9 状態）: モデルが1つなので起点とモデルの交絡が無い
  - 12 か月を除く（11 状態）: 時間軸の端の1点に引かれていないか

維持されなければ、その事実を本文に書く。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import mantel_s1 as MS  # noqa: E402

OUT = HERE / "results" / "round5"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rng = np.random.default_rng(MS.SEED)
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in MS.GENES.read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    print(f"遺伝子 {len(D)}")

    tm = pd.read_csv(HERE / "time_map.tsv", sep="\t")
    hours = tm.set_index("state")["hours"].astype(float).to_dict()
    model_states = [s for s in D.columns if np.isfinite(hours.get(s, np.nan))]
    iri = [s for s in model_states if s.startswith("IRI_")]
    no12 = [s for s in model_states if s != "IRI_12mo"]
    iri_no12 = [s for s in iri if s != "IRI_12mo"]

    out = []
    for states, label in [(model_states, "12 model states (as published)"),
                          (iri, "IRI only (9 states)"),
                          (no12, "12 model states minus IRI 12 mo (11 states)"),
                          (iri_no12, "IRI only minus 12 mo (8 states)")]:
        R, conf, _, _ = MS.run(D, states, label, hours, rng)
        out.append(R)
        print(f"\n=== {label}: n = {len(states)} ===")
        d = R[R.predictor == "elapsed time |Δlog10 h|"]
        for _, r in d.iterrows():
            print(f"  {r.dissimilarity:28s} rho {r.mantel_rho:+.3f}  p = {r.p_perm:.4f}")
    T = pd.concat(out, ignore_index=True)
    T.round(4).to_csv(OUT / "mantel_subsets.tsv", sep="\t", index=False)

    print("\n--- 経過時間 x 方向 (D_cos) と 振幅 (D_nrm) の比 ---")
    for lab in T.subset.unique():
        d = T[(T.subset == lab) & (T.predictor == "elapsed time |Δlog10 h|")]
        c = d[d.dissimilarity == "D_cos (1-cos)"]
        n = d[d.dissimilarity == "D_nrm (|log2 norm ratio|)"]
        if len(c) and len(n):
            cv, nv = float(c.mantel_rho.iloc[0]), float(n.mantel_rho.iloc[0])
            ratio = cv / nv if nv != 0 else np.nan
            print(f"  {lab:46s} 方向 {cv:+.3f} (p={float(c.p_perm.iloc[0]):.4f})  "
                  f"振幅 {nv:+.3f} (p={float(n.p_perm.iloc[0]):.4f})  比 {ratio:+.2f}")
    print(f"\n書き出し: {OUT / 'mantel_subsets.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
