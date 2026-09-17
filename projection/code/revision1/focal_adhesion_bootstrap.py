"""Part D. focal adhesion の3セルに個体ブートストラップの区間を付ける。

Figure 3A の帯は「同じ大きさのランダム遺伝子集合を選び直したときのスコアのばらつき」で
あって、特定の経路について動物を取り直したときの推定誤差ではない。本文が名指しする
3 セル（IRI 7 d / 72 h / 28 d の focal adhesion）について、後者を出して並べる。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pair_uncertainty as PU  # noqa: E402
import build_delta_matrix as B  # noqa: E402
from pathway_cos import read_gmt, COLLECTIONS  # noqa: E402

OUT = HERE / "results" / "round5"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260917
N_BOOT = int(os.environ.get("FA_BOOT", 2000))
REF = "cat_CKD34"
STATES = ["IRI_7d", "IRI_72h", "IRI_28d"]
PATHWAY = "Focal adhesion"


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    gl = [g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()]
    genes = sorted(set(gl) & set(D.index))
    genes = [g for g in genes if D.loc[g].notna().all()]

    sets = {}
    for coll, fn in COLLECTIONS.items():
        for name, members in read_gmt(B.REF / fn).items():
            sets[(coll, name)] = members
    key = [k for k in sets if k[1] == PATHWAY]
    assert key, f"{PATHWAY} が見つからない"
    coll = key[0][0]
    members = sorted(set(sets[key[0]]) & set(genes))
    print(f"{PATHWAY} ({coll}): Group A 空間で {len(members)} 遺伝子")
    idx = [genes.index(g) for g in members]

    S = PU.build_states(genes)
    for st in [REF] + STATES:
        S[st]["case_pos"] = {i: j for j, i in enumerate(S[st]["case_ids"])}
        S[st]["ctrl_pos"] = {i: j for j, i in enumerate(S[st]["ctrl_ids"])}

    def obs_delta(st):
        return PU.delta(None, S[st], {"case": S[st]["case_ids"], "ctrl": S[st]["ctrl_ids"]})

    rng = np.random.default_rng(SEED)
    rows = []
    NUL = pd.read_csv(HERE / "results" / "revision1" / "pathway_reactome"
                      / "random_set_null.tsv", sep="\t")
    for st in STATES:
        groups = {("case", REF): S[REF]["case_ids"], ("ctrl", REF): S[REF]["ctrl_ids"],
                  ("case", st): S[st]["case_ids"], ("ctrl", st): S[st]["ctrl_ids"]}
        grouping = PU.pools(groups)
        vals = np.full(N_BOOT, np.nan)
        undef = 0
        for i in range(N_BOOT):
            drawn = {}
            for pool in grouping:
                ids = sorted({x for k in pool for x in groups[k]})
                take = [ids[j] for j in rng.integers(0, len(ids), len(ids))]
                for k in pool:
                    member = set(groups[k])
                    drawn[k] = [x for x in take if x in member]
            if any(len(v) == 0 for v in drawn.values()):
                undef += 1
                continue
            a = PU.delta(None, S[REF], {"case": drawn[("case", REF)], "ctrl": drawn[("ctrl", REF)]})
            v = PU.delta(None, S[st], {"case": drawn[("case", st)], "ctrl": drawn[("ctrl", st)]})
            vals[i] = PU.cosine(v[idx], a[idx])
        v = vals[np.isfinite(vals)]
        a0, v0 = obs_delta(REF), obs_delta(st)
        point = PU.cosine(v0[idx], a0[idx])
        sizes = sorted(NUL.n_genes.unique(), key=lambda x: abs(x - len(members)))
        n = NUL[(NUL.state == st) & (NUL.n_genes == sizes[0])]
        hi95 = float(n.hi95.iloc[0]) if len(n) else np.nan
        rows.append(dict(state=st, pathway=PATHWAY, n_genes=len(members),
                         cos_point=round(point, 4),
                         boot_median=round(float(np.median(v)), 4),
                         boot_lo95=round(float(np.percentile(v, 2.5)), 4),
                         boot_hi95=round(float(np.percentile(v, 97.5)), 4),
                         frac_undefined=round(undef / N_BOOT, 4),
                         random_set_hi95_matched_size=round(hi95, 4),
                         boot_lo95_above_random_hi95=bool(
                             np.percentile(v, 2.5) > hi95) if np.isfinite(hi95) else None))
        print(f"  {st}: 点 {point:.3f}  個体ブートストラップ中央 {rows[-1]['boot_median']:.3f} "
              f"[{rows[-1]['boot_lo95']:.3f}, {rows[-1]['boot_hi95']:.3f}]  "
              f"同サイズランダム集合の 97.5 パーセンタイル {hi95:.3f}")
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "focal_adhesion_animal_bootstrap.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'focal_adhesion_animal_bootstrap.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
