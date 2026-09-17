"""Part B-2. 順位逆転が個体の入れ替わりに対して安定かを調べる。

Figure 1 の区間は遺伝子ブートストラップで、遺伝子の入れ替わりに対する不確実性しか
表していない。α = ν cos θ という代数から逆転が起こりうることと、このコホートで
安定に観測されることは別なので、個体単位で再標本化して維持率を出す。

設計の保持:
  - 基準軸（猫皮質 CKD3/4）と皮質 CKD1/2 は同じ 6 頭の対照を共有する。
  - 皮質と髄質は同一個体から採られている。
  したがって ID 集合が交わる群はひとつのプールとして一度だけ抽き、各群は
  そのうち自分に属する ID を取る（pair_uncertainty.pools と同じ扱い）。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import pair_uncertainty as PU  # noqa: E402

OUT = HERE / "results" / "round5"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260916
N_BOOT = int(os.environ.get("RV_BOOT", 2000))

REF = "cat_CKD34"
CORTEX = "cat_CKD12"
MEDULLA = "cat_med_CKD12"


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    gl = [g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()]
    genes = sorted(set(gl) & set(D.index))
    genes = [g for g in genes if D.loc[g].notna().all()]
    print(f"遺伝子 {len(genes)}")

    S = PU.build_states(genes)
    for st in (REF, CORTEX, MEDULLA):
        S[st]["case_pos"] = {i: j for j, i in enumerate(S[st]["case_ids"])}
        S[st]["ctrl_pos"] = {i: j for j, i in enumerate(S[st]["ctrl_ids"])}
        obs = PU.delta(None, S[st], {"case": S[st]["case_ids"], "ctrl": S[st]["ctrl_ids"]})
        assert np.nanmax(np.abs(obs - D.loc[genes, st].to_numpy())) < 1e-6, st
    print("Δ の再現を確認（3 状態）")

    groups = {}
    for st in (REF, CORTEX, MEDULLA):
        groups[("case", st)] = S[st]["case_ids"]
        groups[("ctrl", st)] = S[st]["ctrl_ids"]
    grouping = PU.pools(groups)
    print("同じプールにまとめた群:")
    for pool in grouping:
        print("   ", sorted(f"{k[1]}/{k[0]}" for k in pool),
              f"  個体 {len(sorted({x for k in pool for x in groups[k]}))} 頭")

    def stats(dref, dv):
        n2 = float(dref @ dref)
        alpha = float(dv @ dref) / n2
        nrm = float(np.linalg.norm(dv) / np.linalg.norm(dref))
        cos = PU.cosine(dv, dref)
        return alpha, cos, nrm

    # 観測値
    obs = {st: PU.delta(None, S[st], {"case": S[st]["case_ids"], "ctrl": S[st]["ctrl_ids"]})
           for st in (REF, CORTEX, MEDULLA)}
    oc = stats(obs[REF], obs[CORTEX])
    om = stats(obs[REF], obs[MEDULLA])
    print(f"\n観測: 皮質 CKD1/2  α={oc[0]:.3f} cos={oc[1]:.3f} ν={oc[2]:.3f}")
    print(f"      髄質 CKD1/2  α={om[0]:.3f} cos={om[1]:.3f} ν={om[2]:.3f}")

    rng = np.random.default_rng(SEED)
    rec = {k: np.full(N_BOOT, np.nan) for k in
           ("a_cortex", "a_med", "c_cortex", "c_med", "n_cortex", "n_med")}
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
        dd = {st: PU.delta(None, S[st], {"case": drawn[("case", st)], "ctrl": drawn[("ctrl", st)]})
              for st in (REF, CORTEX, MEDULLA)}
        ac, cc, nc = stats(dd[REF], dd[CORTEX])
        am, cm, nm = stats(dd[REF], dd[MEDULLA])
        rec["a_cortex"][i], rec["c_cortex"][i], rec["n_cortex"][i] = ac, cc, nc
        rec["a_med"][i], rec["c_med"][i], rec["n_med"][i] = am, cm, nm

    R = pd.DataFrame(rec).dropna()
    both = ((R.a_med > R.a_cortex) & (R.c_med < R.c_cortex)).mean()
    only_a = (R.a_med > R.a_cortex).mean()
    only_c = (R.c_med < R.c_cortex).mean()
    above1 = (R.a_med > 1.0).mean()
    rows = [
        ("replicates drawn", N_BOOT),
        ("replicates undefined (some group empty)", n_undef),
        ("replicates used", len(R)),
        ("alpha medulla > alpha cortex", round(float(only_a), 4)),
        ("cos medulla < cos cortex", round(float(only_c), 4)),
        ("both, i.e. the reversal holds", round(float(both), 4)),
        ("alpha medulla > 1.000 (passes the reference)", round(float(above1), 4)),
    ]
    for q, col in [("alpha cortex", "a_cortex"), ("alpha medulla", "a_med"),
                   ("cos cortex", "c_cortex"), ("cos medulla", "c_med"),
                   ("nu cortex", "n_cortex"), ("nu medulla", "n_med")]:
        rows += [(f"{q} median", round(float(R[col].median()), 4)),
                 (f"{q} 2.5th percentile", round(float(np.percentile(R[col], 2.5)), 4)),
                 (f"{q} 97.5th percentile", round(float(np.percentile(R[col], 97.5)), 4))]
    rows.append(("alpha medulla - alpha cortex, median",
                 round(float((R.a_med - R.a_cortex).median()), 4)))
    rows.append(("alpha medulla - alpha cortex, 2.5th percentile",
                 round(float(np.percentile(R.a_med - R.a_cortex, 2.5)), 4)))
    rows.append(("cos cortex - cos medulla, median",
                 round(float((R.c_cortex - R.c_med).median()), 4)))
    rows.append(("cos cortex - cos medulla, 2.5th percentile",
                 round(float(np.percentile(R.c_cortex - R.c_med, 2.5)), 4)))
    T = pd.DataFrame(rows, columns=["quantity", "value"])
    T.to_csv(OUT / "reversal_animal_bootstrap.tsv", sep="\t", index=False)
    R.round(5).to_csv(OUT / "reversal_animal_bootstrap_replicates.tsv", sep="\t", index=False)
    print("\n" + T.to_string(index=False))
    print(f"\n書き出し: {OUT / 'reversal_animal_bootstrap.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
