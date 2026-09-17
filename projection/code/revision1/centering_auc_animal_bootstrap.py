"""Part 5-3. 中心化による経路集約後 AUC の変化（0.868 → 0.511）が、動物を取り直しても
どの程度再現するかを調べる。

中心結果として扱う以上、遺伝子の入れ替わりではなく個体の入れ替わりに対して安定かを
見ておく必要がある。設計は保持する: 対照群を共有する状態は同じ抽出を受け、猫の皮質と
髄質は同一個体なので一緒に抽く（pair_uncertainty.pools と同じ扱い）。

各反復で 16 状態の Δ を作り直し、対照を共有しない 87 ペアについて
  - 中心化なしの経路平均どうしの cos θ
  - 状態ごとに遺伝子方向で中心化してから経路平均した cos θ
を計算し、種内 39 対 種間 48 の AUC をそれぞれ出す。
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
from pathway_cos import read_gmt, MIN_GENES  # noqa: E402

# pair_uncertainty を先に import すると projection/pathway_cos.py（3 コレクション、193 セット）が
# 優先されてしまう。本文の数値は Reactome 込みの 422 セットなので、ここで明示する。
COLLECTIONS = {"Reactome": "geneset_reactome.gmt", "GO_BP": "geneset_gobp.gmt",
               "KEGG": "geneset_kegg.gmt", "Hallmark": "geneset_hallmark.gmt"}

OUT = HERE / "results" / "round7"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260917
N_BOOT = int(os.environ.get("CA_BOOT", 1000))


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def main() -> int:
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    gl = [g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()]
    genes = sorted(set(gl) & set(D.index))
    genes = [g for g in genes if D.loc[g].notna().all()]
    print(f"遺伝子 {len(genes)}")

    S = PU.build_states(genes)
    states = list(D.columns)
    for st in states:
        S[st]["case_pos"] = {i: j for j, i in enumerate(S[st]["case_ids"])}
        S[st]["ctrl_pos"] = {i: j for j, i in enumerate(S[st]["ctrl_ids"])}

    # 経路のメンバーシップ行列
    sets = {}
    for coll, fn in COLLECTIONS.items():
        for name, members in read_gmt(B.REF / fn).items():
            sets[(coll, name)] = members
    gi = {g: i for i, g in enumerate(genes)}
    rows = []
    for key, members in sets.items():
        ix = [gi[g] for g in members if g in gi]
        if len(ix) >= MIN_GENES:
            rows.append(ix)
    print(f"経路セット {len(rows)}")

    pairs = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    keep = pairs[~pairs.shared_control].reset_index(drop=True)
    within = (keep["class"] != "cross_species").to_numpy()
    print(f"対照非共有ペア {len(keep)}（種内 {int(within.sum())} / 種間 {int((~within).sum())}）")

    groups = {}
    for st in states:
        groups[("case", st)] = S[st]["case_ids"]
        groups[("ctrl", st)] = S[st]["ctrl_ids"]
    grouping = PU.pools(groups)
    print(f"プール {len(grouping)}")

    def aggregate(X):
        """状態 x 遺伝子の行列から、中心化あり／なしの経路平均行列を作る。"""
        out = {}
        Z = (X - X.mean(axis=0)) / X.std(axis=0)
        for lab, mat in (("raw", X), ("centred", Z)):
            out[lab] = np.array([mat[ix, :].mean(axis=0) for ix in rows])
        return out

    def aucs(X):
        M = aggregate(X)
        res = {}
        for lab, Mt in M.items():
            v = np.array([cosine(Mt[:, states.index(r.a)], Mt[:, states.index(r.b)])
                          for _, r in keep.iterrows()])
            res[lab] = auc(v[within], v[~within])
        return res

    X0 = D.loc[genes, states].to_numpy(float)
    obs = aucs(X0)
    print(f"\n観測: 中心化なし {obs['raw']:.4f}  中心化あり {obs['centred']:.4f}  "
          f"差 {obs['raw'] - obs['centred']:.4f}")

    rng = np.random.default_rng(SEED)
    rec = {k: np.full(N_BOOT, np.nan) for k in ("raw", "centred", "drop")}
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
        X = np.column_stack([
            PU.delta(None, S[st], {"case": drawn[("case", st)], "ctrl": drawn[("ctrl", st)]})
            for st in states])
        r = aucs(X)
        rec["raw"][i], rec["centred"][i] = r["raw"], r["centred"]
        rec["drop"][i] = r["raw"] - r["centred"]
        if (i + 1) % 200 == 0:
            print(f"  {i + 1}/{N_BOOT}", flush=True)

    Rb = pd.DataFrame(rec).dropna()
    out = [
        ("replicates drawn", N_BOOT),
        ("replicates undefined", n_undef),
        ("replicates used", len(Rb)),
        ("observed AUC, uncentred pathway means", round(obs["raw"], 4)),
        ("observed AUC, centred pathway means", round(obs["centred"], 4)),
        ("observed drop", round(obs["raw"] - obs["centred"], 4)),
    ]
    for lab, col in (("uncentred", "raw"), ("centred", "centred"), ("drop", "drop")):
        out += [(f"bootstrap {lab} median", round(float(Rb[col].median()), 4)),
                (f"bootstrap {lab} 2.5th percentile", round(float(np.percentile(Rb[col], 2.5)), 4)),
                (f"bootstrap {lab} 97.5th percentile", round(float(np.percentile(Rb[col], 97.5)), 4))]
    out += [
        ("replicates with uncentred AUC above centred", round(float((Rb["drop"] > 0).mean()), 4)),
        ("replicates with drop at least 0.20", round(float((Rb["drop"] >= 0.20).mean()), 4)),
        ("replicates with centred AUC below 0.60", round(float((Rb["centred"] < 0.60).mean()), 4)),
        ("replicates with uncentred AUC above 0.80", round(float((Rb["raw"] > 0.80).mean()), 4)),
    ]
    T = pd.DataFrame(out, columns=["quantity", "value"])
    T.to_csv(OUT / "centering_auc_animal_bootstrap.tsv", sep="\t", index=False)
    Rb.round(5).to_csv(OUT / "centering_auc_animal_bootstrap_replicates.tsv", sep="\t", index=False)
    print("\n" + T.to_string(index=False))
    print(f"\n書き出し: {OUT / 'centering_auc_animal_bootstrap.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
