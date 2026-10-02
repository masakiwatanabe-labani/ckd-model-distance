# -*- coding: utf-8 -*-
"""AI: 参照実装（R）に渡す入力と、比較対象になる Python 側のスコア行列を書き出す。

R 側には AUC を計算させない。スコア行列だけを作らせ、比較と AUC は Python の
既存の関数で行う。実装の違いだけを見たいので、経路集合は Python が選んだ 422 を
そのまま GMT にして渡す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402
import pathway_scoring_methods as PSM  # noqa: E402

OUT = HERE / "results" / "roundR" / "r_input"


def main() -> int:
    D, P = PCS.load()
    paths = PCS.collect(D, PCS.COLL_FULL)
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)
    genes = list(D.index)
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. Δ 行列と、状態ごとに中心化した行列
    pd.DataFrame(X, index=D.index, columns=D.columns).to_csv(
        OUT / "delta_uncentered.tsv", sep="\t", float_format="%.10g")
    pd.DataFrame(C, index=D.index, columns=D.columns).to_csv(
        OUT / "delta_centered.tsv", sep="\t", float_format="%.10g")

    # 2. 422 集合をそのまま GMT に。集合名は Python 側の (collection, name) を保つ
    with open(OUT / "pathways_422.gmt", "w") as fh:
        for c, n, gi in paths:
            name = f"{c}__{n}".replace("\t", " ")
            fh.write(name + "\tna\t" + "\t".join(genes[i] for i in gi) + "\n")
    names = [f"{c}__{n}" for c, n, _g in paths]

    # 3. 比較対象になる Python 側のスコア行列
    for tag, M in (("uncentered", X), ("centered", C)):
        groups = [gi for _c, _n, gi in paths]
        pd.DataFrame(PSM.score_singscore(M, groups), index=names,
                     columns=D.columns).to_csv(
            OUT / f"py_singscore_{tag}.tsv", sep="\t", float_format="%.12g")
        pd.DataFrame(PSM.score_plage(M, groups), index=names,
                     columns=D.columns).to_csv(
            OUT / f"py_plage_{tag}.tsv", sep="\t", float_format="%.12g")

    print(f"書き出し: {OUT}")
    print(f"  delta_uncentered.tsv / delta_centered.tsv  {X.shape[0]} x {X.shape[1]}")
    print(f"  pathways_422.gmt  {len(paths)} 集合（最小 {min(len(g) for _c,_n,g in paths)} 遺伝子）")
    print(f"  py_singscore_*.tsv / py_plage_*.tsv  {len(names)} x {X.shape[1]}")
    assert len(set(names)) == len(names), "集合名が一意でない。GMT が壊れる"
    return 0


if __name__ == "__main__":
    sys.exit(main())
