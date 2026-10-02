# -*- coding: utf-8 -*-
"""再生成したファイルと、リポジトリに置かれていた既存ファイルをセル単位で突き合わせる。

既存ファイルは答え合わせにのみ使う。差が出たら差として報告する（実装を寄せない）。

  python compare_round10.py <既存ファイルを相対パスのまま置いたディレクトリ>
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
NEW = HERE / "results" / "round10"
PAIRS = [
    ("results/pathway/cos_before_after.tsv", NEW / "pathway" / "cos_before_after.tsv"),
    ("results/pathway/aggregation_steps.tsv", NEW / "pathway" / "aggregation_steps.tsv"),
    ("results/revision1/pathway_reactome/cos_before_after.tsv",
     NEW / "pathway_reactome" / "cos_before_after.tsv"),
]


def main() -> int:
    if len(sys.argv) < 2:
        print("既存ファイルのあるディレクトリを渡すこと", file=sys.stderr)
        return 2
    old_root = Path(sys.argv[1])
    rows, bad = [], 0
    for rel, new_p in PAIRS:
        old_p = old_root / rel
        if not old_p.exists():
            print(f"見つからない: {old_p}", file=sys.stderr)
            return 2
        o, n = pd.read_csv(old_p, sep="\t"), pd.read_csv(new_p, sep="\t")
        same_shape, same_cols = o.shape == n.shape, list(o.columns) == list(n.columns)
        if same_shape and same_cols:
            num = list(o.select_dtypes(include=[np.number]).columns)
            d = o[num].to_numpy(float) - n[num].to_numpy(float)
            maxabs = float(np.nanmax(np.abs(d))) if d.size else 0.0
            ncell, ndiff = int(d.size), int(np.sum(np.abs(d) > 0))
            same_txt = all((o[c] == n[c]).all() for c in o.columns if c not in num)
        else:
            maxabs, ncell, ndiff, same_txt = float("nan"), 0, -1, False
        ok = same_shape and same_cols and same_txt and ndiff == 0
        bad += 0 if ok else 1
        rows.append({"file": rel, "shape_old": str(o.shape), "shape_new": str(n.shape),
                     "columns_identical": same_cols, "text_columns_identical": same_txt,
                     "numeric_cells": ncell, "cells_differing": ndiff,
                     "max_abs_difference": maxabs, "identical": ok})
        print(f"{'一致 ' if ok else '相違 '} {rel}: {o.shape} vs {n.shape}, "
              f"数値セル {ncell}、相違 {ndiff}、最大差 {maxabs:.2e}")
    NEW.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(NEW / "regeneration_comparison.tsv", sep="\t", index=False)
    print(f"\n書き出し: {NEW / 'regeneration_comparison.tsv'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
