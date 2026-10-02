# -*- coding: utf-8 -*-
"""AF2: 公開表と、パイプラインが使った元表の値そのものを比べる。

Group A が 1 遺伝子ずれた理由を、表の性質として特定する。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE), )
sys.path.insert(0, str(HERE / "code" / "revision1"))
import build_delta_matrix as BDM          # noqa: E402
import af_reconstruct_groupA as AF        # noqa: E402

OUT = HERE / "results" / "roundR"


def main() -> int:
    M, _ = AF.load_published()
    src = BDM.read_int("mouse_prot")
    common = sorted(set(M.index) & set(src.index))
    A = src.loc[common, AF.SAMPLES].to_numpy(float)
    B = M.loc[common, AF.SAMPLES].to_numpy(float)
    both = np.isfinite(A) & np.isfinite(B)

    d = B[both] - A[both]
    r = float(np.corrcoef(A[both], B[both])[0, 1])
    # 元表で欠測、公開表で値があるセル
    only_pub = np.isnan(A) & np.isfinite(B)
    # そのセルの公開値が、その遺伝子の観測値の中で低い側か高い側か
    lo = []
    for i, j in zip(*np.where(only_pub)):
        row = B[i][np.isfinite(B[i])]
        lo.append(float((row < B[i, j]).mean()))
    rows = [
        {"quantity": "gene symbols in both tables", "value": len(common)},
        {"quantity": "cells with a value in both", "value": int(both.sum())},
        {"quantity": "Pearson r of shared cells", "value": r},
        {"quantity": "median published minus source", "value": float(np.median(d))},
        {"quantity": "IQR of the difference", "value": float(np.subtract(*np.percentile(d, [75, 25])))},
        {"quantity": "cells missing in source but present in published",
         "value": int(only_pub.sum())},
        {"quantity": "median within-gene percentile of those published values",
         "value": float(np.median(lo)) if lo else np.nan},
        {"quantity": "cells present in source but missing in published",
         "value": int((np.isfinite(A) & np.isnan(B)).sum())},
    ]
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "af_published_vs_source.tsv", sep="\t", index=False)
    for r_ in rows:
        v = r_["value"]
        print(f"  {r_['quantity']:56s} {v:,.4f}" if isinstance(v, float)
              else f"  {r_['quantity']:56s} {v:,}")
    print(f"\n書き出し: {OUT / 'af_published_vs_source.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
