# -*- coding: utf-8 -*-
"""AA: GSE98622 の単位の判定根拠を結果ファイルに残す（Methods §4.1 が引く値の出典）。

これまで報告書にしか無かったので、機械照合の対象にできるよう書き出す。
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
OUT = HERE / "results" / "roundR"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    src = ROOT / "data" / "external" / "GSE98622" / "GSE98622_mouse-iri-master.xlsx"
    if not src.exists():
        print(f"IRI ファイルが無い: {src}", file=sys.stderr)
        return 2
    ws = openpyxl.load_workbook(src, read_only=True, data_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    # 数値列を特定する（最初の数列は注釈）
    rows, ncol = [], len(hdr)
    for r in it:
        rows.append(r)
    A = pd.DataFrame(rows, columns=[str(h) for h in hdr])
    num = A.apply(pd.to_numeric, errors="coerce")
    keep = [c for c in num.columns if num[c].notna().mean() > 0.9]
    X = num[keep].to_numpy(float)
    X = X[np.isfinite(X).all(axis=1)]
    nz = X[X > 0]
    # FPKM は小数を持つ。整数判定は厳密一致で行う（isclose では丸めが入る）
    frac_int = float(np.mean(nz == np.round(nz)))
    colsum = X.sum(axis=0)
    rec = [
        {"quantity": "file name", "value": src.name},
        {"quantity": "sha256", "value": sha256(src)},
        {"quantity": "rows in the sheet", "value": len(A)},
        {"quantity": "numeric sample columns", "value": len(keep)},
        {"quantity": "genes with finite values in every column", "value": X.shape[0]},
        {"quantity": "share of non-zero entries that are integers", "value": frac_int},
        {"quantity": "share of non-zero entries that are non-integer", "value": 1 - frac_int},
        {"quantity": "smallest non-zero value", "value": float(nz.min())},
        {"quantity": "smallest column sum", "value": float(colsum.min())},
        {"quantity": "largest column sum", "value": float(colsum.max())},
        {"quantity": "coefficient of variation of column sums",
         "value": float(colsum.std(ddof=1) / colsum.mean())},
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rec).to_csv(OUT / "aa_unit_evidence.tsv", sep="\t", index=False)
    for r in rec:
        v = r["value"]
        print(f"  {r['quantity']:48s} {v if isinstance(v, str) else f'{v:,.6g}'}")
    print(f"\n書き出し: {OUT / 'aa_unit_evidence.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
