# -*- coding: utf-8 -*-
"""Part 1: 各状態の「遺伝子 universe 全体にわたる平均Δ」を出す。

2.1 の分解の段落は「state の遺伝子平均という共通成分が 422 座標すべてに同符号で入り、
その水準が猫とマウスで違う」と述べているが、その共通成分の値が本文にも表にも無かった。
Δ 行列から直接計算できるので、16 状態分を 1 列にして出す。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "code" / "revision1"))
import pathway_centering_steps as PCS  # noqa: E402

OUT = HERE / "results" / "roundF"
LABEL = {"cat_CKD12": "Feline cortex, CKD1/2", "cat_CKD34": "Feline cortex, CKD3/4",
         "cat_med_CKD12": "Feline medulla, CKD1/2", "cat_med_CKD34": "Feline medulla, CKD3/4",
         "mouse_5D": "Pod-TRECK day 5", "mouse_2W": "Pod-TRECK day 14",
         "mouse_3W": "Pod-TRECK day 21", "IRI_2h": "IRI 2 h", "IRI_4h": "IRI 4 h",
         "IRI_24h": "IRI 24 h", "IRI_48h": "IRI 48 h", "IRI_72h": "IRI 72 h",
         "IRI_7d": "IRI 7 d", "IRI_14d": "IRI 14 d", "IRI_28d": "IRI 28 d",
         "IRI_12mo": "IRI 12 mo"}


def main() -> int:
    D, _P = PCS.load()
    X = D.to_numpy(float)
    rows = []
    for s in D.columns:
        v = X[:, list(D.columns).index(s)]
        rows.append({"state": LABEL.get(s, s),
                     "species": "feline" if s.startswith("cat") else "mouse",
                     "mean_delta": float(v.mean()), "sd_delta": float(v.std(ddof=1)),
                     "median_delta": float(np.median(v)),
                     "frac_positive": float((v > 0).mean())})
    T = pd.DataFrame(rows).sort_values(["species", "mean_delta"], ascending=[True, False])
    OUT.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(OUT / "state_mean_shift.tsv", sep="\t", index=False)

    cat = T[T.species == "feline"].mean_delta
    mo = T[T.species == "mouse"].mean_delta
    sep = cat.max() < mo.min()
    print(f"Group A {X.shape[0]} 遺伝子 x {X.shape[1]} 状態")
    for r in T.itertuples():
        print(f"  {r.species:7s} {r.state:24s} {r.mean_delta:+.4f}")
    print(f"\n猫 4 状態:     {cat.min():+.4f} 〜 {cat.max():+.4f}")
    print(f"マウス 12 状態: {mo.min():+.4f} 〜 {mo.max():+.4f}")
    print(f"水準は分かれているか: {'はい（重なりなし）' if sep else 'いいえ（重なる）'}"
          f"  余白 {mo.min() - cat.max():+.4f}")
    print(f"\n書き出し: {OUT / 'state_mean_shift.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
