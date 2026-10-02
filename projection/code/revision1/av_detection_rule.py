# -*- coding: utf-8 -*-
"""検出判定を、二つのプロテオーム表の実装どおりに数える（AV1 の解析部分）。

実装を読んで確定した事実:
  Pod-TRECK（Supplementary Data S1 として添付した定量表）
    空欄は 0 個、定量されなかった検体は数値の 0。読み込み時に 0 を欠測に変換する
    （src/01_load.py の load_mouse_prot: df.replace(0, np.nan)）。
    9 検体の 50% 以上 → 5 検体以上で正の値があれば検出。
  ネコ皮質・髄質プロテオーム（Li らの補足 S5/S6）
    欠測は空欄（数値として読めないセル）。ちょうど 0 の値は皮質 5 件・髄質 4 件あり、
    これは値として残す。皮質 23 検体の 12 検体以上、髄質 19 検体の 10 検体以上。

本文への反映は apply_AV1.py が行う（原稿の文言を持つので公開しない）。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
M = HERE / "manuscript_C"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "code" / "revision1"))
import build_delta_matrix as BDM   # noqa: E402
import af7_attach_source_table as AF7  # noqa: E402
OUT = HERE / "results" / "roundR" / "detection_rule.tsv"




def facts() -> dict:
    import openpyxl
    ws = openpyxl.load_workbook(AF7.DEST / AF7.ATTACH, read_only=True,
                               data_only=True)["解析結果"]
    it = ws.iter_rows(values_only=True)
    next(it)
    cells = zeros = blanks = rows_n = dash_n = 0
    for r in it:
        if r[5] is None or not str(r[5]).strip():
            continue
        rows_n += 1
        if str(r[5]).strip() in BDM.NO_SYMBOL:
            dash_n += 1
        for i in range(10, 19):
            cells += 1
            v = r[i]
            if v is None or (isinstance(v, str) and not str(v).strip()):
                blanks += 1
            elif float(v) == 0:
                zeros += 1
    f = {"pod_treck cells": cells, "pod_treck blank cells": blanks,
         "pod_treck zero cells": zeros, "pod_treck data rows": rows_n,
         "pod_treck rows without a gene symbol": dash_n}
    for nm, lab in (("mouse_prot", "pod_treck"), ("cat_prot_ctx", "feline cortex"),
                    ("cat_prot_med", "feline medulla")):
        d = BDM.read_int(nm)
        # 記号が付与されていない行は §4.2 で除外しているので、集計からも除く
        d = d.loc[[g for g in d.index if str(g).strip() not in BDM.NO_SYMBOL]]
        f[f"{lab} samples"] = d.shape[1]
        f[f"{lab} threshold"] = int(np.ceil(d.shape[1] * 0.5))
        f[f"{lab} proteins"] = d.shape[0]
        f[f"{lab} detected"] = int((d.notna().mean(axis=1) >= 0.5).sum())
        f[f"{lab} exact zeros kept"] = int((d == 0).sum().sum())
        f[f"{lab} missing cells"] = int(d.isna().sum().sum())
    return f



def main() -> int:
    F = facts()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"quantity": k, "value": v} for k, v in F.items()]).to_csv(
        OUT, sep="\t", index=False)
    for k, v in F.items():
        print(f"  {k:32s} {v:,}")
    print(f"書き出し: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
