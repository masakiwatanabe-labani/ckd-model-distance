# -*- coding: utf-8 -*-
"""AF7: Group A の検出コールに用いた定量表を supplementary として添付する。

ファイルは加工しない。元のバイト列のまま複製し、チェックサムを記録する。
添付したファイルから検出コールを計算し直し、Group A（2,240 / complete case 2,016）が
再現することを確かめる。ここが合わなくなれば、添付ファイルと本文が食い違っている。
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import build_delta_matrix as BDM  # noqa: E402

OUT = HERE / "results" / "roundR"
DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
ATTACH = "Supplementary_Data_S1_PodTRECK_proteome.xlsx"
FRAC = 0.5


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def matrix_from_workbook(p: Path) -> pd.DataFrame:
    """src/01_load.py の load_mouse_prot と同じ手順で 遺伝子 x 9 サンプルを作る。"""
    ws = openpyxl.load_workbook(p, read_only=True, data_only=True)["解析結果"]
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    samples = hdr[10:19]
    acc: dict[str, list] = {}
    for r in it:
        if r[5] is None or not str(r[5]).strip():
            continue
        vals = []
        for i in range(10, 19):
            try:
                vals.append(float(r[i]))
            except (TypeError, ValueError):
                vals.append(np.nan)
        acc.setdefault(str(r[5]).strip(), []).append(vals)
    df = pd.DataFrame.from_dict(
        {g: list(np.nanmax(np.array(v), axis=0)) for g, v in acc.items()},
        orient="index", columns=samples)
    df = df.replace(0, np.nan)
    log = np.log2(df)
    return log - log.median() + log.median().mean()      # 列中央値正規化


def main() -> int:
    cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    src = ROOT / cfg["paths"]["raw"] / cfg["files"]["mouse_prot"]
    if not src.exists():
        print(f"元表が無い: {src}", file=sys.stderr)
        return 2
    # 添付ファイルは ak_strip_docprops.py が作る（著者でない個人名と受託先の内部パスを除く）。
    # ここでは作らない。あるものを検証するだけ。
    dst = DEST / ATTACH
    if not dst.exists():
        print(f"添付が無い: {dst}\n  先に code/revision1/ak_strip_docprops.py を実行する",
              file=sys.stderr)
        return 2
    h_dst = sha256(dst)
    print(f"添付: {dst.name}  {dst.stat().st_size:,} bytes")
    print(f"  sha256 {h_dst}")
    print(f"  元ファイル {src.name}（{src.stat().st_size:,} bytes、"
          f"sha256 {sha256(src)}）は投稿フォルダの外に置いたまま")

    # 添付したファイルから検出コールを作り直す
    M = matrix_from_workbook(dst)
    ref = BDM.read_int("mouse_prot")
    common = sorted(set(M.index) & set(ref.index))
    A = M.loc[common, ref.columns].to_numpy(float)
    B = ref.loc[common, ref.columns].to_numpy(float)
    same_shape = M.shape == ref.shape and set(M.index) == set(ref.index)
    both = np.isfinite(A) & np.isfinite(B)
    maxdiff = float(np.abs(A[both] - B[both]).max()) if both.any() else np.nan
    same_na = bool((np.isnan(A) == np.isnan(B)).all())
    print(f"\n添付ファイルから再計算した行列: {M.shape}  パイプラインの行列: {ref.shape}")
    print(f"  同じ遺伝子集合 {same_shape} / 欠測の位置が同一 {same_na} / 値の最大差 {maxdiff:.3e}")

    det = set(M.index[M.notna().mean(axis=1) >= FRAC])
    s = pd.Series(1.0, index=sorted(det))
    hm, _ = BDM.to_human(s, "mouse")

    def cat_det(name):
        d = BDM.read_int(name)
        t = pd.Series(1.0, index=d.index[d.notna().mean(axis=1) >= FRAC])
        h, _ = BDM.to_human(t, "cat")
        return set(h.index)

    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    ga = sorted(((cat_det("cat_prot_ctx") & cat_det("cat_prot_med")) & set(hm.index))
                & set(D.index))
    complete = set(D.dropna(axis=0, how="any").index)
    ga_c = [g for g in ga if g in complete]
    now = set((HERE / "groupA_intersection.txt").read_text().split())
    print(f"\n添付ファイルから再構成した Group A: {len(ga)}（うち 16 状態で有限 {len(ga_c)}）")
    print(f"  現行の groupA_intersection.txt と一致: {set(ga) == now}")

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"quantity": "attached file name", "value": ATTACH},
        {"quantity": "attached file bytes", "value": dst.stat().st_size},
        {"quantity": "attached file sha256", "value": h_dst},
        {"quantity": "gene symbols in the attached table", "value": M.shape[0]},
        {"quantity": "sample columns", "value": M.shape[1]},
        {"quantity": "matrix matches the pipeline input", "value": bool(same_shape and same_na)},
        {"quantity": "largest value difference", "value": maxdiff},
        {"quantity": "detected mouse symbols", "value": len(det)},
        {"quantity": "Group A from the attached table", "value": len(ga)},
        {"quantity": "Group A complete in 16 states", "value": len(ga_c)},
        {"quantity": "matches groupA_intersection.txt", "value": bool(set(ga) == now)},
    ]).to_csv(OUT / "af7_attached_table.tsv", sep="\t", index=False)
    print(f"\n書き出し: {OUT / 'af7_attached_table.tsv'}")
    # 元の表の名前とチェックサムは受託先の案件番号を含むので公開しない。
    # 添付（Supplementary Data S1）との対応だけ、公開しない場所に残す。
    loc = OUT / "local"
    loc.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"quantity": "original file name", "value": src.name},
        {"quantity": "original file sha256", "value": sha256(src)},
        {"quantity": "attached file name", "value": ATTACH},
        {"quantity": "attached file sha256", "value": h_dst},
    ]).to_csv(loc / "attachment_source.tsv", sep="\t", index=False)
    print(f"  対応（非公開）: {loc / 'attachment_source.tsv'}")
    return 0 if set(ga) == now else 1


if __name__ == "__main__":
    sys.exit(main())
