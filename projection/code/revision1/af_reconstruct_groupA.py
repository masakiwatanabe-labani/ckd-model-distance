# -*- coding: utf-8 -*-
"""AF1/AF2: 公開された Pod-TRECK プロテオーム Table S1 から検出コールを再構成する。

判定したいこと:
  1. 公開表にサンプルごとの欠測が分かる形で値が入っているか
  2. 入っているなら「9 サンプルの 50% 以上で非欠測」という現行の判定を再現できるか
  3. 再現できるなら、現行の Group A（2,240 / complete case 2,016）と一致するか

パイプライン側（build_delta_matrix.group_a_sets）の手順:
  元の定量表（受託解析の表。名前は公開しない。local_source.py 参照）
  → 遺伝子シンボルごとに nanmax → 0 を NaN → log2 → 列中央値正規化
  → 非欠測が 9 サンプルの 50% 以上 → ヒト空間へ写す
公開表は log2 intensity が既に入っているので、検出の判定（非欠測かどうか）だけを比べる。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
OUT = HERE / "results" / "roundR"
# 公開稿の Supplementary Table S1。手元の置き場所は環境変数で与える
# （既定はリポジトリの data/external/）。絶対パスを書かない。
import os
XLSX = Path(os.environ.get("CKD_PUBLISHED_TABLE_S1",
                           ROOT / "data" / "external" / "Supplementary_Table_S1.xlsx"))
import local_source as _LS  # noqa: E402
SAMPLES = [f"{_LS.sample_prefix()}-{i:02d}_{g}" for i, g in
           zip(range(1, 10), ["Ctrl-1", "Ctrl-2", "Ctrl-3", "Day14-1", "Day14-2", "Day14-3",
                              "Day21-1", "Day21-2", "Day21-3"])]
FRAC = 0.5


def load_published() -> pd.DataFrame:
    """公開表の 3 シートを結合し、遺伝子 x 9 サンプルの log2 強度行列にする。"""
    import openpyxl
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    per_gene: dict[str, dict[str, float]] = {}
    n_rows = {}
    for sheet in [s for s in wb.sheetnames if s != "README"]:
        ws = wb[sheet]
        it = ws.iter_rows(values_only=True)
        hdr = [str(c) if c is not None else "" for c in next(it)]
        cols = {h.split("|")[1].strip(): i for i, h in enumerate(hdr)
                if h.startswith("Log2 protein intensity")}
        gi = hdr.index("Gene Symbol")
        n = 0
        for r in it:
            if r is None or gi >= len(r):
                continue
            g = r[gi]
            if g is None or not str(g).strip():
                continue
            n += 1
            d = per_gene.setdefault(str(g).strip(), {})
            for s, i in cols.items():
                v = r[i]
                if v is None or (isinstance(v, str) and not v.strip()):
                    continue
                try:
                    x = float(v)
                except (TypeError, ValueError):
                    continue
                if not np.isfinite(x):
                    continue
                d[s] = max(d.get(s, -np.inf), x)      # パイプラインと同じく重複は最大値
        n_rows[sheet] = n
    M = pd.DataFrame.from_dict(per_gene, orient="index").reindex(columns=SAMPLES)
    return M, n_rows


def main() -> int:
    if not XLSX.exists():
        print(f"公開表が無い: {XLSX}", file=sys.stderr)
        return 2
    M, n_rows = load_published()
    print(f"公開表: シートの行数 {n_rows}")
    print(f"  遺伝子シンボル {M.shape[0]} / サンプル列 {M.shape[1]}")
    filled = M.notna().sum()
    print("  サンプルごとの非欠測数:")
    for s in SAMPLES:
        print(f"    {s:24s} {int(filled[s])}")
    miss = int(M.isna().sum().sum())
    print(f"  欠測セル {miss} / {M.size}（{miss / M.size:.1%}）")

    det_pub = set(M.index[M.notna().mean(axis=1) >= FRAC])
    print(f"  50% 以上で非欠測: {len(det_pub)} 遺伝子シンボル")
    npub = M.notna().sum(axis=1)
    print("  非欠測サンプル数の分布: " +
          ", ".join(f"{k}列={v}" for k, v in sorted(npub.value_counts().items())))

    # パイプラインが使った元表からの検出コール
    import build_delta_matrix as BDM
    src = BDM.read_int("mouse_prot")
    det_src = set(src.index[src.notna().mean(axis=1) >= FRAC])
    print(f"\n元表: {src.shape[0]} 遺伝子 x {src.shape[1]} サンプル")
    print(f"  50% 以上で非欠測: {len(det_src)} 遺伝子シンボル")
    nsrc = src.notna().sum(axis=1)
    print("  非欠測サンプル数の分布: " +
          ", ".join(f"{k}列={v}" for k, v in sorted(nsrc.value_counts().items())))
    print(f"  欠測を1つ以上もつ遺伝子: {int((nsrc < 9).sum())} / {len(nsrc)}")
    print(f"  元表にあって公開表に無い: {len(set(src.index) - set(M.index))}")
    print(f"  公開表にあって元表に無い: {len(set(M.index) - set(src.index))}")

    only_pub, only_src = det_pub - det_src, det_src - det_pub
    print(f"\nマウス側の検出コールの一致: 共通 {len(det_pub & det_src)} / "
          f"公開のみ {len(only_pub)} / 元表のみ {len(only_src)}")

    # ヒト空間へ写して Group A を作り直す
    def to_h(names):
        s = pd.Series(1.0, index=sorted(names))
        h, _ = BDM.to_human(s, "mouse")
        return set(h.index)

    ga_now = set((HERE / "groupA_intersection.txt").read_text().split())
    hp, hs = to_h(det_pub), to_h(det_src)
    print(f"\nヒト空間: 公開表 {len(hp)} / 元表 {len(hs)}")

    ctx = BDM.read_int("cat_prot_ctx")
    med = BDM.read_int("cat_prot_med")

    def cat_det(d, sp="cat"):
        s = pd.Series(1.0, index=d.index[d.notna().mean(axis=1) >= FRAC])
        h, _ = BDM.to_human(s, sp)
        return set(h.index)

    c, m = cat_det(ctx), cat_det(med)
    delta_idx = set(pd.read_csv(HERE / "delta_matrix.tsv", sep="\t",
                                index_col=0, usecols=[0]).index)
    ga_pub = sorted(((c & m) & hp) & delta_idx)
    ga_src = sorted(((c & m) & hs) & delta_idx)
    print(f"\nGroup A（cortex ∩ medulla ∩ mouse prot）")
    print(f"  現行のファイル      {len(ga_now)}")
    print(f"  元表から再計算      {len(ga_src)}  一致 {set(ga_src) == ga_now}")
    print(f"  公開表から再構成    {len(ga_pub)}  一致 {set(ga_pub) == ga_now}")
    d1, d2 = set(ga_pub) - ga_now, ga_now - set(ga_pub)
    print(f"  差分: 公開のみ {len(d1)} / 現行のみ {len(d2)}")
    if d1:
        print(f"    公開のみの例: {sorted(d1)[:10]}")
    if d2:
        print(f"    現行のみの例: {sorted(d2)[:10]}")

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"quantity": "published table, gene symbols", "value": M.shape[0]},
        {"quantity": "published table, sample columns", "value": M.shape[1]},
        {"quantity": "published table, missing cells", "value": miss},
        {"quantity": "published table, detected symbols", "value": len(det_pub)},
        {"quantity": "source table, gene symbols", "value": src.shape[0]},
        {"quantity": "source table, detected symbols", "value": len(det_src)},
        {"quantity": "source table, symbols with any missing sample",
         "value": int((nsrc < src.shape[1]).sum())},
        {"quantity": "one-gene difference", "value": sorted(d1)[0] if len(d1) == 1 else ""},
        {"quantity": "detected symbols in common", "value": len(det_pub & det_src)},
        {"quantity": "detected only in published", "value": len(only_pub)},
        {"quantity": "detected only in source", "value": len(only_src)},
        {"quantity": "Group A, current file", "value": len(ga_now)},
        {"quantity": "Group A, recomputed from source", "value": len(ga_src)},
        {"quantity": "Group A, reconstructed from published", "value": len(ga_pub)},
        {"quantity": "Group A, published only", "value": len(d1)},
        {"quantity": "Group A, current only", "value": len(d2)},
    ]).to_csv(OUT / "af_groupA_reconstruction.tsv", sep="\t", index=False)
    pd.Series(sorted(ga_pub)).to_csv(OUT / "af_groupA_from_published.txt",
                                     index=False, header=False)
    print(f"\n書き出し: {OUT / 'af_groupA_reconstruction.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
