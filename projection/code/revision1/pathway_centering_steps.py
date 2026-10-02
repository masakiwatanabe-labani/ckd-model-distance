# -*- coding: utf-8 -*-
"""§2.3 の四条件（中心化の有無 × 経路平均の有無）をペア単位で書き出す。

第10ラウンドで書き直したもの。この計算の出力である
  results/pathway/cos_before_after.tsv        （GO BP + KEGG + Hallmark、193 セット）
  results/pathway/aggregation_steps.tsv       （同上の中央値と四分位）
  results/revision1/pathway_reactome/cos_before_after.tsv （Reactome 込み 422 セット）
は、読むコード（verify_numbers.py、make_figures.py）はあるのに書くコードが
リポジトリにもその履歴にも無い状態だった。本文の 0.813 / 0.759 / 0.511 / 0.868 と
Figure 4D を供給しているのはこの3ファイルなので、入力から作り直せるようにする。

入力は Δ 行列・Group A の遺伝子リスト・ペア表・GMT だけで、既存の結果ファイルは読まない。

四条件（X は Δ 行列、列が状態、行が遺伝子）:
  1. cos_raw          cos θ の X そのまま
  2. cos_centred      各状態から全遺伝子平均を引いた C = X - mean_g(X) の cos θ
                      （= 生の Δ の Pearson 相関）
  3. cos_agg_centred  C を経路ごとに平均した 422（193）次元ベクトルの cos θ
  4. cos_agg_raw      X を経路ごとに平均した同上

標準化のうち SD で割る操作は列ごとの正のスカラー倍なので cos θ に効かない。
経路平均のあとも各状態の共通因子として残るため、そこでも消える（§4.13）。
したがって上の C からの計算は、z 化してから平均した場合と一致する。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
from lib_stats import load_config  # noqa: E402

# 参照ファイルの置き場は config から取る。build_delta_matrix を import すると
# オルソログ表まで読み込みにいくので、GMT の場所だけが要るここでは import しない。
REF = ROOT / load_config(ROOT / "config" / "config.yaml")["paths"]["ref"]
MIN_GENES = 30                      # pathway_cos.py と同じ。集合の最小遺伝子数

COLL_BASE = {"GO_BP": "geneset_gobp.gmt", "KEGG": "geneset_kegg.gmt",
             "Hallmark": "geneset_hallmark.gmt"}
COLL_FULL = dict(COLL_BASE, Reactome="geneset_reactome.gmt")

COND = [("cos_raw", "1. cos of Δ as analysed"),
        ("cos_centred", "2. Δ centred per state (= Pearson r of Δ)"),
        ("cos_agg_centred", "3. pathway means of centred Δ"),
        ("cos_agg_raw", "4. pathway means of uncentred Δ")]


def read_gmt(path):
    """GMT を {集合名: 遺伝子の集合} に読む（pathway_cos.py と同一の実装）。"""
    out = {}
    for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            out[f[0]] = {g.strip().upper() for g in f[2:] if g.strip()}
    return out


def load():
    """Δ 行列（Group A、全 16 状態で有限）と、ペア表。"""
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    P = pd.read_csv(HERE / "results" / "control_analysis" / "pairs.tsv", sep="\t")
    P["wi"] = (P.a.str.startswith("cat")) == (P.b.str.startswith("cat"))
    return D, P


def collect(D, colls):
    """各コレクションから、Group A 内で MIN_GENES 以上ある集合の遺伝子位置を取る。"""
    genes = set(D.index)
    gpos = {g: i for i, g in enumerate(D.index)}
    out = []
    for coll, fn in colls.items():
        for name, members in read_gmt(REF / fn).items():
            gi = np.array(sorted(gpos[g] for g in members & genes))
            if len(gi) >= MIN_GENES:
                out.append((coll, name, gi))
    return out


def conditions(D, P, paths):
    """対照を共有しない 87 ペアについて、四条件の cos θ を返す。"""
    X = D.to_numpy(float)
    C = X - X.mean(axis=0)                       # 状態ごとに全遺伝子平均を引く
    agg_c = np.vstack([C[gi].mean(axis=0) for _c, _n, gi in paths])
    agg_r = np.vstack([X[gi].mean(axis=0) for _c, _n, gi in paths])

    idx = {s: i for i, s in enumerate(D.columns)}
    ia, ib = P.a.map(idx).to_numpy(), P.b.map(idx).to_numpy()
    m = (~P.shared_control).to_numpy()

    def cosmat(M):
        Mn = M / np.linalg.norm(M, axis=0, keepdims=True)
        return (Mn.T @ Mn)[ia, ib]

    out = pd.DataFrame({"cos_raw": cosmat(X), "cos_centred": cosmat(C),
                        "cos_agg_centred": cosmat(agg_c), "cos_agg_raw": cosmat(agg_r),
                        "within": P.wi.astype(int).to_numpy()})
    return out.loc[m].reset_index(drop=True)


def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


def steps_table(BA):
    """条件 × クラスの件数・中央値・四分位。"""
    rows = []
    for col, label in COND:
        for cls, w in (("within", 1), ("cross", 0)):
            v = BA[BA.within == w][col]
            rows.append({"condition": label, "cls": cls, "n": len(v),
                         "median": v.median(), "q1": v.quantile(.25), "q3": v.quantile(.75)})
    return pd.DataFrame(rows)


def report(BA, label):
    print(f"  {label}: {len(BA)} ペア（within {int(BA.within.sum())} / "
          f"cross {int((1 - BA.within).sum())}）")
    for col, _lab in COND:
        a, b = BA[BA.within == 1][col], BA[BA.within == 0][col]
        print(f"    {col:16s} AUC {auc(a, b):.4f}   within {a.median():.4f} / "
              f"cross {b.median():.4f}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=None,
                    help="既定は results/ の本来の場所。突き合わせ用に別の場所へ出せる。")
    a = ap.parse_args()

    missing = [fn for fn in COLL_FULL.values() if not (REF / fn).exists()]
    if missing:
        print(f"遺伝子セットの GMT が {REF} に無い: {missing}\n"
              "再配布していないので、先に `python src/00_fetch_refs.py` で取得すること。",
              file=sys.stderr)
        return 2

    D, P = load()
    print(f"Δ 行列: {D.shape[0]} 遺伝子 x {D.shape[1]} 状態 / ペア表 {len(P)} 組")

    base, full = collect(D, COLL_BASE), collect(D, COLL_FULL)
    print(f"セット数: 既存 3 コレクション {len(base)} / Reactome 込み {len(full)}")

    BA_base, BA_full = conditions(D, P, base), conditions(D, P, full)
    report(BA_base, f"GO+KEGG+Hallmark ({len(base)})")
    report(BA_full, f"+ Reactome ({len(full)})")

    if a.out_dir:
        d1 = d2 = Path(a.out_dir)
        (d1 / "pathway").mkdir(parents=True, exist_ok=True)
        (d1 / "pathway_reactome").mkdir(parents=True, exist_ok=True)
        p_ba, p_st = d1 / "pathway" / "cos_before_after.tsv", d1 / "pathway" / "aggregation_steps.tsv"
        p_ba2 = d2 / "pathway_reactome" / "cos_before_after.tsv"
    else:
        p_ba = HERE / "results" / "pathway" / "cos_before_after.tsv"
        p_st = HERE / "results" / "pathway" / "aggregation_steps.tsv"
        p_ba2 = HERE / "results" / "revision1" / "pathway_reactome" / "cos_before_after.tsv"
        for p in (p_ba, p_st, p_ba2):
            p.parent.mkdir(parents=True, exist_ok=True)

    BA_base.round(4).to_csv(p_ba, sep="\t", index=False)
    steps_table(BA_base).round(4).to_csv(p_st, sep="\t", index=False)
    BA_full.round(4).to_csv(p_ba2, sep="\t", index=False)
    for p in (p_ba, p_st, p_ba2):
        print(f"書き出し: {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
