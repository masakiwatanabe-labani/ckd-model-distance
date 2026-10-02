# -*- coding: utf-8 -*-
"""どのマウス状態がいくつの経路で首位に立つかを書き出す（§2.3、Figure 4A）。

第10ラウンドで書き直したもの。出力の
  results/revision1/leads_by_state_breakdown.tsv
  results/revision1/leads_IRI_2h_top10.tsv
は、本文の「IRI 28 d 22%、IRI 12 mo 21%、IRI 2 h 11%」「48 のうち 47 が Reactome」
「上位 10 は 30〜33 遺伝子のプロテアソーム・ユビキチン分解系」を供給しているのに、
書くコードがリポジトリにもその履歴にも無かった。

入力は results/revision1/pathway_reactome/pathway_cos.tsv（pathway_cos.py が作る、
422 セット x 16 状態の cos θ）だけ。首位はマウス 12 状態の中でとる。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
SRC = HERE / "results" / "revision1" / "pathway_reactome" / "pathway_cos.tsv"
KEY = ["collection", "pathway"]


def leads() -> pd.DataFrame:
    d = pd.read_csv(SRC, sep="\t")
    mouse = d[~d.state.str.startswith("cat")]
    return mouse.loc[mouse.groupby(KEY)["cos"].idxmax()]


def breakdown(lead: pd.DataFrame) -> pd.DataFrame:
    n = len(lead)
    base = (lead.groupby("state").size().rename("n_leads").reset_index()
            .assign(share=lambda t: t.n_leads / n))
    per_coll = (lead.pivot_table(index="state", columns="collection", values="cos",
                                 aggfunc="size").fillna(0).astype(int).reset_index())
    cols = [c for c in ["Reactome", "GO_BP", "KEGG", "Hallmark"] if c in per_coll.columns]
    out = base.merge(per_coll[["state"] + cols], on="state", how="left")
    return out.sort_values(["n_leads", "state"], ascending=[False, True]).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--state", default="IRI_2h", help="上位 10 を出す状態")
    a = ap.parse_args()

    lead = leads()
    br = breakdown(lead)
    top = (lead[lead.state == a.state].nlargest(10, "cos")[KEY + ["n_genes", "cos"]]
           .reset_index(drop=True))

    out = Path(a.out_dir) if a.out_dir else (HERE / "results" / "revision1")
    out.mkdir(parents=True, exist_ok=True)
    br.round(4).to_csv(out / "leads_by_state_breakdown.tsv", sep="\t", index=False)
    top.round(4).to_csv(out / "leads_IRI_2h_top10.tsv", sep="\t", index=False)

    print(f"経路 {len(lead)} 件で首位を判定（マウス {lead.state.nunique()} 状態）")
    for r in br.head(3).itertuples():
        reac = getattr(r, "Reactome", 0)
        print(f"  {r.state:9s} {r.n_leads:3d} 経路 ({r.share:.1%})  うち Reactome {reac}")
    print(f"  {a.state} の上位 10: 遺伝子数 {top.n_genes.min()}–{top.n_genes.max()}, "
          f"コレクション {sorted(top.collection.unique())}")
    for p in ("leads_by_state_breakdown.tsv", "leads_IRI_2h_top10.tsv"):
        print(f"書き出し: {out / p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
