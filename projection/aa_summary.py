"""AA. 遺伝子空間を替えても主要結論が再現するかを 1 表にまとめる。

比較する空間:
  groupA        Group A intersection（本文の主解析）
  ortholog_all  Ensembl 1:1 orthologue 全体（Group A 制限なし）
  groupB_matched 対照群発現量を Group A にマッチさせた Group B

各空間で §2.1（alpha の反例）、§2.3（時間）、§2.4（cat-mouse と天井）、
§2.5（経路単位 vs 集約）の 4 点を再計算した結果を読む。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
R = HERE / "results"
OUT = R / "aa_genespace"
OUT.mkdir(parents=True, exist_ok=True)

SPACES = [("groupA", "", "alpha_groupA", "Group A intersection（主解析）"),
          ("ortholog_all", "_ortholog_all", "alpha_ortholog_all", "全 1:1 orthologue"),
          ("groupB_matched", "_groupB_matched", "alpha_groupB_matched", "発現量マッチ Group B")]


def row(space, tag, alpha_dir, label):
    r = {"空間": label}
    pr = pd.read_csv(R / alpha_dir / "projection.tsv", sep="\t", index_col=0)
    r["n_genes"] = int(pr.n_genes.iloc[0])
    m = pr.loc["cat_med_CKD12"]
    r["猫髄質CKD1/2 alpha"] = round(float(m.alpha), 3)
    r["同 cos"] = round(float(m["cos"]), 3)
    r["同 nrm"] = round(float(m.nrm), 3)
    r["alpha>1.000"] = bool(m.alpha > 1.0)

    ta = pd.read_csv(R / alpha_dir / "time_association.tsv", sep="\t")
    ta = ta[(ta.subset == "all_model_states") & (ta.term == "cos x log_time")].iloc[0]
    r["cos x log_time"] = round(float(ta.rho), 3)
    r["同 CI"] = f"[{ta.boot_lo:+.2f}, {ta.boot_hi:+.2f}]"

    mt = pd.read_csv(R / f"mantel_s1{tag}" / "mantel_s1.tsv", sep="\t")
    mt = mt[(mt.subset.str.startswith("12")) & (mt.dissimilarity.str.startswith("D_cos"))
            & (mt.predictor == "elapsed time |Δlog10 h|")].iloc[0]
    r["Mantel 時間 D_cos"] = round(float(mt.mantel_rho), 3)
    r["Mantel p"] = round(float(mt.p_perm), 4)

    cc = pd.read_csv(R / f"reliability{tag}" / "ceiling_comparison.tsv", sep="\t").set_index("class")
    cs = cc.loc["cross_species"]
    r["cat-mouse cos 中央"] = round(float(cs.cos_median), 3)
    r["cat-mouse cos 最大"] = round(float(cs.cos_max), 3)
    r["天井 raw 最小"] = round(float(cs.ceiling_raw_min), 3)
    r["天井超え(raw)"] = int(cs.n_exceed_ceiling_raw)

    sp = pd.read_csv(R / f"pathway{tag}" / "separation_summary.tsv", sep="\t")
    sp = sp[sp.restriction == "exclude shared controls"].iloc[0]
    r["経路単位 AUC 中央"] = round(float(sp.pathway_auc_median), 3)
    r["集約後 AUC"] = round(float(sp.aggregated_auc), 3)

    ap = pd.read_csv(R / f"control_analysis{tag}" / "auc_permutation.tsv", sep="\t")
    ap = ap[ap.level.str.contains("状態レベル")].iloc[0]
    r["全パネル AUC"] = round(float(ap.auc_obs), 4)
    r["AUC p"] = round(float(ap.p), 4)
    return r


def main():
    rows = [row(*s) for s in SPACES]
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "genespace_comparison.tsv", sep="\t", index=False)
    print(T.to_string(index=False))
    print()
    a = T.set_index("空間")
    base = a.index[0]
    print("== 主解析に対する再現の判定 ==")
    for k in ["alpha>1.000"]:
        print(f"  {k}: " + " / ".join(f"{i}={a.loc[i,k]}" for i in a.index))
    for k in ["cos x log_time", "Mantel 時間 D_cos", "cat-mouse cos 最大", "天井 raw 最小",
              "経路単位 AUC 中央", "集約後 AUC"]:
        print(f"  {k}: " + " / ".join(f"{i}={a.loc[i,k]}" for i in a.index))


if __name__ == "__main__":
    main()
