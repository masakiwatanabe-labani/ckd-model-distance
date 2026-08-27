"""候補遺伝子のヒト側（KPMP regional proteomics）での方向一致検証。

距離行列と違い個別分子の話なので、区画解像度の非対称性の影響を受けない。

注意:
  KPMP の CKD は糖尿病性腎症・高血圧性腎硬化症が主体で、ネコ自然発症CKDとは
  病因が異なる。一致しなかった場合は「否定された」ではなく
  「病因の異なるヒトCKDでは確認できなかった」と記述すること。

検定:
  個別遺伝子の有意性ではなく「方向一致そのもの」を所見として扱う。
    - 3種一致率  : 帰無 P=1/4（独立な符号なら 2*(1/2)^3）
    - ペアワイズ  : 帰無 P=1/2
  遺伝子は共発現により独立でないため、二項検定の p は反保守的（楽観的）。
  効果量（一致率）と併記して読むこと。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary, alias_lookup  # noqa: E402

log = get_logger("10_human_validation")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
EXT = ROOT / "data" / "external"

WATCHLIST = ["FAM3D", "PTN", "ANXA3", "ANGPTL2", "NPNT", "THBS1", "SULF1", "ITGB6",
             "SPP1", "MGP", "FBLN1", "SLC14A2", "AQP2", "NR3C2", "FKBP5", "SGK1",
             "TSC22D3", "SCNN1B", "SCNN1G", "NEDD4L", "HSD11B2", "FOXJ1", "DNAH5",
             "HAVCR1", "LCN2"]

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))
CAT2ALIAS = dict(zip(OMAP.cat_symbol, OMAP.mouse_aliases))
CAT2HUMAN = None  # ネコ→ヒトはシンボル一致で足りる（KPMP は HGNC シンボル）


def kpmp_table(comp: str) -> pd.DataFrame:
    p = EXT / "KPMP" / "DataLake_DEPs.txt"
    df = pd.read_csv(p, sep="\t")
    d = df[df.Comparison == f"CKD.vs.HRT.in.{comp}"].copy()
    gn = d.Gene_name.astype(str).str.strip()
    keep = ~(d.Gene_name.isna() | (gn == "") | (gn.str.lower() == "nan"))
    d, gn = d[keep], gn[keep]
    d = d.assign(_s=gn.str.upper(), _a=d.LogFC.abs())
    d = d.sort_values(["_s", "Adj_pvalue", "_a"], ascending=[True, True, False])
    d = d.drop_duplicates("_s")
    return d.set_index("_s")[["LogFC", "Adj_pvalue"]]


def build_table(genes, cat_tab) -> pd.DataFrame:
    TI, G = kpmp_table("TI"), kpmp_table("G")
    mp21 = DE.loc["m_prot_d21"]
    rows = []
    for g in genes:
        rec = {"gene": g}
        if g in cat_tab.index:
            rec["cat_ctx_late_lfc"] = round(float(cat_tab.loc[g, "lfc"]), 3)
            rec["cat_ctx_late_q"] = float(cat_tab.loc[g, "q"])
        else:
            rec["cat_ctx_late_lfc"] = np.nan
            rec["cat_ctx_late_q"] = np.nan
        hit = alias_lookup(mp21.index, CAT2ALIAS.get(g), CAT2MOUSE.get(g))
        if hit is None:
            hit = alias_lookup(mp21.index, None, g)
        if hit is not None:
            rec["mouse_symbol_used"] = hit
            rec["mouse_prot_d21_lfc"] = round(float(mp21.loc[hit, "lfc"]), 3)
            rec["mouse_prot_d21_q"] = float(mp21.loc[hit, "q"])
        else:
            rec["mouse_symbol_used"] = "未検出"
            rec["mouse_prot_d21_lfc"] = np.nan
            rec["mouse_prot_d21_q"] = np.nan
        for tag, tab in [("TI", TI), ("G", G)]:
            if g in tab.index:
                rec[f"kpmp_{tag}_logFC"] = round(float(tab.loc[g, "LogFC"]), 3)
                rec[f"kpmp_{tag}_adjP"] = float(tab.loc[g, "Adj_pvalue"])
            else:
                rec[f"kpmp_{tag}_logFC"] = np.nan
                rec[f"kpmp_{tag}_adjP"] = np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)

    def sgn(v):
        return None if pd.isna(v) or v == 0 else ("+" if v > 0 else "-")

    def verdict(r):
        s = [sgn(r.cat_ctx_late_lfc), sgn(r.mouse_prot_d21_lfc), sgn(r.kpmp_TI_logFC)]
        if any(x is None for x in s):
            miss = [n for n, x in zip(["ネコ", "マウス", "ヒトTI"], s) if x is None]
            return f"判定不可（{'/'.join(miss)}で未検出）"
        return f"3種一致（{s[0]}）" if len(set(s)) == 1 else \
               f"不一致（ネコ{s[0]} / マウス{s[1]} / ヒトTI{s[2]}）"

    t["direction_verdict"] = t.apply(verdict, axis=1)
    t["kpmp_TI_significant"] = t.kpmp_TI_adjP.apply(lambda v: bool(v < 0.05) if pd.notna(v) else None)
    t["kpmp_G_significant"] = t.kpmp_G_adjP.apply(lambda v: bool(v < 0.05) if pd.notna(v) else None)
    return t


def sign_tests(t: pd.DataFrame, label: str) -> dict:
    """方向一致そのものを検定する。"""
    out = {"集合": label}
    d = t.dropna(subset=["kpmp_TI_logFC"])
    out["ヒトTIで検出"] = f"{len(d)} / {len(t)}"

    # 3種一致（ネコRNA・マウス蛋白・ヒトTI蛋白がすべて同符号）
    tri = d.dropna(subset=["cat_ctx_late_lfc", "mouse_prot_d21_lfc"])
    if len(tri):
        s = np.sign(tri[["cat_ctx_late_lfc", "mouse_prot_d21_lfc", "kpmp_TI_logFC"]].values)
        agree = int((np.abs(s.sum(axis=1)) == 3).sum())
        p = stats.binomtest(agree, len(tri), 0.25, alternative="greater").pvalue
        out["3種一致"] = (f"{agree}/{len(tri)} = {agree/len(tri):.1%} "
                          f"(帰無25%, 二項検定 p={p:.4g})")

    # ペアワイズ（帰無 50%）
    pairs = [("cat_ctx_late_lfc", "kpmp_TI_logFC", "ネコRNA × ヒトTI蛋白"),
             ("mouse_prot_d21_lfc", "kpmp_TI_logFC", "マウス蛋白 × ヒトTI蛋白"),
             ("cat_ctx_late_lfc", "mouse_prot_d21_lfc", "ネコRNA × マウス蛋白")]
    for a, b, name in pairs:
        dd = d.dropna(subset=[a, b])
        if not len(dd):
            continue
        ag = int((np.sign(dd[a]) == np.sign(dd[b])).sum())
        p = stats.binomtest(ag, len(dd), 0.5, alternative="greater").pvalue
        out[name] = f"{ag}/{len(dd)} = {ag/len(dd):.1%} (帰無50%, p={p:.4g})"
    return out


def main():
    cat_tab = DE.loc["cat_rna_ctx_late"]

    t1 = build_table(WATCHLIST, cat_tab)
    t1.to_csv(RES / "candidates_human_validation.csv", index=False)
    log.info("watchlist %d 遺伝子 -> results/candidates_human_validation.csv", len(t1))

    stage = pd.read_csv(RES / "stage_classes_ctx.csv", index_col=0)
    late = stage[stage["class"] == "late_onset"].index.tolist()
    t2 = build_table(late, cat_tab)
    t2.to_csv(RES / "late_onset_human_validation.csv", index=False)
    log.info("後期起動 %d 遺伝子 -> results/late_onset_human_validation.csv", len(t2))

    r1 = sign_tests(t1, f"watchlist {len(t1)}遺伝子")
    r2 = sign_tests(t2, f"後期起動 {len(t2)}遺伝子")
    for r in (r1, r2):
        log.info("%s", r)

    append_summary("10_human_validation / 方向一致の符号検定", {
        **{f"[watchlist] {k}": v for k, v in r1.items()},
        **{f"[後期起動] {k}": v for k, v in r2.items()},
        "注意1": ("KPMP の CKD は糖尿病性腎症・高血圧性腎硬化症が主体で、"
                  "ネコ自然発症CKDとは病因が異なる。一致しなかった遺伝子は"
                  "「否定された」ではなく「病因の異なるヒトCKDでは確認できなかった」。"),
        "注意2": ("遺伝子は共発現により独立でないため、二項検定の p 値は反保守的。"
                  "一致率（効果量）と併せて読むこと。"),
        "注意3": "個別遺伝子の有意性ではなく、方向一致そのものを所見として扱っている。",
    }, cfg)


if __name__ == "__main__":
    main()
