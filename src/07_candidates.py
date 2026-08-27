"""尿バイオマーカー候補の選抜。

選抜ロジック（すべて事前定義。後付けで基準を動かさないこと）:
  1. ネコで「後期起動」かつ上昇（進行の指標になりうる）
  2. UniProt でシグナルペプチド保有 or 分泌／細胞外に注釈（＝遊離蛋白として尿に出る筋）
  3. マウス腎プロテオームで検出され、方向が一致
  4. マウス腎で局所転写されている（血漿由来の受動蓄積と区別する）
     ※ 併走論文の plasma-protein-handling signature と混同しないための必須チェック
  5. 皮質・髄質の両方で動く（区画依存性が低いほど尿での再現性が高いと期待）
  6. LOO で符号が安定

膜貫通型（UT-A、AQP6、KCNJ1 等）は遊離蛋白では測れないため、
uEV（尿中細胞外小胞）トラックとして別建てで出力する。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary, alias_lookup  # noqa: E402

log = get_logger("07_candidates")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
REF = ROOT / cfg["paths"]["ref"]
RES = ROOT / cfg["paths"]["results"]
TH = cfg["thresholds"]

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))
CAT2ALIAS = dict(zip(OMAP.cat_symbol, OMAP.mouse_aliases))


def tab(ds):
    return DE.loc[ds]


def load_secretome():
    frames = []
    for sp in ["mouse", "human", "cat"]:
        p = REF / f"secretome_{sp}.tsv"
        if p.exists():
            d = pd.read_csv(p, sep="\t")
            d["species"] = sp
            frames.append(d)
    if not frames:
        log.error("secretome_*.tsv がありません。00_fetch_refs.py を実行してください")
        return pd.DataFrame(columns=["symbol", "has_signal_peptide", "is_secreted",
                                     "n_transmembrane", "species"])
    return pd.concat(frames, ignore_index=True)


def main():
    stage = pd.read_csv(RES / "stage_classes_ctx.csv", index_col=0)
    stage_med = pd.read_csv(RES / "stage_classes_med.csv", index_col=0)
    late_up = stage[(stage["class"] == "late_onset") & (stage.lfc_late > 0)].copy()
    log.info("ネコ皮質 後期起動・上昇: %d 遺伝子", len(late_up))

    sec = load_secretome()
    sec_sym = sec.groupby(sec.symbol.str.upper()).agg(
        signal=("has_signal_peptide", "max"),
        secreted=("is_secreted", "max"),
        tm=("n_transmembrane", "max"),
    )

    mp21, mp14 = tab("m_prot_d21"), tab("m_prot_d14")
    mr2w = tab("m_rna_2w")
    mouse_rna_mat = pd.read_parquet(INT / "mouse_rna.parquet")
    mgrp = pd.read_parquet(INT / "mouse_rna_grp.parquet").squeeze()
    ctrl_fpkm = mouse_rna_mat[mgrp[mgrp == "Ctrl"].index].mean(axis=1)
    d2w_fpkm = mouse_rna_mat[mgrp[mgrp == "DT100ng-2W"].index].mean(axis=1)

    rows = []
    for cat_gene, r in late_up.iterrows():
        mouse_gene = CAT2MOUSE.get(cat_gene)
        alias = CAT2ALIAS.get(cat_gene)
        # 各突合先はそれぞれ別世代のシンボルを使っているため、必ず別名を総当たりする。
        sec_key = alias_lookup(sec_sym.index, alias, mouse_gene or cat_gene)
        ann = sec_sym.loc[sec_key] if sec_key is not None else None
        med_r = stage_med.loc[cat_gene] if cat_gene in stage_med.index else None
        rec = {
            "cat_symbol": cat_gene,
            "mouse_symbol": mouse_gene,
            "cat_ctx_early": round(r.lfc_early, 3),
            "cat_ctx_late": round(r.lfc_late, 3),
            "cat_ctx_q_late": r.q_late,
            "cat_med_late": round(med_r.lfc_late, 3) if med_r is not None else np.nan,
            "cat_med_class": med_r["class"] if med_r is not None else None,
            "has_signal_peptide": bool(ann["signal"]) if ann is not None else None,
            "is_secreted": bool(ann["secreted"]) if ann is not None else None,
            "n_transmembrane": float(ann["tm"]) if ann is not None else np.nan,
        }
        p21 = alias_lookup(mp21.index, alias, mouse_gene)
        if p21 is not None:
            p14 = alias_lookup(mp14.index, alias, mouse_gene)
            rec["mouse_prot_d14"] = round(float(mp14.loc[p14, "lfc"]), 3) if p14 is not None else np.nan
            rec["mouse_prot_d21"] = round(float(mp21.loc[p21, "lfc"]), 3)
            rec["mouse_prot_q"] = float(mp21.loc[p21, "q"])
            rec["mouse_prot_sign_stable"] = bool(mp21.loc[p21, "sign_stable"])
        r2w = alias_lookup(mr2w.index, alias, mouse_gene)
        if r2w is not None:
            rec["mouse_rna_2w"] = round(float(mr2w.loc[r2w, "lfc"]), 3)
        fk = alias_lookup(ctrl_fpkm.index, alias, mouse_gene)
        if fk is not None:
            rec["mouse_fpkm_ctrl"] = round(float(ctrl_fpkm[fk]), 2)
            rec["mouse_fpkm_2w"] = round(float(d2w_fpkm[fk]), 2)
        rows.append(rec)

    cand = pd.DataFrame(rows)

    # ---------- スコアリング ----------
    def score(row):
        s = 0
        if row.get("is_secreted") or row.get("has_signal_peptide"):
            s += 2
        if pd.notna(row.get("n_transmembrane")) and row["n_transmembrane"] == 0:
            s += 1
        if pd.notna(row.get("mouse_prot_d21")) and row["mouse_prot_d21"] > 0:
            s += 2
            if row.get("mouse_prot_q", 1) < 0.05:
                s += 1
            if row.get("mouse_prot_sign_stable"):
                s += 1
        if pd.notna(row.get("mouse_rna_2w")) and row["mouse_rna_2w"] > 0.5:
            s += 2   # 局所転写あり＝血漿由来の受動蓄積ではない
        if pd.notna(row.get("cat_med_late")) and row["cat_med_late"] > 0:
            s += 1   # 両区画で上昇
        return s

    cand["score"] = cand.apply(score, axis=1)
    cand["track"] = np.where(
        (cand.n_transmembrane.fillna(0) > 0) & ~cand.is_secreted.fillna(False),
        "uEV", "free_urine")
    cand = cand.sort_values("score", ascending=False)
    cand.to_csv(RES / "urinary_candidates.csv", index=False)

    top = cand[(cand.track == "free_urine") & (cand.score >= 6)]
    log.info("遊離尿蛋白トラック 上位候補:\n%s",
             top.head(20)[["cat_symbol", "mouse_symbol", "cat_ctx_late",
                           "mouse_prot_d21", "mouse_fpkm_ctrl", "mouse_fpkm_2w", "score"]]
             .to_string(index=False))

    # ---------- 個別確認が必要な既知候補 ----------
    # 末尾3つは Lourenço et al. 2021（虚血誘発ネコCKD）が提示した候補
    watch = ["FAM3D", "PTN", "ANXA3", "NPNT", "ANGPTL2", "THBS1", "SULF1", "ITGB6",
             "SLC14A2", "AQP2", "AQP6", "KCNJ1", "SCNN1G", "FKBP5", "SGK1", "NR3C2",
             "FBLN1", "SPP1", "MGP"]
    w = cand[cand.cat_symbol.isin(watch)]
    w.to_csv(RES / "urinary_candidates_watchlist.csv", index=False)

    append_summary("07_candidates", {
        "後期起動・上昇の総数": len(late_up),
        "分泌/シグナルペプチド保有": int(cand.is_secreted.fillna(False).sum()),
        "遊離尿トラック上位（score>=6）": top.cat_symbol.tolist()[:20],
        "uEVトラック": cand.query("track=='uEV'").cat_symbol.tolist()[:20],
        "注意": ("FAM3D のマウスオルソログは Oit1。オルソログ表を使わないと必ず落ちる。"
                 "また FAM3D は消化管由来で栄養状態により血中濃度が変動するため、"
                 "尿マーカー化する場合は採尿タイミングの影響を先に評価すること。"),
    }, cfg)


if __name__ == "__main__":
    main()
