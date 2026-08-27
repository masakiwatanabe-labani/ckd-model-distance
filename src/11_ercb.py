"""ERCB 尿細管間質トランスクリプトーム（GSE104954）での候補遺伝子検証。

KPMP プロテオームで FAM3D/PTN が測定されていなかったため、転写レベルの
独立したヒトコホートで確認する。ERCB は腎生検の尿細管間質分画で、
病因別の診断名が付いているのが利点（KPMP は DKD/高血圧性が混在）。

重要な制約:
  1. マイクロアレイなので RNA-seq / プロテオームとは正規化が異なる。
     距離行列には入れず、個別遺伝子の検証にのみ使う。
  2. プラットフォームが2つあり、対照(LD)の配置が極端に偏っている
       GPL22945: LD 18, TN 3 / RPGN 21, DN 7, MCD 4, FSGS&MCD 4, FSGS 3
       GPL24120: LD  3, TN 2 / SLE 32, IgA 25, HT 20, MGN 18, DN 10, FSGS 10, MCD 9, TMD 6
     プラットフォームをまたぐ対照プールは batch と交絡するため禁止。
     必ず同一プラットフォーム内でコントラストを組む。
  3. プローブセット ID は "<EntrezID>_at"（Brainarray ENTREZG カスタム CDF）。
     1プローブセット = 1遺伝子で、複数プローブの集約は不要。
     両プラットフォームとも同一の 12,074 プローブセットを持つ。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, moderated_ttest, append_summary  # noqa: E402

log = get_logger("11_ercb")
cfg = load_config()
ROOT = Path(cfg["_root"])
REF = ROOT / cfg["paths"]["ref"]
RES = ROOT / cfg["paths"]["results"]
EXT = ROOT / "data" / "external" / "GSE104954"

TARGETS = ["FAM3D", "PTN", "ANXA3", "TXNIP", "ANGPTL2"]
# 参考として一緒に見る（後期起動で3種一致した分子＋既知マーカー）
EXTRA = ["CRYAB", "H6PD", "SPP1", "HAVCR1", "LCN2", "FBLN1", "THBS1", "SULF1", "NPNT"]

# タイトル末尾のトークン -> 病因ラベル
GROUP_MAP = {
    "LD": "control_LD",        # living donor（健常対照）
    "TN": "control_TN",        # tumor nephrectomy（腫瘍腎摘の非腫瘍部）
    "DN": "DN",                # 糖尿病性腎症（糸球体入口）
    "HT": "HT",                # 高血圧性腎硬化症（血管性・尿細管間質線維化）
    "RPGN": "RPGN",            # ANCA関連血管炎／半月体形成性
    "IgA": "IgAN",
    "SLE": "SLE",
    "MGN": "MGN",
    "FSGS": "FSGS",
    "FSGS&MCD": "FSGS_MCD",
    "MCD": "MCD",
    "TMD": "TMD",
}


def load_platform(stem: str):
    # GEO からの自動取得では "<stem>_series_matrix.txt" という名前で降ってくるが、
    # 手動で置く場合は "<stem>.txt" にしていることもある。両方を受ける。
    cands = [EXT / f"{stem}.txt", EXT / f"{stem}_series_matrix.txt"]
    cands += sorted(EXT.glob(f"{stem}*series_matrix.txt"))
    path = next((c for c in cands if c.exists()), None)
    if path is None:
        raise FileNotFoundError(
            f"{stem} の series matrix が {EXT} にありません。"
            "src/00b_fetch_external.py を実行するか data/README.md を参照してください。")
    lines = path.read_text(encoding="utf-8", errors="ignore").split("\n")

    def meta(pref):
        for L in lines:
            if L.startswith(pref):
                return [x.strip('"') for x in L.split("\t")[1:]]
        return []

    gsm, title = meta("!Sample_geo_accession"), meta("!Sample_title")
    mat = pd.read_csv(path, sep="\t", comment="!", index_col=0, low_memory=False)
    mat = mat[mat.index.notna()]
    mat.columns = [str(c).strip('"') for c in mat.columns]

    import re
    grp = {}
    for g, t in zip(gsm, title):
        tok = re.sub(r"\d+$", "", t.split("-")[-1])
        grp[g] = GROUP_MAP.get(tok, f"unknown_{tok}")
    grp = pd.Series(grp, name="group")
    grp = grp[grp.index.isin(mat.columns)]
    return mat[grp.index], grp


def entrez_to_symbol() -> dict:
    p = REF / "entrez2symbol_human.tsv"
    if not p.exists():
        log.error("%s がありません。00_fetch_refs.py を実行してください", p)
        return {}
    d = pd.read_csv(p, sep="\t")
    return dict(zip(d.entrez_id.astype(str), d.symbol.astype(str)))


def contrast(mat, grp, case_labels, ctrl_labels, label):
    A = mat[grp[grp.isin(case_labels)].index]
    B = mat[grp[grp.isin(ctrl_labels)].index]
    if A.shape[1] < 3 or B.shape[1] < 3:
        log.warning("%s: n が足りません (case=%d ctrl=%d)。スキップ", label, A.shape[1], B.shape[1])
        return None
    res = moderated_ttest(A, B)
    res["contrast"] = label
    res["n_case"], res["n_ctrl"] = A.shape[1], B.shape[1]
    log.info("%s: case n=%d vs ctrl n=%d, q<0.05: %d",
             label, A.shape[1], B.shape[1], int((res.q < 0.05).sum()))
    return res


def main():
    e2s = entrez_to_symbol()
    out = []
    platforms = {}
    for stem in ["GSE104954-GPL22945", "GSE104954-GPL24120"]:
        mat, grp = load_platform(stem)
        sym = [e2s.get(str(i).replace("_at", ""), f"ENTREZ:{str(i).replace('_at','')}")
               for i in mat.index]
        mat = mat.set_axis(sym, axis=0)
        mat = mat.groupby(level=0).mean()   # Entrez CDF なので通常は1対1、保険
        platforms[stem] = (mat, grp)
        log.info("%s: %s, 群構成 %s", stem, mat.shape, grp.value_counts().to_dict())

    # ---- 同一プラットフォーム内のコントラストのみ ----
    plans = [
        # (platform, case, ctrl, label)
        ("GSE104954-GPL22945", ["DN"], ["control_LD"], "GPL22945: DN vs LD"),
        ("GSE104954-GPL22945", ["RPGN"], ["control_LD"], "GPL22945: RPGN(ANCA) vs LD"),
        ("GSE104954-GPL22945", ["FSGS", "FSGS_MCD", "MCD"], ["control_LD"],
         "GPL22945: 足細胞疾患(FSGS/MCD) vs LD"),
        ("GSE104954-GPL22945", ["DN", "RPGN", "FSGS", "FSGS_MCD", "MCD"], ["control_LD"],
         "GPL22945: 全CKD vs LD"),
        ("GSE104954-GPL24120", ["HT"], ["control_LD", "control_TN"], "GPL24120: HT vs LD+TN"),
        ("GSE104954-GPL24120", ["DN"], ["control_LD", "control_TN"], "GPL24120: DN vs LD+TN"),
        ("GSE104954-GPL24120", ["IgAN"], ["control_LD", "control_TN"], "GPL24120: IgAN vs LD+TN"),
        ("GSE104954-GPL24120", ["SLE"], ["control_LD", "control_TN"], "GPL24120: SLE vs LD+TN"),
        ("GSE104954-GPL24120",
         ["HT", "DN", "IgAN", "SLE", "MGN", "FSGS", "MCD", "TMD"],
         ["control_LD", "control_TN"], "GPL24120: 全CKD vs LD+TN"),
    ]
    for stem, case, ctrl, label in plans:
        mat, grp = platforms[stem]
        r = contrast(mat, grp, case, ctrl, label)
        if r is not None:
            out.append(r)

    allres = pd.concat(out)
    allres.to_csv(RES / "ercb_all_contrasts.csv.gz", compression="gzip")

    genes = TARGETS + EXTRA
    rows = []
    for r in out:
        lab = r["contrast"].iloc[0]
        for g in genes:
            if g in r.index:
                rows.append({"contrast": lab, "gene": g,
                             "n_case": int(r["n_case"].iloc[0]), "n_ctrl": int(r["n_ctrl"].iloc[0]),
                             "lfc": round(float(r.loc[g, "lfc"]), 3),
                             "p": float(r.loc[g, "p"]), "q": float(r.loc[g, "q"]),
                             "sig_q05": bool(r.loc[g, "q"] < 0.05)})
            else:
                rows.append({"contrast": lab, "gene": g, "n_case": int(r["n_case"].iloc[0]),
                             "n_ctrl": int(r["n_ctrl"].iloc[0]), "lfc": np.nan,
                             "p": np.nan, "q": np.nan, "sig_q05": None, "note": "未検出"})
    t = pd.DataFrame(rows)
    t.to_csv(RES / "ercb_candidate_validation.csv", index=False)
    log.info("saved results/ercb_candidate_validation.csv")

    mat0 = platforms["GSE104954-GPL22945"][0]
    detected = {g: ("検出" if g in mat0.index else "未検出（このアレイに非搭載）") for g in TARGETS}
    log.info("目的遺伝子の検出状況: %s", detected)

    append_summary("11_ercb / ERCB尿細管間質での候補検証", {
        "データ": "GSE104954 (ERCB tubulointerstitium, 腎生検)",
        "プラットフォーム構成": {s: platforms[s][1].value_counts().to_dict() for s in platforms},
        "目的遺伝子の検出": detected,
        "プローブ集約": ("Brainarray ENTREZG カスタムCDF。プローブセットID='<EntrezID>_at' で "
                         "1プローブセット=1遺伝子。複数プローブの集約は発生しない"
                         "（保険として同一シンボルは平均を取る実装）。両プラットフォームとも"
                         "同一の12,074プローブセット。"),
        "制約1（正規化）": ("マイクロアレイなので KPMP プロテオーム・RNA-seq とは正規化が異なる。"
                            "距離行列には入れず、個別遺伝子の検証にのみ使用。"),
        "制約2（プラットフォーム交絡）": ("対照(LD)は GPL22945 に18検体、GPL24120 に3検体と偏在。"
                                          "プラットフォームをまたぐ対照プールは batch と交絡するため、"
                                          "全コントラストを同一プラットフォーム内で組んでいる。"
                                          "GPL24120 側は対照 n=5（LD3+TN2）と少ない。"),
        "結果": "results/ercb_candidate_validation.csv 参照",
    }, cfg)


if __name__ == "__main__":
    main()
