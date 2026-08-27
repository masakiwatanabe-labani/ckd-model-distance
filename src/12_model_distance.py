"""マウスモデル間距離の「分布」に対して、種間距離がどこに位置するかを測る。

中心的な問い:
    種間距離（ネコ自然発症 × マウス各モデル）は、
    マウスモデル間距離の分布の「中に埋もれる」か「外にある」か。

d_entry という単一の点推定の比較ではなく、分布に対する位置（パーセンタイル）で答える。
点推定の比較は時点の取り方で 1.01〜1.32 に振れることが分かっているため。

モデル構成（すべて RNA、各モデルは自身の対照に対する log2FC）:
    PodTRECK : 5D / 2W / 3W    糸球体入口（足細胞特異的傷害）  対照 = 非誘導 Ctrl
    IRI      : 6mo / 12mo      尿細管入口（虚血）              対照 = 同一プラットフォームの sham/normal
    UUO      : 2D / 8D         尿細管入口（閉塞）              対照 = sham 4腎

対照の取り方（距離に効くので明記する）:
  - UUO (GSE79443): 右尿管結紮。**閉塞腎（右）のみが寄託されており、UUO個体の対側左腎は
    データセットに含まれない**ため、対側腎対照は選べない。sham 手術腎を対照とする。
    sham は個体A・Bの左右4腎。左右差は無視できる（rho 0.992、|差|中央値 0.19 log2、
    個体内左右 0.987 ≒ 個体間同側 0.989）ため4腎を対照に用いる。
    ただし 2個体×2腎の擬似反復である点は制約として残る。
    側一致 sham 右のみ（n=2）を感度解析として併記する。
  - IRI (GSE98622): 全腎。同一プラットフォーム・同週齢の sham を対照（09と同じ定義）。
  - PodTRECK: 非誘導対照。

バッチ交絡:
  3データセットは施設・時期が異なる。log2FC は各データセット内のコントラストとして
  計算しているため、バッチの主効果は相殺されるが、バッチ×病態の交互作用は残る。
  距離の絶対値ではなく分布内での相対位置で解釈すること。
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, moderated_ttest, append_summary  # noqa: E402

log = get_logger("12_model_distance")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
EXT = ROOT / "data" / "external"

import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("m09", Path(__file__).parent / "09_dataset_matrix.py")
m09 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m09)


def uuo_lfc(ctrl_cols=None, tag="uuo", day="8D") -> pd.Series | None:
    p = EXT / "GSE79443" / "GSE79443_SO_2D_8D_norm_counts.txt"
    if not p.exists():
        log.error("GSE79443 がありません: %s", p)
        return None
    d = pd.read_csv(p, sep="\t")
    mat = d.set_index("gene_name")[["SO_AL1", "SO_AR1", "SO_BL1", "SO_BR1",
                                     "2D_AR1", "2D_BR1", "2D_CR1",
                                     "8D_AR1", "8D_BR1", "8D_CR1"]]
    mat = mat.groupby(level=0).sum()
    ctrl_cols = ctrl_cols or ["SO_AL1", "SO_AR1", "SO_BL1", "SO_BR1"]
    case_cols = [c for c in mat.columns if c.startswith(day)]
    keep = mat[ctrl_cols + case_cols].max(axis=1) >= 10
    m = np.log2(mat[keep] + 1)
    A, B = m[case_cols], m[ctrl_cols]
    res = moderated_ttest(A, B)
    log.info("UUO %s: case n=%d vs sham n=%d, 遺伝子 %d, q<0.05 %d",
             day, len(case_cols), len(ctrl_cols), len(res), int((res.q < 0.05).sum()))
    return res["lfc"]


def build_states() -> tuple[dict, dict]:
    """モデル状態 -> ヒト空間 log2FC 系列、および 状態 -> モデル名。"""
    S, model = {}, {}
    for key, ds in [("PodTRECK_5D", "m_rna_5d"), ("PodTRECK_2W", "m_rna_2w"),
                    ("PodTRECK_3W", "m_rna_3w")]:
        S[key] = m09.to_human_space(m09.internal_lfc(ds), "mouse")
        model[key] = "PodTRECK"
    v = m09.mouse_iri_lfc(tag="late")
    if v is not None:
        S["IRI_12mo"] = m09.to_human_space(v, "mouse"); model["IRI_12mo"] = "IRI"
    v = m09.mouse_iri_lfc(m09.IRI_LATE_SENS, m09.IRI_CTRL_SENS, tag="late_sens")
    if v is not None:
        S["IRI_6mo"] = m09.to_human_space(v, "mouse"); model["IRI_6mo"] = "IRI"
    for day, key in [("2D", "UUO_2D"), ("8D", "UUO_8D")]:
        v = uuo_lfc(day=day)
        if v is not None:
            S[key] = m09.to_human_space(v, "mouse"); model[key] = "UUO"
    return S, model


def dist(x, y):
    idx = x.index.intersection(y.index)
    if len(idx) < 50:
        return np.nan, len(idx)
    return 1 - stats.spearmanr(x.loc[idx], y.loc[idx])[0], len(idx)


def main():
    S, model = build_states()
    log.info("マウスモデル状態: %s", {k: model[k] for k in S})

    cats = {"cat_natural_ctx": m09.to_human_space(m09.internal_lfc("cat_rna_ctx_late"), "cat"),
            "cat_natural_med": m09.to_human_space(m09.internal_lfc("cat_rna_med_late"), "cat")}

    rows = []
    for a, b in combinations(list(S), 2):
        d, n = dist(S[a], S[b])
        rows.append({"a": a, "b": b, "kind": "within_model" if model[a] == model[b] else "between_model",
                     "distance": d, "n_genes": n})
    for c, cs in cats.items():
        for k, v in S.items():
            d, n = dist(cs, v)
            rows.append({"a": c, "b": k, "kind": "cross_species", "distance": d, "n_genes": n})
    d_cc, n_cc = dist(cats["cat_natural_ctx"], cats["cat_natural_med"])
    rows.append({"a": "cat_natural_ctx", "b": "cat_natural_med",
                 "kind": "within_cat", "distance": d_cc, "n_genes": n_cc})
    t = pd.DataFrame(rows)
    t.to_csv(RES / "model_distance_pairs.csv", index=False)

    between = t[t.kind == "between_model"].distance.dropna()
    within = t[t.kind == "within_model"].distance.dropna()
    log.info("モデル間距離 (n=%d): min %.3f / 中央 %.3f / max %.3f",
             len(between), between.min(), between.median(), between.max())
    log.info("モデル内距離 (n=%d): min %.3f / 中央 %.3f / max %.3f",
             len(within), within.min(), within.median(), within.max())

    out = {}
    for c in cats:
        sub = t[(t.kind == "cross_species") & (t.a == c)]
        recs = []
        for _, r in sub.iterrows():
            pct = float((between < r.distance).mean() * 100)
            recs.append({"mouse_state": r.b, "d_species": round(r.distance, 4),
                         "percentile_in_between_model": round(pct, 1),
                         "n_genes": int(r.n_genes),
                         "判定": ("分布の外（モデル間より遠い）" if r.distance > between.max()
                                  else "分布の中（モデル間距離に埋もれる）" if r.distance >= between.min()
                                  else "分布の外（モデル間より近い）")})
        out[c] = pd.DataFrame(recs).sort_values("d_species")
        log.info("\n=== %s ===\n%s", c, out[c].to_string(index=False))

    pd.concat([v.assign(cat=k) for k, v in out.items()]).to_csv(
        RES / "species_vs_model_distance.csv", index=False)

    # 感度解析: UUO の対照を側一致 sham 右のみ（n=2）にする
    sens = {}
    v = uuo_lfc(ctrl_cols=["SO_AR1", "SO_BR1"], day="8D")
    if v is not None:
        u = m09.to_human_space(v, "mouse")
        d_all, _ = dist(cats["cat_natural_ctx"], S["UUO_8D"])
        d_r, _ = dist(cats["cat_natural_ctx"], u)
        sens = {"UUO_8D 対照=sham4腎": round(d_all, 4), "UUO_8D 対照=sham右のみ(n=2)": round(d_r, 4)}
        log.info("UUO 対照の取り方による感度: %s", sens)

    append_summary("12_model_distance / 種間距離はモデル間距離分布の中か外か", {
        "マウスモデル状態": {k: model[k] for k in S},
        "モデル間距離 (between_model)": (
            f"n={len(between)}, min={between.min():.3f}, 中央={between.median():.3f}, "
            f"max={between.max():.3f}, IQR={between.quantile(.25):.3f}-{between.quantile(.75):.3f}"),
        "モデル内距離 (within_model, 参考)": (
            f"n={len(within)}, min={within.min():.3f}, 中央={within.median():.3f}, max={within.max():.3f}"),
        "ネコ皮質 × 各マウスモデル": out["cat_natural_ctx"].to_dict("records"),
        "ネコ髄質 × 各マウスモデル": out["cat_natural_med"].to_dict("records"),
        "ネコ皮質×髄質（種内参考）": round(d_cc, 4),
        "UUO対照の感度": sens,
        "制約1（UUOの対照）": ("GSE79443 は右尿管結紮で、**閉塞腎（右）のみ寄託**されており "
                               "UUO個体の対側左腎は存在しない。したがって対側腎対照は選べず、"
                               "sham手術腎を対照とした。sham は2個体×左右4腎で、左右差は"
                               "無視できる（rho 0.992、|差|中央値 0.19 log2）ため4腎を用いたが、"
                               "2個体の擬似反復である点は残る。"),
        "制約2（バッチ交絡）": ("PodTRECK / GSE98622 / GSE79443 は施設・時期が異なる。"
                                "log2FC を各データセット内のコントラストとして計算しているため"
                                "バッチ主効果は相殺されるが、バッチ×病態の交互作用は残る。"
                                "距離の絶対値ではなく分布内の相対位置で解釈すること。"),
        "制約3（Gubra論文）": ("Marstrand-Jørgensen et al., Nephron 2024 (10.1159/000535918) は"
                               "同一施設で UUO/uIRI/ADI 3モデルが揃い理想的だったが、"
                               "Data Availability Statement により公開寄託されていない"
                               "（corresponding author への請求のみ）。GEO/SRA/ArrayExpress/"
                               "BioProject を検索したが該当なし。代替として GSE79443 を使用。"),
    }, cfg)


if __name__ == "__main__":
    main()
