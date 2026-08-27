"""Results 3.4: 16状態（マウス14＋ネコ2）で距離分布を組み直す。

3.3 と状態集合を統一する。3.3 は14状態版（マウスのみ）、
13番の分布は7状態版由来で食い違っていたため、ここで16状態に揃える。

ペアの分類:
  between_model : 異なるモデル同士のマウスペア（PodTRECK / IRI / UUO）
  within_model  : 同一モデルの異なる時点のペア（IRI9時点あるので大半がここ）
  cross_species : ネコ2状態 × マウス14状態 = 28ペア
  （ネコ皮質×髄質の1ペアは内部基準として別扱い）

注意:
  IRI が9時点あるため within_model の大半（36/40）が IRI-IRI 間になる。
  この偏りは分布の解釈に影響するので、IRI を含む版と除く版の両方を出す。
  本文の主張は between_model と cross_species の比較なので within_model は参考値。
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_figure import apply_style, savefig as _savefig  # noqa: E402
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("19_distance_distributions")
apply_style()
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]
FIG = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m16", Path(__file__).parent / "16_entry_model_checks.py")
m16 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m16)
m12, m09 = m16.m12, m16.m09
S, DAYS, MODEL, ENTRY, STATES = m16.S, m16.DAYS, m16.MODEL, m16.ENTRY, m16.STATES
CAT = m16.CAT

LABEL = {"PodTRECK_5D": "PodTRECK 5d", "PodTRECK_2W": "PodTRECK 14d", "PodTRECK_3W": "PodTRECK 21d",
         "UUO_2D": "UUO 2d", "UUO_8D": "UUO 8d", "IRI2h": "IRI 2h", "IRI4h": "IRI 4h",
         "IRI24h": "IRI 24h", "IRI48h": "IRI 48h", "IRI72h": "IRI 72h", "IRI7d": "IRI 7d",
         "IRI14d": "IRI 14d", "IRI28d": "IRI 28d", "IRI12m": "IRI 12mo",
         "cat_natural_ctx": "cat cortex", "cat_natural_med": "cat medulla"}


def desc(x, name):
    x = np.asarray(x)
    return {"group": name, "n": len(x), "min": round(float(np.min(x)), 4),
            "Q1": round(float(np.percentile(x, 25)), 4),
            "median": round(float(np.median(x)), 4),
            "Q3": round(float(np.percentile(x, 75)), 4),
            "max": round(float(np.max(x)), 4)}


def main():
    pool = {**S, **CAT}
    model = {**MODEL, "cat_natural_ctx": "cat", "cat_natural_med": "cat"}
    allst = STATES + list(CAT)

    rows = []
    for a, b in combinations(allst, 2):
        d = m16.d(a, b, pool)
        ca, cb = model[a] == "cat", model[b] == "cat"
        if ca and cb:
            kind = "within_cat"
        elif ca or cb:
            kind = "cross_species"
        elif model[a] == model[b]:
            kind = "within_model"
        else:
            kind = "between_model"
        rows.append({"a": a, "b": b, "label_a": LABEL[a], "label_b": LABEL[b],
                     "model_a": model[a], "model_b": model[b], "kind": kind, "distance": d})
    t = pd.DataFrame(rows)
    t.to_csv(RES / "distance_distributions_16state.csv", index=False)
    log.info("総ペア %d: %s", len(t), t.kind.value_counts().to_dict())

    bet = t[t.kind == "between_model"].distance.values
    wit = t[t.kind == "within_model"]
    wit_all = wit.distance.values
    wit_noiri = wit[~((wit.model_a == "IRI") & (wit.model_b == "IRI"))].distance.values
    crs = t[t.kind == "cross_species"].distance.values
    catref = float(t[t.kind == "within_cat"].distance.iloc[0])

    summary = [desc(bet, "between_model"), desc(wit_all, "within_model (IRI込み)"),
               desc(wit_noiri, "within_model (IRI-IRI除く)"), desc(crs, "cross_species")]
    st = pd.DataFrame(summary)
    st.to_csv(RES / "distance_distribution_summary.csv", index=False)
    log.info("\n%s", st.to_string(index=False))
    log.info("ネコ皮質×髄質（内部基準）= %.4f", catref)

    # ---- 主張の確認 ----
    inside = int(((crs >= bet.min()) & (crs <= bet.max())).sum())
    pct = [float((bet < c).mean() * 100) for c in crs]
    cs = t[t.kind == "cross_species"].copy()
    cs["percentile_in_between_model"] = np.round(pct, 1)
    cs = cs.sort_values("distance")
    cs.to_csv(RES / "cross_species_percentiles_16state.csv", index=False)
    u, p_mw = stats.mannwhitneyu(crs, bet, alternative="greater")
    far_b = t[t.kind == "between_model"].nlargest(1, "distance").iloc[0]
    far_c = cs.nlargest(1, "distance").iloc[0]

    claim = {
        "種間がモデル間レンジ内に収まる数": f"{inside} / {len(crs)}",
        "モデル間レンジ": f"{bet.min():.4f} - {bet.max():.4f}",
        "種間レンジ": f"{crs.min():.4f} - {crs.max():.4f}",
        "Mann-Whitney(種間>モデル間, 片側) p": round(float(p_mw), 4),
        "中央値": f"種間 {np.median(crs):.4f} vs モデル間 {np.median(bet):.4f}",
        "モデル間の最遠ペア": f"{far_b.label_a} x {far_b.label_b} = {far_b.distance:.4f}",
        "種間の最遠ペア": f"{far_c.label_a} x {far_c.label_b} = {far_c.distance:.4f}",
        "種間パーセンタイル範囲": f"{min(pct):.1f} - {max(pct):.1f}",
    }
    for k, v in claim.items():
        log.info("  %s: %s", k, v)
    log.info("\n種間28ペアのパーセンタイル:\n%s",
             cs[["label_a", "label_b", "distance", "percentile_in_between_model"]]
             .to_string(index=False))

    # ---- 7状態版との比較 ----
    old = pd.read_csv(RES / "model_distance_pairs.csv")
    ob = old[old.kind == "between_model"].distance.dropna().values
    oc = old[old.kind == "cross_species"].distance.dropna().values
    o_inside = int(((oc >= ob.min()) & (oc <= ob.max())).sum())
    _, op = stats.mannwhitneyu(oc, ob, alternative="greater")
    comp = pd.DataFrame([
        {"版": "7状態(13番)", "between_n": len(ob), "between_median": round(float(np.median(ob)), 4),
         "between_max": round(float(ob.max()), 4), "cross_n": len(oc),
         "cross_median": round(float(np.median(oc)), 4), "cross_max": round(float(oc.max()), 4),
         "レンジ内": f"{o_inside}/{len(oc)}", "MW_p": round(float(op), 4)},
        {"版": "16状態(本節)", "between_n": len(bet), "between_median": round(float(np.median(bet)), 4),
         "between_max": round(float(bet.max()), 4), "cross_n": len(crs),
         "cross_median": round(float(np.median(crs)), 4), "cross_max": round(float(crs.max()), 4),
         "レンジ内": f"{inside}/{len(crs)}", "MW_p": round(float(p_mw), 4)}])
    comp.to_csv(RES / "distance_distribution_7v16.csv", index=False)
    log.info("\n7状態 vs 16状態:\n%s", comp.to_string(index=False))

    # ================= Fig 4 =================
    groups = [("between-model\n(mouse)", bet, "#4C72B0"),
              ("within-model\n(mouse)", wit_all, "#55A868"),
              ("cross-species\n(cat x mouse)", crs, "#C44E52")]
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    parts = ax.violinplot([g[1] for g in groups], showextrema=False, widths=.82)
    for pc, g in zip(parts["bodies"], groups):
        pc.set_facecolor(g[2]); pc.set_alpha(.30)
    bp = ax.boxplot([g[1] for g in groups], widths=.20, patch_artist=True, showfliers=False)
    for b, g in zip(bp["boxes"], groups):
        b.set_facecolor(g[2]); b.set_alpha(.85)
    for m in bp["medians"]:
        m.set_color("black"); m.set_linewidth(1.3)
    rng = np.random.default_rng(SEED)
    for i, g in enumerate(groups, start=1):
        ax.scatter(rng.normal(i, .052, len(g[1])), g[1], s=13, color=g[2],
                   alpha=.62, zorder=3, edgecolor="none")
    ax.axhline(catref, ls="--", lw=1.4, color="0.25", zorder=4)
    ax.text(3.46, catref, f"cat cortex vs medulla\n(internal reference) = {catref:.3f}",
            fontsize=7.8, va="center", ha="right", color="0.25")
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"{g[0]}\nn={len(g[1])}" for g in groups], fontsize=9)
    ax.set_ylabel("distance (1 - Spearman rho)", fontsize=10)
    ax.set_title("Cross-species distances fall inside the between-model distribution",
                 fontsize=11.5, loc="left")
    ax.grid(axis="y", alpha=.25)
    fig.tight_layout()
    _savefig(fig, FIG / "Fig4_distance_distributions")
    log.info("saved %s", FIG / "Fig4_distance_distributions.png")

    append_summary("19_distance_distributions / Results3.4", {
        "状態集合": "16状態 = マウス14（PodTRECK3 + UUO2 + IRI9）+ ネコ2。3.3と統一。",
        "ペア内訳": t.kind.value_counts().to_dict(),
        "分布": st.to_dict("records"),
        "主張の確認": claim,
        "7状態版との比較": comp.to_dict("records"),
        "内部基準（ネコ皮質×髄質）": {
            "値": round(catref, 4),
            "空間": "ヒトシンボル空間、to_human_space（16状態版と同一定義）",
            "遺伝子数": 14605,
            "デザイン": ("晩期コントラスト（CKD3/4 vs Control）の log2FC 同士の相関であり、"
                         "17頭の対応ありデザインではない。"),
            "06_compartmentとの関係": ("06 の『晩期 皮質×髄質 rho 0.763』(距離0.237) は"
                                       "ネコシンボル空間・対応ありデザインで計算した別の値。"
                                       "06 の主要数値『進行軸 rho 0.204』は進行軸同士の相関で"
                                       "さらに別物。3.3/3.4 では 0.2401 を使うこと。")},
        "within_modelの偏り": (f"IRI が9時点あるため within_model {len(wit_all)}ペア中 "
                               f"{len(wit_all)-len(wit_noiri)} が IRI-IRI 間。"
                               "IRI-IRI を除くと n=%d。参考値扱い。" % len(wit_noiri)),
        "図": "results/figures/Fig4_distance_distributions.png (300 dpi)",
    }, cfg)


if __name__ == "__main__":
    main()
