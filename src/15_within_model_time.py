"""(1) IRI を完全に除いた版で時間効果の符号・大きさを見る。
(2) 各モデル内での「時点間距離 vs 経過時間差」を個別に描き、
    対数関係が3モデルすべてで成立するかを確認する。

注意:
  ネコには経過日数を割り当てられない（13_time_axis.py で確定: 自然発症で原点が
  定義できず、横断デザインで、IRISステージは腎機能の指標であって時間の指標ではない）。
  したがってネコは時間Mantelには入れられない。距離の参照として併記するに留める。

  IRI除外版は状態5個・ペア10個で検出力が無い。rho の符号と大きさのみを見る。
  UUO は2時点しか無いのでペアが1つだけで、モデル内の関係は評価できない。
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
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("15_within_model_time")
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]
# 本文の図に採用しなかった補足図はここに置く（Fig1-6 と混ぜない）。
SUP = RES / "supplementary"; SUP.mkdir(parents=True, exist_ok=True)
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m14", Path(__file__).parent / "14_time_axis_full.py")
m14 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m14)
m12, m09 = m14.m12, m14.m09


def main():
    # 全状態（IRI 9時点は GPL13112 内、6mo は別プラットフォームなので除く）
    S, days, model, entry = m14.build(include_6mo=False)

    # ---------- (1) IRI を完全に除いた版 ----------
    keep = [s for s in S if model[s] != "IRI"]
    Sk = {s: S[s] for s in keep}
    states, D, T, E, M = m14.matrices(Sk, days, model, entry)
    iu = np.triu_indices(len(states), 1)
    r_t = stats.spearmanr(D[iu], T[iu])[0]
    r_e = stats.spearmanr(D[iu], E[iu])[0]
    r_te, p_te = m14.mantel(D, T, seed=SEED, partial=E)
    log.info("[IRI除外] 状態 %d / ペア %d", len(states), len(iu[0]))
    log.info("  時間 rho=%+.3f / 入口 rho=%+.3f / 偏(時間|入口) rho=%+.3f (p=%.3f 参考)",
             r_t, r_e, r_te, p_te)
    noiri = pd.DataFrame({"a": [states[i] for i in iu[0]], "b": [states[j] for j in iu[1]],
                          "distance": D[iu], "dlog10days": T[iu],
                          "entry_mismatch": E[iu].astype(int)})
    noiri.to_csv(RES / "time_axis_noIRI_pairs.csv", index=False)
    log.info("\n%s", noiri.sort_values("distance").to_string(index=False))

    # ネコは時間軸に載せられないので、参照として距離のみ併記
    cats = {"cat_ctx_late": m09.to_human_space(m09.internal_lfc("cat_rna_ctx_late"), "cat"),
            "cat_med_late": m09.to_human_space(m09.internal_lfc("cat_rna_med_late"), "cat")}
    catrows = [{"cat": c, "mouse_state": s, "days": days[s],
                "distance": round(m12.dist(cs, Sk[s])[0], 4)}
               for c, cs in cats.items() for s in states]
    pd.DataFrame(catrows).to_csv(RES / "time_axis_noIRI_cat_ref.csv", index=False)

    # ---------- (2) モデル内の 距離 vs 経過時間差 ----------
    within = []
    for m in ["PodTRECK", "IRI", "UUO"]:
        ss = [s for s in S if model[s] == m]
        for a, b in combinations(ss, 2):
            d, n = m12.dist(S[a], S[b])
            within.append({"model": m, "a": a, "b": b, "distance": d,
                           "dlog10days": abs(np.log10(days[a]) - np.log10(days[b])),
                           "ddays": abs(days[a] - days[b]), "n_genes": n})
    w = pd.DataFrame(within)
    w.to_csv(RES / "within_model_time_pairs.csv", index=False)

    stats_rows = []
    for m in ["PodTRECK", "IRI", "UUO"]:
        d = w[w.model == m]
        if len(d) >= 3:
            rl = stats.spearmanr(d.distance, d.dlog10days)[0]
            rn = stats.spearmanr(d.distance, d.ddays)[0]
            pl = stats.spearmanr(d.distance, d.dlog10days)[1]
        else:
            rl = rn = pl = np.nan
        stats_rows.append({"model": m, "n_states": len([s for s in S if model[s] == m]),
                           "n_pairs": len(d),
                           "rho_vs_log10days": None if np.isnan(rl) else round(float(rl), 3),
                           "p_vs_log10days": None if np.isnan(pl) else round(float(pl), 4),
                           "rho_vs_lineardays": None if np.isnan(rn) else round(float(rn), 3),
                           "評価可能か": "可" if len(d) >= 3 else f"不可（ペア{len(d)}個）"})
    st = pd.DataFrame(stats_rows)
    log.info("\n%s", st.to_string(index=False))

    # ---------- プロット ----------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, m in zip(axes, ["PodTRECK", "IRI", "UUO"]):
        d = w[w.model == m]
        ax.scatter(d.dlog10days, d.distance, s=48, alpha=.8,
                   color={"PodTRECK": "#4C72B0", "IRI": "#DD8452", "UUO": "#55A868"}[m],
                   edgecolor="white", zorder=3)
        if len(d) >= 3:
            z = np.polyfit(d.dlog10days, d.distance, 1)
            xs = np.linspace(d.dlog10days.min(), d.dlog10days.max(), 50)
            ax.plot(xs, np.polyval(z, xs), "--", color="0.35", lw=1.2, zorder=2)
            r = stats.spearmanr(d.distance, d.dlog10days)
            sub = f"Spearman rho={r[0]:+.2f}, p={r[1]:.3f}  (n={len(d)} pairs)"
        else:
            sub = f"n={len(d)} pair — not evaluable"
        ax.set_title(f"{m}\n{sub}", fontsize=10)
        ax.set_xlabel("|$\\Delta$ log10(days post-injury)|")
        ax.set_ylabel("distance (1 - Spearman rho)")
        ax.grid(alpha=.25, zorder=0)
    fig.suptitle("Within-model: transcriptome distance vs elapsed-time difference", fontsize=12)
    fig.tight_layout()
    fig.savefig(SUP / "within_model_time.png", dpi=160)
    log.info("saved results/supplementary/within_model_time.png")

    # IRI の時間軸に沿った軌跡（対数関係の可視化）
    iri = [s for s in S if model[s] == "IRI"]
    iri_sorted = sorted(iri, key=lambda s: days[s])
    base = iri_sorted[0]
    traj = [{"state": s, "days": days[s], "distance_from_earliest": m12.dist(S[base], S[s])[0]}
            for s in iri_sorted]
    tr = pd.DataFrame(traj)
    tr.to_csv(RES / "iri_trajectory.csv", index=False)
    fig2, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(tr.days, tr.distance_from_earliest, "o-", color="#DD8452")
    ax.set_xscale("log")
    for _, r in tr.iterrows():
        ax.annotate(r.state.replace("IRI", ""), (r.days, r.distance_from_earliest),
                    textcoords="offset points", xytext=(4, 5), fontsize=8)
    ax.set_xlabel("days post-injury (log scale)")
    ax.set_ylabel(f"distance from {base}")
    ax.set_title("IRI time course: distance from earliest timepoint")
    ax.grid(alpha=.25)
    fig2.tight_layout()
    fig2.savefig(SUP / "iri_trajectory.png", dpi=160)
    log.info("saved results/supplementary/iri_trajectory.png")

    append_summary("15_within_model_time / IRI除外版とモデル内の時間関係", {
        "IRI除外版（PodTRECK3 + UUO2 = 5状態, 10ペア）": {
            "時間 rho": round(float(r_t), 3), "入口 rho": round(float(r_e), 3),
            "偏(時間|入口) rho": round(float(r_te), 3),
            "注記": "検出力が無いので p 値は参考。符号と大きさのみ見る。"},
        "ネコの扱い": ("ネコには経過日数を割り当てられないため時間Mantelには含めない。"
                       "距離の参照として results/time_axis_noIRI_cat_ref.csv に併記。"),
        "モデル内の 距離 vs 経過時間差": st.to_dict("records"),
        "対数関係は3モデルで成立するか": ("IRI のみ評価可能（9時点36ペア）。PodTRECK は3時点3ペアで"
                                          "参考値、UUO は2時点1ペアで評価不可。"
                                          "『3モデルすべてで成立』は現データでは検証できない。"),
        "図": "results/supplementary/within_model_time.png, results/supplementary/iri_trajectory.png",
    }, cfg)


if __name__ == "__main__":
    main()
