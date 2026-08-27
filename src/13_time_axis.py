"""時間軸を明示的な変数として扱い、距離を「経過日数」と「入口」で説明する。

12_model_distance.py で副産物として気づいた
「入口より時間軸の方が距離を支配しているように見える」を、正面から検定する。

方法:
  各マウス状態に傷害後経過日数を数値で与え、状態間の
      T = |Δ log10(days)|   時間差（対数）… 2日〜365日と2桁以上開くため対数
      E = 入口不一致 (0/1)  糸球体 vs 尿細管
      M = モデル不一致 (0/1) 参考
  を距離行列 D と突き合わせる。

  ペアワイズ距離は独立でない（1状態が複数ペアに現れる）ため、
  通常の相関検定の p は使えない。状態ラベルの並べ替えによる
  Mantel 検定・偏 Mantel 検定で評価する。

交絡に注意:
  糸球体入口は PodTRECK のみで 5-21日、尿細管入口は UUO(2-8日) と IRI(180-365日)。
  入口と時間は直交していないので、交絡の程度を必ず定量して併記する。
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("13_time_axis")
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m12", Path(__file__).parent / "12_model_distance.py")
m12 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(m12)
m09 = m12.m09

# 傷害後経過日数。PodTRECK は DT 誘導後、IRI は虚血再灌流後、UUO は結紮後。
DAYS = {"PodTRECK_5D": 5, "PodTRECK_2W": 14, "PodTRECK_3W": 21,
        "IRI_6mo": 180, "IRI_12mo": 365, "UUO_2D": 2, "UUO_8D": 8}
ENTRY = {"PodTRECK_5D": "glomerular", "PodTRECK_2W": "glomerular", "PodTRECK_3W": "glomerular",
         "IRI_6mo": "tubular", "IRI_12mo": "tubular", "UUO_2D": "tubular", "UUO_8D": "tubular"}


def mantel(D, X, n_perm=9999, seed=0, partial=None):
    """Mantel（または偏 Mantel）検定。D, X, partial は対称行列（正方 ndarray）。

    偏 Mantel は D と X をそれぞれ partial で回帰した残差の相関を取る。
    """
    rng = np.random.default_rng(seed)
    k = D.shape[0]
    iu = np.triu_indices(k, 1)

    def vec(M):
        return M[iu]

    def resid(y, z):
        A = np.column_stack([np.ones_like(z), z])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        return y - A @ beta

    d, x = vec(D), vec(X)
    if partial is not None:
        z = vec(partial)
        d, x = resid(d, z), resid(x, z)
    r_obs = stats.spearmanr(d, x)[0]

    cnt = 0
    for _ in range(n_perm):
        p = rng.permutation(k)
        Dp = D[np.ix_(p, p)]
        dp = vec(Dp)
        if partial is not None:
            dp = resid(dp, vec(partial))
        rp = stats.spearmanr(dp, x)[0]
        if abs(rp) >= abs(r_obs) - 1e-12:
            cnt += 1
    return r_obs, (cnt + 1) / (n_perm + 1)


def main():
    S, model = m12.build_states()
    states = [s for s in DAYS if s in S]
    k = len(states)
    log.info("状態 %d: %s", k, states)

    D = np.zeros((k, k))
    for i, j in combinations(range(k), 2):
        d, _ = m12.dist(S[states[i]], S[states[j]])
        D[i, j] = D[j, i] = d

    T = np.zeros((k, k))   # |Δ log10 days|
    Tlin = np.zeros((k, k))
    E = np.zeros((k, k))
    M = np.zeros((k, k))
    for i, j in combinations(range(k), 2):
        a, b = states[i], states[j]
        T[i, j] = T[j, i] = abs(np.log10(DAYS[a]) - np.log10(DAYS[b]))
        Tlin[i, j] = Tlin[j, i] = abs(DAYS[a] - DAYS[b])
        E[i, j] = E[j, i] = float(ENTRY[a] != ENTRY[b])
        M[i, j] = M[j, i] = float(model[a] != model[b])

    iu = np.triu_indices(k, 1)
    pairs = pd.DataFrame({
        "a": [states[i] for i in iu[0]], "b": [states[j] for j in iu[1]],
        "distance": D[iu], "dlog10days": T[iu], "ddays": Tlin[iu],
        "entry_mismatch": E[iu].astype(int), "model_mismatch": M[iu].astype(int)})
    pairs.to_csv(RES / "time_axis_pairs.csv", index=False)

    # 交絡の程度
    conf = stats.spearmanr(T[iu], E[iu])[0]
    log.info("交絡: |Δlog10日数| と 入口不一致 の相関 rho=%.3f", conf)

    res = {}
    for name, X in [("時間差 |Δlog10(日数)|", T), ("時間差 |Δ日数|(線形)", Tlin),
                    ("入口不一致", E), ("モデル不一致", M)]:
        r, p = mantel(D, X, seed=SEED)
        res[name] = {"mantel_rho": round(float(r), 3), "p": round(float(p), 4)}
        log.info("Mantel %-24s rho=%+.3f p=%.4f", name, r, p)

    # 偏 Mantel（互いを打ち消す）
    r_te, p_te = mantel(D, T, seed=SEED, partial=E)
    r_et, p_et = mantel(D, E, seed=SEED, partial=T)
    res["時間差 | 入口を統制"] = {"mantel_rho": round(float(r_te), 3), "p": round(float(p_te), 4)}
    res["入口 | 時間差を統制"] = {"mantel_rho": round(float(r_et), 3), "p": round(float(p_et), 4)}
    log.info("偏Mantel 時間|入口統制 rho=%+.3f p=%.4f", r_te, p_te)
    log.info("偏Mantel 入口|時間統制 rho=%+.3f p=%.4f", r_et, p_et)

    verdict = ("時間軸が距離をより説明する" if abs(r_te) > abs(r_et)
               else "入口が距離をより説明する")
    log.info("判定: %s", verdict)

    # ---- ネコを時間軸に載せられるか（経験的チェック）----
    cat_e = m09.to_human_space(m09.internal_lfc("cat_rna_ctx_early"), "cat")
    cat_l = m09.to_human_space(m09.internal_lfc("cat_rna_ctx_late"), "cat")
    catrows = []
    for lab, cs in [("ネコ皮質 早期(IRIS1/2)", cat_e), ("ネコ皮質 晩期(IRIS3/4)", cat_l)]:
        ds = {s: m12.dist(cs, S[s])[0] for s in states}
        best = min(ds, key=ds.get)
        # 距離の逆数で重みづけした「見かけの経過日数」
        w = np.array([1 / ds[s] for s in states])
        wd = float(np.exp(np.sum(w * np.log([DAYS[s] for s in states])) / w.sum()))
        catrows.append({"cat_state": lab, "最近傍マウス状態": best,
                        "その距離": round(ds[best], 4),
                        "距離重みづけ見かけ日数(幾何平均)": round(wd, 1),
                        **{f"d_{s}": round(ds[s], 4) for s in states}})
    catdf = pd.DataFrame(catrows)
    catdf.to_csv(RES / "cat_on_time_axis.csv", index=False)
    log.info("\n%s", catdf[["cat_state", "最近傍マウス状態", "その距離",
                             "距離重みづけ見かけ日数(幾何平均)"]].to_string(index=False))

    append_summary("13_time_axis / 距離を説明するのは時間軸か入口か", {
        "経過日数の割り当て": DAYS,
        "入口の割り当て": ENTRY,
        "Mantel検定（状態ラベル並べ替え, n_perm=9999）": res,
        "判定": verdict,
        "交絡の程度": (f"|Δlog10日数| と 入口不一致 の相関 rho={conf:.3f}。事前には交絡を懸念したが、"
                       "実測ではほぼ直交していた。理由: 尿細管入口が最短(UUO 2-8日)と"
                       "最長(IRI 180-365日)の両端を占め、糸球体入口(PodTRECK 5-21日)が"
                       "中間に入るため、時間と入口が分離している。"
                       "したがって偏Mantel の値は素直に解釈できる。"),
        "ネコを時間軸に載せた場合": catdf.to_dict("records"),
        "ネコは経過時間の代理に使えるか": (
            "使えない。理由: (1) 自然発症で明確な傷害イベントが無く、発症時点が不明。"
            "マウスの『傷害後n日』に対応する原点が定義できない。"
            "(2) 横断デザインで、各病期は別個体。同一個体の時間推移を見ていない。"
            "(3) IRISステージは血清クレアチニン（腎機能）による定義であって時間の定義ではない。"
            "残存ネフロン量の指標であり、進行速度は個体差が大きい。"
            "(4) 実測では、病期の進行が時間軸上の前進として現れない。むしろ逆向きに動く: "
            "早期(IRIS1/2)の最近傍は IRI_12mo（365日）だが、晩期(IRIS3/4)の最近傍は "
            "PodTRECK_2W（14日）。IRI 6mo/12mo への距離は早期→晩期でむしろ増加する "
            "(0.536→0.624, 0.460→0.535)。病期が進むほど慢性マウス状態から離れる。"
            "→ 病期は『進行度』の順序尺度としては使えるが、経過日数の代理には使えない。"),
        "制約": ("状態7個・ペア21個と小さい。Mantel の検出力は低く、"
                 "rho の点推定も不安定。方向性の議論に留めること。"),
    }, cfg)


if __name__ == "__main__":
    main()
