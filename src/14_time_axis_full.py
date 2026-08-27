"""GSE98622 の全時点を個別状態にして、時間軸 vs 入口を再検定する。

13_time_axis.py は状態7個・ペア21個と小さかった。IRI の全10時点を個別状態に
すると 14-15 状態・91-105 ペアになり検出力が上がる。IRI 内は同一データセット
なのでバッチ交絡が入らない。

対照の取り方（結果を左右するので明示する）:
  GPL13112 の sham は SHAM4h(3), SHAM24h(3), SHAM12m(3) の計9。
  **全shamをプールしてはいけない。** SHAM12m は12か月齢、SHAM4h/24h は若齢で、
  プールすると IRI12m だけが「加齢 + 傷害」対「若齢」の比較になり、
  加齢シグナルが時間差として現れて時間効果を過大評価する（仮説に有利な方向のバイアス）。
  したがって:
      IRI 2h..28d -> SHAM4h + SHAM24h (n=6, 若齢)
      IRI 12mo    -> SHAM12m (n=3, 同週齢)
      IRI 6mo     -> NORM3m/9m/15m (n=9, GPL19057。別プラットフォーム)
  全shamプール版も感度解析として出し、選択の影響を定量する。

IRI 6mo は別プラットフォーム(GPL19057)なので、含める版と除く版の両方を出す。
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

log = get_logger("14_time_axis_full")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m12", Path(__file__).parent / "12_model_distance.py")
m12 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m12)
m09 = m12.m09

YOUNG_SHAM = ["SHAM4h", "SHAM24h"]
IRI_TIMEPOINTS = {          # ラベル: (経過日数, 対照ラベル群, プラットフォーム)
    "IRI2h":  (2 / 24,  YOUNG_SHAM, "GPL13112"),
    "IRI4h":  (4 / 24,  YOUNG_SHAM, "GPL13112"),
    "IRI24h": (1.0,     YOUNG_SHAM, "GPL13112"),
    "IRI48h": (2.0,     YOUNG_SHAM, "GPL13112"),
    "IRI72h": (3.0,     YOUNG_SHAM, "GPL13112"),
    "IRI7d":  (7.0,     YOUNG_SHAM, "GPL13112"),
    "IRI14d": (14.0,    YOUNG_SHAM, "GPL13112"),
    "IRI28d": (28.0,    YOUNG_SHAM, "GPL13112"),
    "IRI12m": (365.0,   ["SHAM12m"], "GPL13112"),
    "IRI6mN": (180.0,   ["NORM3m", "NORM9m", "NORM15m"], "GPL19057"),
}
OTHER = {"PodTRECK_5D": (5.0, "PodTRECK", "glomerular"),
         "PodTRECK_2W": (14.0, "PodTRECK", "glomerular"),
         "PodTRECK_3W": (21.0, "PodTRECK", "glomerular"),
         "UUO_2D": (2.0, "UUO", "tubular"),
         "UUO_8D": (8.0, "UUO", "tubular")}


def iri_states(pool_all_sham=False):
    mat = pd.read_parquet(INT / "mouse_iri_matrix.parquet")
    grp = pd.read_parquet(INT / "mouse_iri_grp.parquet").squeeze()
    out = {}
    for lab, (days, ctrl, plat) in IRI_TIMEPOINTS.items():
        c = ctrl
        if pool_all_sham and plat == "GPL13112":
            c = ["SHAM4h", "SHAM24h", "SHAM12m"]
        A = mat[grp[grp == lab].index]
        B = mat[grp[grp.isin(c)].index]
        A, B = np.log2(A + 1), np.log2(B + 1)
        keep = pd.concat([A, B], axis=1).max(axis=1) > 1
        res = moderated_ttest(A[keep], B[keep])
        out[lab] = m09.to_human_space(res["lfc"], "mouse")
    return out


def build(include_6mo=True, pool_all_sham=False):
    S, days, model, entry = {}, {}, {}, {}
    for k, ds in [("PodTRECK_5D", "m_rna_5d"), ("PodTRECK_2W", "m_rna_2w"),
                  ("PodTRECK_3W", "m_rna_3w")]:
        S[k] = m09.to_human_space(m09.internal_lfc(ds), "mouse")
    for k in ["UUO_2D", "UUO_8D"]:
        v = m12.uuo_lfc(day=k.split("_")[1])
        S[k] = m09.to_human_space(v, "mouse")
    for k, v in OTHER.items():
        days[k], model[k], entry[k] = v
    for k, v in iri_states(pool_all_sham).items():
        if k == "IRI6mN" and not include_6mo:
            continue
        S[k] = v
        days[k], model[k], entry[k] = IRI_TIMEPOINTS[k][0], "IRI", "tubular"
    return S, days, model, entry


def matrices(S, days, model, entry):
    states = list(S)
    k = len(states)
    D = np.zeros((k, k)); T = np.zeros((k, k)); E = np.zeros((k, k)); M = np.zeros((k, k))
    for i, j in combinations(range(k), 2):
        a, b = states[i], states[j]
        d, _ = m12.dist(S[a], S[b])
        D[i, j] = D[j, i] = d
        T[i, j] = T[j, i] = abs(np.log10(days[a]) - np.log10(days[b]))
        E[i, j] = E[j, i] = float(entry[a] != entry[b])
        M[i, j] = M[j, i] = float(model[a] != model[b])
    return states, D, T, E, M


def mantel(D, X, n_perm=9999, seed=0, partial=None, blocks=None):
    """Mantel / 偏Mantel。blocks を与えるとブロック内に限定した制限付き並べ替え。

    blocks（モデル名の配列）を使うと、モデル同一性を保ったまま状態を入れ替えるので、
    「モデルの違いでは説明できない時間の効果」を検定できる。
    """
    rng = np.random.default_rng(seed)
    k = D.shape[0]
    iu = np.triu_indices(k, 1)
    vec = lambda M: M[iu]  # noqa: E731

    def resid(y, z):
        A = np.column_stack([np.ones_like(z), z])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        return y - A @ beta

    d, x = vec(D), vec(X)
    if partial is not None:
        z = vec(partial); d, x = resid(d, z), resid(x, z)
    r_obs = stats.spearmanr(d, x)[0]

    idx_by_block = None
    if blocks is not None:
        blocks = np.asarray(blocks)
        idx_by_block = [np.where(blocks == b)[0] for b in pd.unique(blocks)]

    cnt = 0
    for _ in range(n_perm):
        if idx_by_block is None:
            p = rng.permutation(k)
        else:
            p = np.arange(k)
            for ib in idx_by_block:
                p[ib] = rng.permutation(ib)
        dp = vec(D[np.ix_(p, p)])
        if partial is not None:
            dp = resid(dp, vec(partial))
        if abs(stats.spearmanr(dp, x)[0]) >= abs(r_obs) - 1e-12:
            cnt += 1
    return float(r_obs), (cnt + 1) / (n_perm + 1)


def run(label, include_6mo, pool_all_sham=False):
    S, days, model, entry = build(include_6mo, pool_all_sham)
    states, D, T, E, M = matrices(S, days, model, entry)
    iu = np.triu_indices(len(states), 1)
    conf = stats.spearmanr(T[iu], E[iu])[0]
    log.info("[%s] 状態 %d / ペア %d、交絡 rho(T,E)=%.3f",
             label, len(states), len(iu[0]), conf)

    out = {"状態数": len(states), "ペア数": int(len(iu[0])),
           "交絡 rho(時間, 入口)": round(float(conf), 3)}
    for name, X in [("時間差 |Δlog10(日数)|", T), ("入口不一致", E), ("モデル不一致", M)]:
        r, p = mantel(D, X, seed=SEED)
        out[name] = {"rho": round(r, 3), "p": round(p, 4)}
        log.info("  Mantel %-22s rho=%+.3f p=%.4f", name, r, p)
    r, p = mantel(D, T, seed=SEED, partial=E)
    out["時間 | 入口を統制"] = {"rho": round(r, 3), "p": round(p, 4)}
    log.info("  偏Mantel 時間|入口   rho=%+.3f p=%.4f", r, p)
    r2, p2 = mantel(D, E, seed=SEED, partial=T)
    out["入口 | 時間を統制"] = {"rho": round(r2, 3), "p": round(p2, 4)}
    log.info("  偏Mantel 入口|時間   rho=%+.3f p=%.4f", r2, p2)
    rb, pb = mantel(D, T, seed=SEED, blocks=[model[s] for s in states])
    out["時間（モデル内制限並べ替え）"] = {"rho": round(rb, 3), "p": round(pb, 4)}
    log.info("  制限Mantel 時間(モデル保持) rho=%+.3f p=%.4f", rb, pb)
    out["判定"] = ("時間軸が距離をより説明する" if abs(out["時間 | 入口を統制"]["rho"])
                   > abs(out["入口 | 時間を統制"]["rho"]) else "入口が距離をより説明する")

    pd.DataFrame({"a": [states[i] for i in iu[0]], "b": [states[j] for j in iu[1]],
                  "distance": D[iu], "dlog10days": T[iu],
                  "entry_mismatch": E[iu].astype(int),
                  "model_mismatch": M[iu].astype(int)}).to_csv(
        RES / f"time_axis_full_pairs_{label}.csv", index=False)
    return out


def main():
    res = {}
    res["A: IRI6mo込み(15状態)"] = run("with6mo", True)
    res["B: IRI6mo除く(14状態, 全て同一プラットフォーム内)"] = run("no6mo", False)
    res["C: 感度 全shamプール(IRI6mo除く)"] = run("no6mo_poolsham", False, pool_all_sham=True)

    append_summary("14_time_axis_full / GSE98622全時点での時間軸 vs 入口", {
        **{k: v for k, v in res.items()},
        "対照の取り方": ("IRI 2h-28d は若齢sham(SHAM4h+SHAM24h, n=6)、IRI 12mo は同週齢"
                         "SHAM12m(n=3)、IRI 6mo は NORM3m/9m/15m(n=9, GPL19057)。"
                         "全shamをプールすると IRI12m だけが加齢+傷害 vs 若齢の比較になり、"
                         "加齢シグナルが時間差として混入して時間効果を過大評価する"
                         "（仮説に有利な方向のバイアス）ため、主解析では避けた。"
                         "感度解析Cがその影響の大きさを示す。"),
        "版の使い分け": ("版Aは IRI6mo を含むが GPL19057 との跨ぎがある。"
                         "版Bは IRI6mo を除き、全状態が各データセット内コントラストのみで、"
                         "IRI 9時点は完全に同一プラットフォーム。主結果は版Bで述べること。"),
        "制限並べ替えの意味": ("モデル同一性を保ったまま状態を入れ替える検定。"
                               "『モデルの違いでは説明できない、時間そのものの効果』を測る。"),
    }, cfg)


if __name__ == "__main__":
    main()
