"""Results 3.3 用の確認。

(1) 入口Mantel を「ネコ込み(16状態)」と「マウスのみ(14状態)」の両方で出す。
    注: 13-15番の時間軸解析にネコは一度も入っていない（14状態は全てマウス）。
    ネコには経過日数を割り当てられないため、ネコ込み版では
    時間を統制した偏Mantel は原理的に計算できない。入口・モデルのみ。

(2) 「粗いサンプリングのときだけモデル効果が見える」の検証。
    7状態版は PodTRECK3 + IRI2 + UUO2 だった。PodTRECK と UUO には
    予備の時点が無いので、変えられるのは「IRI 9時点からどの2点を選ぶか」だけ。
    C(9,2)=36 通りを**全数列挙**して、モデル不一致 rho の分布を出す。
    元の選択（6mo+12mo 相当＝ここでは 12m を含む組）が分布のどこかを見る。
    さらに IRI を k=2,4,6,9 点に増やしたときの rho の推移を見る。
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

log = get_logger("16_entry_model_checks")
cfg = load_config()
ROOT = Path(cfg["_root"])
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]

import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("m14", Path(__file__).parent / "14_time_axis_full.py")
m14 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m14)
m12, m09 = m14.m12, m14.m09

S, DAYS, MODEL, ENTRY = m14.build(include_6mo=False)
STATES = list(S)
# ネコ（経過日数なし。入口は尿細管）
CAT = {"cat_natural_ctx": m09.to_human_space(m09.internal_lfc("cat_rna_ctx_late"), "cat"),
       "cat_natural_med": m09.to_human_space(m09.internal_lfc("cat_rna_med_late"), "cat")}

_DCACHE: dict = {}


def d(a, b, pool):
    k = tuple(sorted([a, b]))
    if k not in _DCACHE:
        _DCACHE[k] = m12.dist(pool[a], pool[b])[0]
    return _DCACHE[k]


def build_mats(states, pool, model, entry, days=None):
    k = len(states)
    D = np.zeros((k, k)); E = np.zeros((k, k)); M = np.zeros((k, k))
    T = np.zeros((k, k)) if days else None
    for i, j in combinations(range(k), 2):
        a, b = states[i], states[j]
        D[i, j] = D[j, i] = d(a, b, pool)
        E[i, j] = E[j, i] = float(entry[a] != entry[b])
        M[i, j] = M[j, i] = float(model[a] != model[b])
        if days:
            T[i, j] = T[j, i] = abs(np.log10(days[a]) - np.log10(days[b]))
    return D, E, M, T


def task1():
    pool = {**S, **CAT}
    model_c = {**MODEL, "cat_natural_ctx": "cat_natural", "cat_natural_med": "cat_natural"}
    entry_c = {**ENTRY, "cat_natural_ctx": "tubular", "cat_natural_med": "tubular"}

    out = {}
    # --- マウスのみ 14状態 ---
    D, E, M, T = build_mats(STATES, pool, model_c, entry_c, DAYS)
    n_pairs = len(np.triu_indices(len(STATES), 1)[0])
    r_e, p_e = m14.mantel(D, E, seed=SEED)
    r_m, p_m = m14.mantel(D, M, seed=SEED)
    r_et, p_et = m14.mantel(D, E, seed=SEED, partial=T)
    r_mt, p_mt = m14.mantel(D, M, seed=SEED, partial=T)
    out["マウスのみ(14状態)"] = {
        "状態数": len(STATES), "ペア数": n_pairs,
        "入口不一致": {"rho": round(r_e, 3), "p": round(p_e, 4)},
        "モデル不一致": {"rho": round(r_m, 3), "p": round(p_m, 4)},
        "入口|時間を統制": {"rho": round(r_et, 3), "p": round(p_et, 4)},
        "モデル|時間を統制": {"rho": round(r_mt, 3), "p": round(p_mt, 4)}}

    # --- ネコ込み 16状態（時間を持たないので偏Mantelは不可）---
    st2 = STATES + list(CAT)
    D2, E2, M2, _ = build_mats(st2, pool, model_c, entry_c, None)
    n_pairs2 = len(np.triu_indices(len(st2), 1)[0])
    r_e2, p_e2 = m14.mantel(D2, E2, seed=SEED)
    r_m2, p_m2 = m14.mantel(D2, M2, seed=SEED)
    out["ネコ込み(16状態)"] = {
        "状態数": len(st2), "ペア数": n_pairs2,
        "入口不一致": {"rho": round(r_e2, 3), "p": round(p_e2, 4)},
        "モデル不一致": {"rho": round(r_m2, 3), "p": round(p_m2, 4)},
        "偏Mantel": "計算不可（ネコに経過日数を割り当てられないため時間行列が作れない）"}

    for k, v in out.items():
        log.info("[%s] %s", k, v)
    return out


def task2():
    iri = [s for s in STATES if MODEL[s] == "IRI"]
    fixed = [s for s in STATES if MODEL[s] != "IRI"]   # PodTRECK3 + UUO2
    pool = S
    rows = []
    for k in [2, 4, 6, 9]:
        combos = list(combinations(iri, k))
        rng = np.random.default_rng(SEED)
        if len(combos) > 200:
            idx = rng.choice(len(combos), 200, replace=False)
            combos = [combos[i] for i in idx]
        for c in combos:
            st = fixed + list(c)
            D, E, M, T = build_mats(st, pool, MODEL, ENTRY, DAYS)
            iu = np.triu_indices(len(st), 1)
            rows.append({"n_iri": k, "n_states": len(st), "n_pairs": len(iu[0]),
                         "iri_set": "+".join(sorted(c, key=lambda s: DAYS[s])),
                         "rho_model": stats.spearmanr(D[iu], M[iu])[0],
                         "rho_entry": stats.spearmanr(D[iu], E[iu])[0],
                         "rho_time": stats.spearmanr(D[iu], T[iu])[0]})
    t = pd.DataFrame(rows)
    t.to_csv(RES / "subsample_model_effect.csv", index=False)

    summ = (t.groupby("n_states")
             .agg(n_combos=("rho_model", "size"),
                  model_med=("rho_model", "median"), model_min=("rho_model", "min"),
                  model_max=("rho_model", "max"),
                  model_frac_pos=("rho_model", lambda x: float((x > 0.2).mean())),
                  time_med=("rho_time", "median"), entry_med=("rho_entry", "median"))
             .round(3).reset_index())
    log.info("\n%s", summ.to_string(index=False))

    # 7状態(=IRI2点)のうち、12mを含む組 vs 含まない組
    k2 = t[t.n_iri == 2].copy()
    k2["has_12m"] = k2.iri_set.str.contains("IRI12m")
    cmp = k2.groupby("has_12m").rho_model.agg(["size", "median", "min", "max"]).round(3)
    log.info("\n7状態版: IRI12m を含むか別の rho_model\n%s", cmp.to_string())
    return t, summ, cmp


def main():
    o1 = task1()
    t, summ, cmp = task2()
    append_summary("16_entry_model_checks / Results3.3の確認", {
        "前提の訂正": ("13-15番の時間軸解析（14状態）は全てマウスで、ネコは含まれていない。"
                       "『14状態にネコが tubular として含まれる』は誤り。"
                       "ネコが入るのは09の距離行列のみ。"),
        "(1) 入口・モデルMantel": o1,
        "(2) 部分抽出の要約（IRI時点数を変えたとき）": summ.to_dict("records"),
        "(2) 7状態版で IRI12m を含むか別": cmp.to_dict(),
        "備考": ("7状態版で変えられるのは IRI 9時点から2点を選ぶ組み合わせのみ"
                 "（PodTRECK と UUO は予備時点が無い）。C(9,2)=36通りを全数列挙した。"),
    }, cfg)


if __name__ == "__main__":
    main()
