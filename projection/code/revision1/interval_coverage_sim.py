"""Part B-1. 個体ブートストラップ区間の被覆率を調べる。

§4.11 の手続き（復元抽出 → 異なる個体が 2 頭未満なら未定義 → 定義できた反復の
2.5/97.5 パーセンタイル）が、想定した 95% の被覆を持つかは未検証だった。
「95% 区間がゼロを除外した 24 組」を主要な根拠にしている以上、確かめる必要がある。

手順:
  1. 真の Δ を 2 本用意する（中心化成分の角度 = true_cos、全遺伝子共通の平均シフト
     を実データの中央値に合わせる）。
  2. 実データと同じ群サイズで標本を作る。ノイズ SD は、その設計での r_half が
     実データの値に合うように較正する。
  3. 推定対象 θ を決める: **真の Δ を固定したまま**個体だけを引き直して多数の
     データを作り、それぞれで gap = sqrt(r_half_X r_half_Y) − cos を計算した
     分布の中央値。個体ブートストラップが表せるのは個体の入れ替わりに対する
     不確実性だけなので、推定対象もそれに揃える（真の Δ まで引き直すと、
     区間が原理的に表せない変動を推定対象に含めてしまう）。
  4. データを N_SIM 組作り、各組で §4.11 の個体ブートストラップを回して
     2.5/97.5 パーセンタイル区間を作り、θ を含む割合を数える。
  5. 併せて「区間がゼロを上回る」割合と、未定義反復の割合を出す。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
OUT = HERE / "results" / "round5"
OUT.mkdir(parents=True, exist_ok=True)

P_GENES = 2016
SEED = 20260916
N_TRUTH = int(os.environ.get("CV_TRUTH", 600))
N_TRUTH_REAL = int(os.environ.get("CV_REAL", 4))
N_SIM = int(os.environ.get("CV_SIM", 400))
N_BOOT = int(os.environ.get("CV_BOOT", 400))
N_INNER = int(os.environ.get("CV_INNER", 12))
MIN_DISTINCT = 2
# CV_SB=1 で Spearman-Brown 補正済み基準を使う。SB は半標本の信頼性を全標本の信頼性に
# 外挿する式なので、基準が半標本・観測が全標本という不一致をまさに補正する操作にあたる。
USE_SB = os.environ.get("CV_SB", "0") == "1"

# 実データの猫 x マウス ペアに出てくる設計（results/reliability/reliability.tsv より）
CONFIGS = [
    # 名前, 猫 (n_case, n_ctrl, r_half 目標), マウス (n_case, n_ctrl, r_half 目標)
    ("cat CKD3/4 cortex x IRI (7v6, 3v6)", (7, 6, 0.841), (3, 6, 0.891)),
    ("cat CKD1/2 cortex x IRI (8v6, 3v6)", (8, 6, 0.568), (3, 6, 0.891)),
    ("cat CKD3/4 medulla x Pod-TRECK (5v6, 3v3)", (5, 6, 0.616), (3, 3, 0.908)),
    ("cat CKD1/2 medulla x IRI 2 h (7v6, 3v6)", (7, 6, 0.851), (3, 6, 0.677)),
]
TRUE_COS = [0.2, 0.35, 0.5]
SHIFT = 0.26          # 実データの状態の中央値 |mean|/sd


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def true_pair(true_cos, shift, rng, p=P_GENES):
    z1 = rng.standard_normal(p)
    z2 = rng.standard_normal(p)
    z2 = z2 - (z2 @ z1) / (z1 @ z1) * z1
    z1 = z1 / z1.std()
    z2 = z2 / z2.std()
    ang = np.arccos(np.clip(true_cos, -1, 1))
    return z1 + shift, np.cos(ang) * z1 + np.sin(ang) * z2 + shift


def sample(delta, n_case, n_ctrl, sigma, rng, p=P_GENES):
    case = delta[:, None] + rng.standard_normal((p, n_case)) * sigma
    ctrl = rng.standard_normal((p, n_ctrl)) * sigma
    return case, ctrl


def r_half(case, ctrl, rng, n_inner=N_INNER):
    """異なる列だけで分割半 cos の中央値。2 列未満なら NaN（§4.11 と同じ規則）。"""
    nc, nk = case.shape[1], ctrl.shape[1]
    if nc < MIN_DISTINCT or nk < MIN_DISTINCT:
        return np.nan
    vals = np.empty(n_inner)
    for i in range(n_inner):
        a = rng.permutation(nc); b = rng.permutation(nk)
        d1 = case[:, a[: nc // 2]].mean(1) - ctrl[:, b[: nk // 2]].mean(1)
        d2 = case[:, a[nc // 2:]].mean(1) - ctrl[:, b[nk // 2:]].mean(1)
        vals[i] = cosine(d1, d2)
    return float(np.nanmedian(vals))


def calibrate(n_case, n_ctrl, target, rng, p=P_GENES):
    """その設計で r_half が target になるノイズ SD を二分探索で決める。"""
    def rh(sigma):
        out = []
        for _ in range(40):
            d = rng.standard_normal(p)
            d = d / d.std() + SHIFT
            c, k = sample(d, n_case, n_ctrl, sigma, rng, p)
            out.append(r_half(c, k, rng, 10))
        return float(np.nanmedian(out))
    lo, hi = 0.02, 8.0
    for _ in range(26):
        mid = np.sqrt(lo * hi)
        if rh(mid) > target:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def benchmark(ra, rb):
    """2 つの r_half から基準値を作る。USE_SB なら Spearman-Brown 補正を挟む。"""
    a, b = max(ra, 0.0), max(rb, 0.0)
    if USE_SB:
        a = 2.0 * a / (1.0 + a) if a > 0 else 0.0
        b = 2.0 * b / (1.0 + b) if b > 0 else 0.0
    return float(np.sqrt(a * b))


def gap_point(ca, ka, cb, kb, rng):
    """1 データセットの点推定 gap = 基準値 − 観測 cos。"""
    da = ca.mean(1) - ka.mean(1)
    db = cb.mean(1) - kb.mean(1)
    ra, rb = r_half(ca, ka, rng), r_half(cb, kb, rng)
    if not (np.isfinite(ra) and np.isfinite(rb)):
        return np.nan
    return benchmark(ra, rb) - cosine(da, db)


def boot_interval(ca, ka, cb, kb, rng, n_boot=N_BOOT):
    """§4.11 の個体ブートストラップ。列 = 個体。未定義反復は捨てる。"""
    gaps = np.full(n_boot, np.nan)
    n_undef = 0
    na, nka, nb, nkb = ca.shape[1], ka.shape[1], cb.shape[1], kb.shape[1]
    for i in range(n_boot):
        ia = rng.integers(0, na, na); ika = rng.integers(0, nka, nka)
        ib = rng.integers(0, nb, nb); ikb = rng.integers(0, nkb, nkb)
        ua, uka = np.unique(ia), np.unique(ika)
        ub, ukb = np.unique(ib), np.unique(ikb)
        if min(len(ua), len(uka), len(ub), len(ukb)) < MIN_DISTINCT:
            n_undef += 1
            continue
        da = ca[:, ia].mean(1) - ka[:, ika].mean(1)
        db = cb[:, ib].mean(1) - kb[:, ikb].mean(1)
        ra = r_half(ca[:, ua], ka[:, uka], rng)
        rb = r_half(cb[:, ub], kb[:, ukb], rng)
        if not (np.isfinite(ra) and np.isfinite(rb)):
            n_undef += 1
            continue
        gaps[i] = benchmark(ra, rb) - cosine(da, db)
    g = gaps[np.isfinite(gaps)]
    if len(g) < 20:
        return np.nan, np.nan, n_undef / n_boot, np.nan
    return (float(np.percentile(g, 2.5)), float(np.percentile(g, 97.5)),
            n_undef / n_boot, float(np.median(g)))


def main() -> int:
    rng = np.random.default_rng(SEED)
    rows = []
    for name, (nca, nka, tra), (ncb, nkb, trb) in CONFIGS:
        sa = calibrate(nca, nka, tra, rng)
        sb = calibrate(ncb, nkb, trb, rng)
        print(f"{name}: ノイズ SD 猫 {sa:.3f} / マウス {sb:.3f}（目標 r_half {tra} / {trb}）",
              flush=True)
        for tc in TRUE_COS:
            cover = excl0 = usable = below = above = 0
            undef, widths, thetas = [], [], []
            # 真の Δ は固定し、個体だけを引き直す。個体ブートストラップが表せるのは
            # 個体の入れ替わりに対する不確実性だけなので、推定対象もそれに揃える。
            for _ in range(N_TRUTH_REAL):
                dx, dy = true_pair(tc, SHIFT, rng)
                truth = []
                for _ in range(N_TRUTH):
                    ca, ka = sample(dx, nca, nka, sa, rng)
                    cb, kb = sample(dy, ncb, nkb, sb, rng)
                    truth.append(gap_point(ca, ka, cb, kb, rng))
                truth = np.array([t for t in truth if np.isfinite(t)])
                theta = float(np.median(truth))
                thetas.append(theta)
                for _ in range(N_SIM // N_TRUTH_REAL):
                    ca, ka = sample(dx, nca, nka, sa, rng)
                    cb, kb = sample(dy, ncb, nkb, sb, rng)
                    lo, hi, fu, med = boot_interval(ca, ka, cb, kb, rng)
                    undef.append(fu)
                    if not np.isfinite(lo):
                        continue
                    usable += 1
                    widths.append(hi - lo)
                    if lo <= theta <= hi:
                        cover += 1
                    elif theta < lo:
                        below += 1
                    else:
                        above += 1
                    if lo > 0:
                        excl0 += 1
            rows.append(dict(
                config=name, true_cos=tc,
                theta_gap=round(float(np.mean(thetas)), 4),
                n_sim=usable,
                coverage_95=round(cover / usable, 4) if usable else np.nan,
                miss_theta_below_interval=round(below / usable, 4) if usable else np.nan,
                miss_theta_above_interval=round(above / usable, 4) if usable else np.nan,
                frac_interval_excludes_zero=round(excl0 / usable, 4) if usable else np.nan,
                mean_frac_undefined=round(float(np.mean(undef)), 4),
                mean_interval_width=round(float(np.mean(widths)), 4) if widths else np.nan))
            print(f"   true cos {tc}: θ={rows[-1]['theta_gap']:.3f}  "
                  f"被覆 {rows[-1]['coverage_95']:.3f}  "
                  f"未定義 {rows[-1]['mean_frac_undefined']:.3f}  "
                  f"幅 {rows[-1]['mean_interval_width']:.3f}  "
                  f"区間が 0 を上回る {rows[-1]['frac_interval_excludes_zero']:.3f}", flush=True)

    T = pd.DataFrame(rows)
    name = "interval_coverage_sim_SB.tsv" if USE_SB else "interval_coverage_sim.tsv"
    T.to_csv(OUT / name, sep="\t", index=False)
    print("\n全条件の被覆率: 平均 %.3f 中央 %.3f 最小 %.3f 最大 %.3f"
          % (T.coverage_95.mean(), T.coverage_95.median(),
             T.coverage_95.min(), T.coverage_95.max()))
    print("未定義反復の割合: 中央 %.3f [%.3f, %.3f]（実データは 0.098-0.235）"
          % (T.mean_frac_undefined.median(), T.mean_frac_undefined.min(),
             T.mean_frac_undefined.max()))
    print("θ が区間より下（区間が高すぎる）: 平均 %.3f / θ が区間より上（低すぎる）: 平均 %.3f"
          % (T.miss_theta_below_interval.mean(), T.miss_theta_above_interval.mean()))
    print("区間が 0 を上回る割合: 平均 %.3f（真の gap は全条件で正）"
          % T.frac_interval_excludes_zero.mean())
    print(f"\n基準: {'Spearman-Brown 補正' if USE_SB else '無補正'}")
    print(f"書き出し: {OUT / name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
