"""Part B-1 の続き。判断に直結する操作特性を測る。

被覆率は「区間が推定対象を覆う割合」だが、本文が使っているのは
「区間の下限がゼロを上回る ＝ 測定誤差だけでは説明できない」という判断である。
真の gap がゼロの条件を作り、その判断が誤って出る割合（偽陽性率）を測る。

真の gap = 0 になる true cos を、設計ごとに二分探索で決める。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from interval_coverage_sim import (CONFIGS, SHIFT, USE_SB, boot_interval,  # noqa: E402
                                   calibrate, gap_point, sample, true_pair)

OUT = HERE / "results" / "round5"
SEED = 20260917
N_THETA = int(os.environ.get("NR_THETA", 300))
N_SIM = int(os.environ.get("NR_SIM", 150))
N_BOOT = int(os.environ.get("NR_BOOT", 250))


def theta_of(tc, nca, nka, sa, ncb, nkb, sb, rng, n=N_THETA):
    dx, dy = true_pair(tc, SHIFT, rng)
    vals = [gap_point(*sample(dx, nca, nka, sa, rng), *sample(dy, ncb, nkb, sb, rng), rng)
            for _ in range(n)]
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.median(vals)), dx, dy


def main() -> int:
    rng = np.random.default_rng(SEED)
    rows = []
    for name, (nca, nka, tra), (ncb, nkb, trb) in CONFIGS:
        sa = calibrate(nca, nka, tra, rng)
        sb = calibrate(ncb, nkb, trb, rng)
        lo, hi = 0.3, 1.0                      # gap は true cos について単調減少
        for _ in range(12):
            mid = 0.5 * (lo + hi)
            th, _, _ = theta_of(mid, nca, nka, sa, ncb, nkb, sb, rng)
            if th > 0:
                lo = mid
            else:
                hi = mid
        tc0 = 0.5 * (lo + hi)
        th0, dx, dy = theta_of(tc0, nca, nka, sa, ncb, nkb, sb, rng, n=600)
        excl = usable = 0
        med = []
        for _ in range(N_SIM):
            ca, ka = sample(dx, nca, nka, sa, rng)
            cb, kb = sample(dy, ncb, nkb, sb, rng)
            l, h, _fu, m = boot_interval(ca, ka, cb, kb, rng, n_boot=N_BOOT)
            if not np.isfinite(l):
                continue
            usable += 1
            med.append(m)
            if l > 0:
                excl += 1
        rows.append(dict(config=name, true_cos_at_zero_gap=round(tc0, 4),
                         theta_gap=round(th0, 4), n_sim=usable,
                         false_positive_rate=round(excl / usable, 4) if usable else np.nan,
                         median_bootstrap_gap=round(float(np.median(med)), 4) if med else np.nan))
        print(f"{name}: 真の gap が 0 になる true cos = {tc0:.3f} (θ={th0:+.4f}) → "
              f"区間が 0 を上回る割合 {rows[-1]['false_positive_rate']:.4f}", flush=True)
    T = pd.DataFrame(rows)
    name = "interval_null_rate_SB.tsv" if USE_SB else "interval_null_rate.tsv"
    T.to_csv(OUT / name, sep="\t", index=False)
    print("\n偽陽性率: 平均 %.4f 最大 %.4f（名目 2.5%% 片側）"
          % (T.false_positive_rate.mean(), T.false_positive_rate.max()))
    print(f"\n基準: {'Spearman-Brown 補正' if USE_SB else '無補正'}")
    print(f"書き出し: {OUT / name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
