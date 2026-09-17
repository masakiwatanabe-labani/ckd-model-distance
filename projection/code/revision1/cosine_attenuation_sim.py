"""Task 2. 非中心化コサインに減衰補正式 sqrt(r_X r_Y) を適用してよいかを確認する。

古典的テスト理論の減衰式は Pearson 相関について導かれたものである。本稿の cos θ は
中心化していないので、同じ関係が厳密には成立しない。全遺伝子に共通する平均シフトは
非中心化コサインを膨らませるので、天井が上限として働かない可能性がある。

問いは「観測 cos θ が sqrt(r_half_X r_half_Y) を超えうるか」である。天井は本来
「2つの状態が同一の真の応答を持つときに観測されるはずの cos」なので、検査は

  (1) 真の応答が同一（true cos = 1）のとき、観測 cos が天井に一致するか（較正）
  (2) 真の応答が異なる（true cos < 1）とき、観測 cos が天井を超える頻度（上限性）

の2つになる。r_half も観測 cos と同じく非中心化コサインとして推定されるので、
平均シフトの影響は両者に乗る。したがって解析手順をそのまま模して標本から作る。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
OUT = HERE / "results" / "revision1"
OUT.mkdir(parents=True, exist_ok=True)
FIGDIR = HERE / "manuscript" / "figures"
P_GENES = 2016
N_REP = 300
N_INNER = 50
SEED = 20260826
TRUE_COS = [1.0, 0.9, 0.7, 0.5, 0.3]
NOISE = [0.3, 0.6, 1.0, 1.5, 2.2]      # 1 標本あたりの誤差 SD（信号 SD = 1）
DESIGNS = {"mouse-like (3 vs 6)": (3, 6), "feline-like (7 vs 6)": (7, 6)}


def observed_shift_range():
    """実データの各状態について |mean(Δ)| / sd(Δ)。非中心化の膨らみの大きさを決める。"""
    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    want = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split() if g.strip()}
    D = D.loc[D.index.intersection(sorted(want))].dropna(axis=0, how="any")
    return (D.mean(axis=0) / D.std(axis=0)).abs()


def true_pair(true_cos, shift, rng, p=P_GENES):
    """中心化成分の角度が true_cos になる 2 本に、全遺伝子共通の平均シフトを足す。"""
    z1 = rng.standard_normal(p)
    z2 = rng.standard_normal(p)
    z2 = z2 - (z2 @ z1) / (z1 @ z1) * z1
    z1 = z1 / z1.std()
    z2 = z2 / z2.std()
    ang = np.arccos(np.clip(true_cos, -1, 1))
    dx = z1 + shift
    dy = np.cos(ang) * z1 + np.sin(ang) * z2 + shift
    return dx, dy


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def one_state(delta, n_case, n_ctrl, sigma, rng, p=P_GENES):
    """標本を作り、全標本の Δ と分割半 r_half を返す（解析と同じ手順）。"""
    case = delta[:, None] + rng.standard_normal((p, n_case)) * sigma
    ctrl = rng.standard_normal((p, n_ctrl)) * sigma
    full = case.mean(1) - ctrl.mean(1)
    vals = np.empty(N_INNER)
    for i in range(N_INNER):
        a = rng.permutation(n_case); b = rng.permutation(n_ctrl)
        d1 = case[:, a[: n_case // 2]].mean(1) - ctrl[:, b[: n_ctrl // 2]].mean(1)
        d2 = case[:, a[n_case // 2:]].mean(1) - ctrl[:, b[n_ctrl // 2:]].mean(1)
        vals[i] = cosine(d1, d2)
    return full, float(np.nanmedian(vals))


def plot(T):
    """図 S1 を描く。シミュレーションを回し直さずに描き直せるよう分離してある。"""
    shifts = sorted(T["shift"].unique())
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from make_figures import S1, S2, S3, INK, INK2, FIGW, clean, savefig
    palette = [S1, S2, S3, INK2, "#7a5cc4"]
    fig, axes = plt.subplots(1, len(shifts), figsize=(FIGW, 2.55), sharey=True)
    for ax, shift in zip(np.atleast_1d(axes), shifts):
        e = T[(T["shift"] == shift) & (T.design == "mouse-like (3 vs 6)")]
        for col, tc in zip(palette, TRUE_COS):
            g = e[e.true_cos == tc].sort_values("mean_r_half")
            ax.plot(g.mean_r_half, g.frac_obs_above_ceiling, marker="o", markersize=4.0,
                    linewidth=1.5, color=col, label=f"{tc:.1f}")
        ax.axhline(0.01, color=INK, linestyle=(0, (4, 3)), linewidth=1.0)
        ax.set_title(f"mean shift = {shift:g}", loc="left", pad=4)
        ax.set_xlabel("estimated split-half reliability")
        clean(ax)
    np.atleast_1d(axes)[0].set_ylabel("fraction of replicates with\nobserved above ceiling")
    # 3 パネル共通の符号化なので、パネル内ではなく図下部に置く（曲線に重なるのを避ける）
    h, l = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(h, l, title=r"true $\cos\theta$", loc="lower center", ncol=5,
               frameon=False, handlelength=1.6, columnspacing=1.4,
               bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.10, 1, 1), w_pad=0.6)
    print("\nFigS1:", savefig(fig, FIGDIR / "FigS1_cosine_attenuation")[0].name)
    plt.close(fig)


def main():
    sh = observed_shift_range()
    print("実データの |mean(Δ)|/sd(Δ): 中央 %.3f, 範囲 [%.3f, %.3f]"
          % (sh.median(), sh.min(), sh.max()))
    shifts = [0.0, float(np.round(sh.median(), 2)), float(np.round(sh.max(), 2))]
    rng = np.random.default_rng(SEED)
    rows = []
    for dname, (n_case, n_ctrl) in DESIGNS.items():
        for shift in shifts:
            for sigma in NOISE:
                for tc in TRUE_COS:
                    obs = np.empty(N_REP); ceil = np.empty(N_REP); rel = np.empty(N_REP)
                    for i in range(N_REP):
                        dx, dy = true_pair(tc, shift, rng)
                        fx, rx = one_state(dx, n_case, n_ctrl, sigma, rng)
                        fy, ry = one_state(dy, n_case, n_ctrl, sigma, rng)
                        obs[i] = cosine(fx, fy)
                        ceil[i] = np.sqrt(max(rx, 0.0) * max(ry, 0.0))
                        rel[i] = 0.5 * (rx + ry)
                    rows.append({"design": dname, "shift": shift, "sigma": sigma,
                                 "true_cos": tc,
                                 "mean_r_half": round(float(np.nanmean(rel)), 4),
                                 "median_obs_cos": round(float(np.nanmedian(obs)), 4),
                                 "median_ceiling": round(float(np.nanmedian(ceil)), 4),
                                 "frac_obs_above_ceiling": float(np.nanmean(obs > ceil))})
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "cosine_attenuation_sim.tsv", sep="\t", index=False)

    print("\n(1) 較正: 真の応答が同一 (true cos = 1) のとき、観測 cos と天井")
    c = T[T.true_cos == 1.0]
    print(c.pivot_table(index=["design", "shift"], columns="sigma",
                        values="frac_obs_above_ceiling").round(3).to_string())
    print("\n(2) 上限性: true cos < 1 のとき観測 cos が天井を超えた割合")
    d = T[T.true_cos < 1.0]
    print(d.pivot_table(index=["design", "shift"], columns="true_cos",
                        values="frac_obs_above_ceiling").round(4).to_string())
    band = T[(T.mean_r_half.between(0.5, 0.9)) & (T["shift"] > 0) & (T.true_cos < 1.0)]
    print("\n実データの信頼性帯 (r_half 0.5-0.9) かつ平均シフトあり、true cos < 1:")
    print("  超過割合 最大 %.4f / 中央 %.4f （%d 条件）"
          % (band.frac_obs_above_ceiling.max(), band.frac_obs_above_ceiling.median(), len(band)))

    plot(T)
    print(f"書き出し: {OUT}")


if __name__ == "__main__":
    main()
