"""Part A-1. 天井超過を「標本数の差」と「非中心化」に分解する。

基準値 C = sqrt(r_half_X r_half_Y) は半標本から推定した量だが、比較する観測値は
全標本で計算している。半標本と全標本では測定誤差の大きさが違うので、観測値が C を
超えることは非中心化とは独立に起こりうる。旧 cosine_attenuation_sim.py は
非中心化の影響しか見ておらず、この2つを分離していない。

3条件を同じ標本の上で評価する。

  (i)   中心化ピアソン相関 × 標本数を揃えた比較
        基準も観測も半標本。理論どおりならほとんど超えないはず。ベースライン。
  (ii)  中心化ピアソン相関 × 半標本基準 vs 全標本観測
        (ii) - (i) が標本数の差だけの寄与。
  (iii) 非中心化コサイン × 半標本基準 vs 全標本観測
        本稿の無補正基準。(iii) - (ii) が非中心化の寄与。
  (iv)  中心化ピアソン × Spearman-Brown 補正基準 vs 全標本観測
        SB 補正は半標本の信頼性を全標本の信頼性に外挿する式なので、(ii) の標本数の
        不一致をまさに補正する操作にあたる。(iv) が (i) に近ければ較正されている。
  (v)   非中心化コサイン × Spearman-Brown 補正基準 vs 全標本観測
        本稿が Table 1 に併記している SB 基準そのもの。(v) - (iv) が非中心化の寄与。

信頼性の範囲は実データの全 state を包含するように取る（r_half 0.568-0.986、
results/reliability/reliability.tsv）。旧版の 0.50-0.84 は Table 1 の
同種・異データセット群（無補正基準の中央値 0.896）を外していた。

真の cos θ と観測 cos θ は別物なので、両方を出力する（Part A-3）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
OUT = HERE / "results" / "round5"
OUT.mkdir(parents=True, exist_ok=True)

P_GENES = 2016
N_REP = 300
N_INNER = 50
SEED = 20260916
TRUE_COS = [1.0, 0.9, 0.7, 0.5, 0.3]
# 実データの r_half は 0.568-0.986。その全域を覆うように誤差 SD を取る。
NOISE = [0.15, 0.25, 0.35, 0.5, 0.7, 1.0, 1.5, 2.2]
DESIGNS = {"mouse-like (3 vs 6)": (3, 6), "feline-like (7 vs 6)": (7, 6)}


def observed_shift_range():
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
    return z1 + shift, np.cos(ang) * z1 + np.sin(ang) * z2 + shift


def cosine(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else np.nan


def pearson(u, v):
    return cosine(u - u.mean(), v - v.mean())


def draw(delta, n_case, n_ctrl, sigma, rng, p=P_GENES):
    case = delta[:, None] + rng.standard_normal((p, n_case)) * sigma
    ctrl = rng.standard_normal((p, n_ctrl)) * sigma
    return case, ctrl


def halves(case, ctrl, rng):
    n_case, n_ctrl = case.shape[1], ctrl.shape[1]
    a = rng.permutation(n_case)
    b = rng.permutation(n_ctrl)
    d1 = case[:, a[: n_case // 2]].mean(1) - ctrl[:, b[: n_ctrl // 2]].mean(1)
    d2 = case[:, a[n_case // 2:]].mean(1) - ctrl[:, b[n_ctrl // 2:]].mean(1)
    return d1, d2


def one_cell(dx, dy, n_case, n_ctrl, sigma, rng):
    """1 反復。3 条件それぞれの (基準値, 観測値) を返す。"""
    cx, tx = draw(dx, n_case, n_ctrl, sigma, rng)
    cy, ty = draw(dy, n_case, n_ctrl, sigma, rng)
    full_x = cx.mean(1) - tx.mean(1)
    full_y = cy.mean(1) - ty.mean(1)

    rP = {"x": np.empty(N_INNER), "y": np.empty(N_INNER)}
    rC = {"x": np.empty(N_INNER), "y": np.empty(N_INNER)}
    matchP = np.empty(N_INNER)
    matchC = np.empty(N_INNER)
    for i in range(N_INNER):
        x1, x2 = halves(cx, tx, rng)
        y1, y2 = halves(cy, ty, rng)
        rP["x"][i], rP["y"][i] = pearson(x1, x2), pearson(y1, y2)
        rC["x"][i], rC["y"][i] = cosine(x1, x2), cosine(y1, y2)
        matchP[i] = pearson(x1, y1)          # 基準と同じ標本数どうし
        matchC[i] = cosine(x1, y1)

    def ceil(r):
        return np.sqrt(max(float(np.nanmedian(r["x"])), 0.0) * max(float(np.nanmedian(r["y"])), 0.0))

    def ceil_sb(r):
        """Spearman-Brown: r_XX = 2 r_half / (1 + r_half)、負は 0 に落とす。"""
        out = []
        for k in ("x", "y"):
            rh = max(float(np.nanmedian(r[k])), 0.0)
            out.append(2.0 * rh / (1.0 + rh) if rh > 0 else 0.0)
        return float(np.sqrt(out[0] * out[1]))

    cP, cC = ceil(rP), ceil(rC)
    return {
        "ceil_pearson": cP, "ceil_cosine": cC,
        "ceil_pearson_sb": ceil_sb(rP), "ceil_cosine_sb": ceil_sb(rC),
        "obs_matched_pearson": float(np.nanmedian(matchP)),
        "obs_full_pearson": pearson(full_x, full_y),
        "obs_full_cosine": cosine(full_x, full_y),
        "r_half_pearson": 0.5 * (float(np.nanmedian(rP["x"])) + float(np.nanmedian(rP["y"]))),
        "r_half_cosine": 0.5 * (float(np.nanmedian(rC["x"])) + float(np.nanmedian(rC["y"]))),
    }


def main() -> int:
    shifts_obs = observed_shift_range()
    shifts = [0.0, round(float(shifts_obs.median()), 2), round(float(shifts_obs.max()), 2)]
    print("実データの平均シフト |mean|/sd: 中央 %.3f 最大 %.3f  → 使う値 %s"
          % (shifts_obs.median(), shifts_obs.max(), shifts))

    rows = []
    rng = np.random.default_rng(SEED)
    for dname, (n_case, n_ctrl) in DESIGNS.items():
        for shift in shifts:
            for sigma in NOISE:
                for tc in TRUE_COS:
                    acc = {k: [] for k in ("i", "ii", "iii", "iv", "v", "rP", "rC",
                                           "obsP", "obsC", "obsM")}
                    for _ in range(N_REP):
                        dx, dy = true_pair(tc, shift, rng)
                        c = one_cell(dx, dy, n_case, n_ctrl, sigma, rng)
                        acc["i"].append(c["obs_matched_pearson"] > c["ceil_pearson"])
                        acc["ii"].append(c["obs_full_pearson"] > c["ceil_pearson"])
                        acc["iii"].append(c["obs_full_cosine"] > c["ceil_cosine"])
                        acc["iv"].append(c["obs_full_pearson"] > c["ceil_pearson_sb"])
                        acc["v"].append(c["obs_full_cosine"] > c["ceil_cosine_sb"])
                        acc["rP"].append(c["r_half_pearson"])
                        acc["rC"].append(c["r_half_cosine"])
                        acc["obsP"].append(c["obs_full_pearson"])
                        acc["obsC"].append(c["obs_full_cosine"])
                        acc["obsM"].append(c["obs_matched_pearson"])
                    rows.append({
                        "design": dname, "shift": shift, "sigma": sigma, "true_cos": tc,
                        "mean_r_half_pearson": round(float(np.nanmean(acc["rP"])), 4),
                        "mean_r_half_cosine": round(float(np.nanmean(acc["rC"])), 4),
                        "median_obs_cos_full": round(float(np.nanmedian(acc["obsC"])), 4),
                        "median_obs_pearson_full": round(float(np.nanmedian(acc["obsP"])), 4),
                        "median_obs_pearson_matched": round(float(np.nanmedian(acc["obsM"])), 4),
                        "frac_i_pearson_matched_n": round(float(np.mean(acc["i"])), 4),
                        "frac_ii_pearson_half_vs_full": round(float(np.mean(acc["ii"])), 4),
                        "frac_iii_cosine_half_vs_full": round(float(np.mean(acc["iii"])), 4),
                        "frac_iv_pearson_SB_vs_full": round(float(np.mean(acc["iv"])), 4),
                        "frac_v_cosine_SB_vs_full": round(float(np.mean(acc["v"])), 4),
                    })
            print(f"  {dname} shift={shift} 完了")

    T = pd.DataFrame(rows)
    T["contrib_sample_size"] = (T.frac_ii_pearson_half_vs_full
                                - T.frac_i_pearson_matched_n).round(4)
    T["contrib_uncentred"] = (T.frac_iii_cosine_half_vs_full
                              - T.frac_ii_pearson_half_vs_full).round(4)
    T["contrib_uncentred_SB"] = (T.frac_v_cosine_SB_vs_full
                                 - T.frac_iv_pearson_SB_vs_full).round(4)
    T["SB_minus_raw_cosine"] = (T.frac_v_cosine_SB_vs_full
                                - T.frac_iii_cosine_half_vs_full).round(4)
    T.to_csv(OUT / "ceiling_decomposition_sim.tsv", sep="\t", index=False)

    # 実データの信頼性帯（全 state を包含）かつ真の cos < 1 の範囲でまとめる
    rel = pd.read_csv(HERE / "results" / "reliability" / "reliability.tsv", sep="\t")
    lo, hi = float(rel.r_half_median.min()), float(rel.r_half_median.max())
    print(f"\n実データの r_half: {lo:.3f} - {hi:.3f}")
    band = T[(T.mean_r_half_cosine.between(lo, hi)) & (T.true_cos < 1.0)]
    print(f"この帯にある条件 {len(band)} 件（真の cos < 1）:")
    for col, lab in [("frac_i_pearson_matched_n", "(i)   中心化 x 標本数そろえ  "),
                     ("frac_ii_pearson_half_vs_full", "(ii)  中心化 x 半標本基準  "),
                     ("frac_iii_cosine_half_vs_full", "(iii) 非中心化 x 半標本基準"),
                     ("frac_iv_pearson_SB_vs_full", "(iv)  中心化 x SB 補正基準  "),
                     ("frac_v_cosine_SB_vs_full", "(v)   非中心化 x SB 補正基準")]:
        print(f"  {lab}: 平均 {band[col].mean():.4f}  中央 {band[col].median():.4f}  "
              f"最大 {band[col].max():.4f}")
    print(f"  標本数の寄与 (ii)-(i):  平均 {band.contrib_sample_size.mean():+.4f}  "
          f"中央 {band.contrib_sample_size.median():+.4f}")
    print(f"  非中心化の寄与 (iii)-(ii): 平均 {band.contrib_uncentred.mean():+.4f}  "
          f"中央 {band.contrib_uncentred.median():+.4f}")
    print(f"  SB 補正の効果 (v)-(iii):   平均 {band.SB_minus_raw_cosine.mean():+.4f}  "
          f"中央 {band.SB_minus_raw_cosine.median():+.4f}")
    print(f"  SB 基準での非中心化の寄与 (v)-(iv): 平均 {band.contrib_uncentred_SB.mean():+.4f}")

    print("\n真の cos = 1（同一応答）のとき、基準値を超える割合:")
    one = T[(T.true_cos == 1.0) & (T.mean_r_half_cosine.between(lo, hi))]
    for col, lab in [("frac_i_pearson_matched_n", "(i)  "),
                     ("frac_ii_pearson_half_vs_full", "(ii) "),
                     ("frac_iii_cosine_half_vs_full", "(iii)"),
                     ("frac_iv_pearson_SB_vs_full", "(iv) "),
                     ("frac_v_cosine_SB_vs_full", "(v)  ")]:
        print(f"  {lab} 平均 {one[col].mean():.4f}")

    print("\n真の cos と観測 cos の対応（mouse-like、shift 中央値）:")
    v = T[(T.design == "mouse-like (3 vs 6)") & (T["shift"] == shifts[1])]
    print(v.pivot_table(index="true_cos", columns="mean_r_half_cosine",
                        values="median_obs_cos_full").round(3).to_string())
    print(f"\n書き出し: {OUT / 'ceiling_decomposition_sim.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
