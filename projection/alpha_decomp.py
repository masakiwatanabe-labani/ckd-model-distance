"""Step 2: 射影指標を「方向 cos」と「応答振幅 nrm」に分解する。

基準軸 a（既定ではネコ CKD3/4 の応答ベクトル）に対し、各状態のベクトル v を

    alpha  = (a·v)/(a·a)          = (|v|/|a|) * cos(theta)   = nrm * cos
    cos    = (a·v)/(|a||v|)
    nrm    = |v|/|a|
    R_perp = |v - alpha*a| / |v|  = sin(theta)                = sqrt(1 - cos^2)

alpha は角度と応答振幅の積であり、R_perp は sin(theta) そのものなので、
alpha と R_perp は独立な2軸ではない。alpha 単独では

    「向きが基準に近い」  と  「応答が大きい」

を区別できず、後者だけで alpha が動く。したがって alpha を「再現度」として
読むことはできない。本スクリプトは alpha を必ず cos と nrm に分けて出す。

時間との関係も同じ理由で3本立てで見る:
    nrm x log_time          振幅が時間で説明されるか（旧稿 Mantel rho の実体）
    cos x log_time          方向が時間で説明されるか
    cos | nrm (partial)     振幅を統制してなお方向が時間と関係するか  <- 結論はここ

状態数は一桁〜十数個しかない。Spearman の CI はきわめて広いので、
点推定の大小ではなく CI を見ること。本スクリプトは有意性の判定を出さず、
CI と、旧稿の Mantel rho = 0.553 を CI が含むかどうかだけを出す。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PRIOR_MANTEL_RHO = 0.553      # 旧稿の Mantel rho（距離 x 経過時間）
TOL = 1e-9


# ------------------------------------------------------------------ 幾何
def decompose(a: np.ndarray, V: np.ndarray) -> pd.DataFrame:
    """基準ベクトル a と状態行列 V（遺伝子 x 状態）から alpha / cos / nrm / R_perp。"""
    aa = float(a @ a)
    na = float(np.sqrt(aa))
    dot = a @ V                                   # (states,)
    nv = np.sqrt((V ** 2).sum(axis=0))            # (states,)
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha = dot / aa
        cos = dot / (na * nv)
        nrm = nv / na
        resid = V - np.outer(a, alpha)
        R_perp = np.sqrt((resid ** 2).sum(axis=0)) / nv
    return pd.DataFrame({"alpha": alpha, "cos": cos, "nrm": nrm, "R_perp": R_perp})


def check_identities(df: pd.DataFrame) -> None:
    """alpha == nrm*cos と R_perp == sin(theta) を検証する。

    ここで落ちたらデータ側の問題（NaN の混入、状態間で遺伝子集合が
    揃っていない等）。assert を外して先に進んではいけない。
    """
    ok = df[["alpha", "cos", "nrm", "R_perp"]].notna().all(axis=1)
    d = df[ok]
    bad = np.abs(d["alpha"] - d["nrm"] * d["cos"]) > TOL * np.maximum(1.0, np.abs(d["alpha"]))
    assert not bad.any(), (
        "alpha != nrm * cos が成立しない状態がある: "
        f"{list(d.index[bad])}\n{d[bad].to_string()}")
    sin_from_cos = np.sqrt(np.clip(1.0 - d["cos"] ** 2, 0.0, None))
    bad2 = np.abs(d["R_perp"] - sin_from_cos) > 1e-7
    assert not bad2.any(), (
        "R_perp != sqrt(1 - cos^2) が成立しない状態がある: "
        f"{list(d.index[bad2])}\n{d[bad2].to_string()}")


# ------------------------------------------------------------------ 相関
def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return np.nan
    return float(stats.spearmanr(x[m], y[m])[0])


def partial_spearman(x, y, z):
    """z を統制した x, y の偏 Spearman（順位に直してから偏 Pearson）。"""
    x, y, z = (np.asarray(v, float) for v in (x, y, z))
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if m.sum() < 4:
        return np.nan
    rx, ry, rz = (stats.rankdata(v[m]) for v in (x, y, z))
    if min(np.ptp(rx), np.ptp(ry), np.ptp(rz)) == 0:
        return np.nan
    rxy, rxz, ryz = (np.corrcoef(p, q)[0, 1] for p, q in ((rx, ry), (rx, rz), (ry, rz)))
    den = np.sqrt(max((1 - rxz ** 2) * (1 - ryz ** 2), 0.0))
    return float((rxy - rxz * ryz) / den) if den > 1e-12 else np.nan


def fisher_ci(r, n, k=0, level=0.95):
    """Fisher-z 近似 CI。k は統制した変数の数。"""
    if not np.isfinite(r) or n - 3 - k <= 0:
        return (np.nan, np.nan)
    z, se = np.arctanh(np.clip(r, -0.999999, 0.999999)), 1 / np.sqrt(n - 3 - k)
    q = stats.norm.ppf(0.5 + level / 2)
    return float(np.tanh(z - q * se)), float(np.tanh(z + q * se))


def boot_ci(fn, n, rng, n_boot, level=0.95):
    """状態をブートストラップした percentile CI。"""
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(idx).size < 4:
            continue
        v = fn(idx)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < 50:
        return (np.nan, np.nan, len(vals))
    lo, hi = np.percentile(vals, [(1 - level) / 2 * 100, (0.5 + level / 2) * 100])
    return float(lo), float(hi), len(vals)


def perm_p(fn_stat, obs, logt, rng, n_perm):
    """時間ラベルの並べ替えによる両側 p。"""
    if not np.isfinite(obs):
        return np.nan
    cnt = 0
    for _ in range(n_perm):
        v = fn_stat(rng.permutation(logt))
        if np.isfinite(v) and abs(v) >= abs(obs) - 1e-12:
            cnt += 1
    return (cnt + 1) / (n_perm + 1)


def time_association(P: pd.DataFrame, hours: pd.Series, rng, n_boot, n_perm,
                     label: str) -> pd.DataFrame:
    """P（状態 x 指標）と経過時間の関係。log10(hours) を時間軸に使う。"""
    use = P.index[hours.reindex(P.index).notna()]
    if len(use) < 4:
        return pd.DataFrame()
    logt = np.log10(hours.reindex(use).to_numpy(float))
    A = P.loc[use, "alpha"].to_numpy(float)
    C = P.loc[use, "cos"].to_numpy(float)
    N = P.loc[use, "nrm"].to_numpy(float)
    R = P.loc[use, "R_perp"].to_numpy(float)
    n = len(use)

    terms = [
        ("alpha x log_time", lambda i: spearman(A[i], logt[i]),
         lambda t: spearman(A, t), spearman(A, logt), 0),
        ("nrm x log_time", lambda i: spearman(N[i], logt[i]),
         lambda t: spearman(N, t), spearman(N, logt), 0),
        ("cos x log_time", lambda i: spearman(C[i], logt[i]),
         lambda t: spearman(C, t), spearman(C, logt), 0),
        ("R_perp x log_time", lambda i: spearman(R[i], logt[i]),
         lambda t: spearman(R, t), spearman(R, logt), 0),
        ("cos | nrm (partial)", lambda i: partial_spearman(C[i], logt[i], N[i]),
         lambda t: partial_spearman(C, t, N), partial_spearman(C, logt, N), 1),
        ("nrm | cos (partial)", lambda i: partial_spearman(N[i], logt[i], C[i]),
         lambda t: partial_spearman(N, t, C), partial_spearman(N, logt, C), 1),
        ("alpha | nrm (partial)", lambda i: partial_spearman(A[i], logt[i], N[i]),
         lambda t: partial_spearman(A, t, N), partial_spearman(A, logt, N), 1),
    ]
    ctrl_col = {"cos | nrm (partial)": spearman(N, logt),
                "nrm | cos (partial)": spearman(C, logt),
                "alpha | nrm (partial)": spearman(N, logt)}
    rows = []
    for name, bfn, pfn, obs, k in terms:
        flo, fhi = fisher_ci(obs, n, k)
        blo, bhi, nb = boot_ci(bfn, n, rng, n_boot)
        p = perm_p(pfn, obs, logt, rng, n_perm)
        rows.append({
            "subset": label, "term": name, "n_states": n, "rho": obs,
            "ctrl_x_log_time": ctrl_col.get(name, np.nan),
            "fisher_lo": flo, "fisher_hi": fhi,
            "boot_lo": blo, "boot_hi": bhi, "n_boot_ok": nb,
            "p_perm": p,
            "fisher_ci_contains_0.553": (bool(flo <= PRIOR_MANTEL_RHO <= fhi)
                                         if np.isfinite(flo) and np.isfinite(fhi) else None),
            "boot_ci_contains_0.553": (bool(blo <= PRIOR_MANTEL_RHO <= bhi)
                                       if np.isfinite(blo) and np.isfinite(bhi) else None),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ 本解析
def run_real(a_args):
    D = pd.read_csv(a_args.delta, sep="\t", index_col=0)
    if a_args.ref not in D.columns:
        raise SystemExit(f"基準状態 {a_args.ref} が {a_args.delta} にありません。列: {list(D.columns)}")
    if a_args.genes:
        want = {g.strip() for g in Path(a_args.genes).read_text().split() if g.strip()}
        D = D.loc[D.index.intersection(sorted(want))]
    n_before = len(D)
    # 全状態で有限な遺伝子のみ。|v| を状態間で比較する以上、状態ごとに
    # 遺伝子集合が違ってはいけない（欠測数の差が振幅差に化ける）。
    D = D.dropna(axis=0, how="any")
    n_genes = len(D)
    if n_genes < 50:
        raise SystemExit(f"complete case の遺伝子が {n_genes} 個しかありません")

    states = [c for c in D.columns]
    a = D[a_args.ref].to_numpy(float)
    V = D[states].to_numpy(float)
    P = decompose(a, V)
    P.index = states
    check_identities(P)

    rng = np.random.default_rng(a_args.seed)
    # 遺伝子ブートストラップ（状態ごとの alpha/cos/nrm の不確かさ）
    acc = {k: np.empty((a_args.n_boot_gene, len(states))) for k in ("alpha", "cos", "nrm")}
    for b in range(a_args.n_boot_gene):
        idx = rng.integers(0, n_genes, n_genes)
        d = decompose(a[idx], V[idx])
        for k in acc:
            acc[k][b] = d[k].to_numpy()
    for k, M in acc.items():
        P[f"{k}_lo95"] = np.percentile(M, 2.5, axis=0)
        P[f"{k}_hi95"] = np.percentile(M, 97.5, axis=0)

    tm = pd.read_csv(a_args.time, sep="\t")
    hours = tm.set_index("state")["hours"].astype(float)
    P.insert(0, "hours", hours.reindex(P.index))
    P.insert(0, "n_genes", n_genes)

    ta = [time_association(P, hours, rng, a_args.n_boot, a_args.n_perm, "all_model_states")]
    for pref, lab in [("IRI_", "IRI_only"), ("mouse_", "PodTRECK_only")]:
        sub = P.loc[[s for s in P.index if s.startswith(pref)]]
        if len(sub) >= 4:
            ta.append(time_association(sub, hours, rng, a_args.n_boot, a_args.n_perm, lab))
    TA = pd.concat([t for t in ta if len(t)], ignore_index=True)

    out = Path(a_args.out); out.mkdir(parents=True, exist_ok=True)
    P.round(6).to_csv(out / "projection.tsv", sep="\t")
    TA.round(6).to_csv(out / "time_association.tsv", sep="\t", index=False)
    write_report(out, a_args, P, TA, n_genes, n_before)
    print(P[["n_genes", "hours", "alpha", "cos", "nrm", "R_perp"]].round(4).to_string())
    print()
    print(TA[["subset", "term", "n_states", "rho", "boot_lo", "boot_hi",
              "p_perm", "boot_ci_contains_0.553"]].round(4).to_string(index=False))


def write_report(out, a_args, P, TA, n_genes, n_before):
    L = [f"# alpha 分解: ref = {a_args.ref}", "",
         f"- 遺伝子リスト: `{a_args.genes or '(なし: Δ行列の全遺伝子)'}`",
         f"- リスト適用後 {n_before} 遺伝子 -> 全状態 complete case {n_genes} 遺伝子",
         f"- 状態: {len(P)}", "",
         "## 状態ごとの分解", "",
         "alpha は必ず cos と nrm と併記する。alpha = nrm x cos なので、",
         "alpha の動きが方向由来か振幅由来かは alpha だけでは判別できない。", "",
         "| 状態 | hours | alpha [95% CI] | cos [95% CI] | nrm [95% CI] | R_perp |",
         "|---|---|---|---|---|---|"]
    for s, r in P.iterrows():
        h = "-" if not np.isfinite(r["hours"]) else f"{r['hours']:g}"
        L.append(f"| {s} | {h} | {r['alpha']:.3f} [{r['alpha_lo95']:.3f}, {r['alpha_hi95']:.3f}] "
                 f"| {r['cos']:.3f} [{r['cos_lo95']:.3f}, {r['cos_hi95']:.3f}] "
                 f"| {r['nrm']:.3f} [{r['nrm_lo95']:.3f}, {r['nrm_hi95']:.3f}] "
                 f"| {r['R_perp']:.3f} |")
    L += ["", "（CI は遺伝子のブートストラップ。状態間の比較可能性は担保しない）", "",
          "## 経過時間との関係", "",
          "状態数が一桁〜十数個しかないため、CI はきわめて広い。点推定の大小ではなく",
          f"CI を読むこと。参照値は旧稿の Mantel rho = {PRIOR_MANTEL_RHO}。", "",
          "| 部分集合 | 項 | 状態数 | rho | bootstrap CI | Fisher CI | p(perm) | CIが0.553を含む |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in TA.iterrows():
        rho_s = "未定義" if not np.isfinite(r["rho"]) else f"{r['rho']:+.3f}"
        L.append(f"| {r['subset']} | {r['term']} | {r['n_states']} | {rho_s} "
                 f"| [{r['boot_lo']:+.3f}, {r['boot_hi']:+.3f}] "
                 f"| [{r['fisher_lo']:+.3f}, {r['fisher_hi']:+.3f}] "
                 f"| {r['p_perm']:.4f} | {r['boot_ci_contains_0.553']} |")
    L += ["", "## 読み方", "",
          "1. `nrm x log_time` が高ければ、旧稿の Mantel rho は振幅の時間依存として再現される。",
          "2. `cos x log_time` が低ければ、方向成分は時間では説明されない。",
          "3. `cos | nrm (partial)` が実質ゼロなら、振幅を統制した時点で方向と時間の関係は残らない。",
          "   「alpha は時間の言い換えではない」と主張できるのはこの偏相関だけ。", "",
          "偏相関が「未定義」と出る場合は、統制変数（nrm など）と log_time の順位相関が",
          "±1 に達していて、時間の分散が統制変数で完全に説明され尽くしていることを意味する。",
          "これは『効果がゼロ』ではなく『この状態集合では分離できない』という結果であり、",
          "`ctrl_x_log_time` 列で確認できる。", "",
          "`ctrl_x_log_time` の絶対値が 1 に近いとき、偏相関の分母 sqrt(1 - r^2) が",
          "小さくなり、わずかな残差相関でも点推定が大きく振れる。この場合は点推定ではなく",
          "bootstrap CI を読むこと（合成データでの検証でもこの増幅を確認している）。", "",
          "## 制約", "",
          "- Δ の単位はデータセット間で完全には揃わない（ネコ RNA はシートの log2 正規化値、",
          "  Pod-TRECK は log2(FPKM+1)、IRI は log2(正規化カウント+1)）。nrm = |v|/|a| は",
          "  生物学的な応答振幅とプラットフォームのダイナミックレンジの両方を含む。",
          "  同一プラットフォーム内の比較（ネコ CKD1/2 vs CKD3/4、IRI の時点間、",
          "  Pod-TRECK の時点間）は解釈できるが、種をまたぐ nrm の比較には",
          "  この交絡が乗る。",
          "- 状態数が少ないため、相関の CI は ±0.5 程度の幅を持つ。CI の重なりから",
          "  「差がある」とも「差がない」とも結論できない。CI をそのまま示すこと。"]
    (out / "report.md").write_text("\n".join(L) + "\n")


# ------------------------------------------------------------------ 自己検証
def run_demo(a_args):
    """既知の答えを持つ合成データで、分解と時間相関の挙動を検証する。

    S1「振幅だけが時間依存」: cos は時間と無関係、nrm は時間に単調増加。
        -> alpha は時間と強く相関するが、それは振幅由来。cos|nrm は実質ゼロ。
    S2「方向だけが時間依存」: cos が時間に単調増加、nrm は時間と無関係。
        -> cos|nrm は強く残る。

    S2 は「偏相関がゼロに出るのは手法が鈍いからではない」ことの対照。
    S1 だけでは、偏相関が何も検出できない手法である可能性を排除できない。
    """
    rng = np.random.default_rng(a_args.seed)
    n_genes, n_states = 800, 9
    hours = np.array([2, 6, 24, 72, 168, 336, 672, 2016, 8760], float)
    a = rng.normal(size=n_genes)
    na = np.linalg.norm(a)
    ahat = a / na

    def build(cos_vals, nrm_vals):
        V = np.empty((n_genes, n_states))
        for i in range(n_states):
            b = rng.normal(size=n_genes)
            b -= (b @ ahat) * ahat
            b /= np.linalg.norm(b)
            sin = np.sqrt(max(1 - cos_vals[i] ** 2, 0.0))
            V[:, i] = nrm_vals[i] * na * (cos_vals[i] * ahat + sin * b)
        return V

    # S1 の cos は「時間とも nrm とも順位が無関係」でなければならない。
    # 時間とだけ無相関にしても、nrm と相関していると偏相関に効いてしまう。
    # さらに統制変数 nrm が時間とほぼ共線なので偏相関の分母が小さくなり、
    # わずかな残差相関が増幅される。決定的な探索で両方ゼロ近傍の並べ替えを選ぶ。
    _logt = np.log10(hours)
    _nrm_design = np.array([0.30, 0.55, 0.70, 1.00, 0.90, 1.30, 1.20, 1.70, 1.60])
    _base = np.array([0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50])
    _r = np.random.default_rng(0)
    shuffled, _best = None, np.inf
    for _ in range(20000):
        cand = _r.permutation(_base)
        sc = max(abs(spearman(cand, _logt)), abs(spearman(cand, _nrm_design)))
        if sc < _best:
            shuffled, _best = cand, sc
    monotone_cos = np.linspace(0.20, 0.90, n_states)
    # 単調増加そのものにすると nrm と時間の順位相関が厳密に 1 になり、
    # 時間の分散が nrm で完全に説明されて偏相関が「ゼロ」ではなく「未定義」になる。
    # 実データの Pod-TRECK（Day14 ピーク・Day21 低下）に倣い、増加基調のまま
    # 順位を少し崩しておく。
    monotone_nrm = _nrm_design
    flat_nrm = np.array([1.05, 0.95, 1.10, 0.98, 1.02, 0.92, 1.08, 1.00, 0.90])

    scen = {
        "S1_amplitude_only": (shuffled, monotone_nrm),
        "S2_direction_only": (monotone_cos, flat_nrm),
    }
    out = Path(a_args.out); out.mkdir(parents=True, exist_ok=True)
    fails, all_P, all_TA = [], [], []

    for name, (cv, nv) in scen.items():
        V = build(cv, nv)
        P = decompose(a, V)
        P.index = [f"S{i}" for i in range(n_states)]
        check_identities(P)                      # 恒等式（ここが本体の assert）
        hs = pd.Series(hours, index=P.index)
        P.insert(0, "hours", hs)

        def chk(cond, msg):
            if not cond:
                fails.append(f"[{name}] {msg}")

        chk(np.allclose(P["cos"], cv, atol=1e-9), "cos が設計値と一致しない")
        chk(np.allclose(P["nrm"], nv, atol=1e-9), "nrm が設計値と一致しない")
        chk(np.allclose(P["alpha"], np.asarray(nv) * np.asarray(cv), atol=1e-9),
            "alpha != nrm*cos")
        chk(np.allclose(P["R_perp"], np.sqrt(1 - np.asarray(cv) ** 2), atol=1e-7),
            "R_perp != sin(theta)")

        TA = time_association(P, hs, rng, a_args.n_boot, a_args.n_perm, name)
        g = TA.set_index("term")["rho"]
        if name == "S1_amplitude_only":
            chk(g["nrm x log_time"] > 0.85, f"nrm x time が低い ({g['nrm x log_time']:.3f})")
            chk(g["alpha x log_time"] > 0.70,
                f"alpha x time が低い ({g['alpha x log_time']:.3f})")
            chk(np.isfinite(g["cos | nrm (partial)"]),
                "cos|nrm が未定義（統制変数が時間と完全共線）")
            chk(abs(g["cos x log_time"]) < 0.25,
                f"cos x time が大きい ({g['cos x log_time']:.3f})")
            chk(abs(g["cos | nrm (partial)"]) < 0.35,
                f"cos|nrm がゼロ近傍でない ({g['cos | nrm (partial)']:.3f})")
        else:
            chk(g["cos x log_time"] > 0.99, f"cos x time が低い ({g['cos x log_time']:.3f})")
            chk(g["cos | nrm (partial)"] > 0.80,
                f"cos|nrm が検出できていない ({g['cos | nrm (partial)']:.3f})")
        all_P.append(P.assign(scenario=name))
        all_TA.append(TA)

    pd.concat(all_P).round(6).to_csv(out / "demo_projection.tsv", sep="\t")
    pd.concat(all_TA).round(6).to_csv(out / "demo_time_association.tsv", sep="\t", index=False)
    print(pd.concat(all_P)[["scenario", "hours", "alpha", "cos", "nrm", "R_perp"]]
          .round(4).to_string())
    print()
    print(pd.concat(all_TA)[["subset", "term", "rho", "boot_lo", "boot_hi", "p_perm"]]
          .round(4).to_string(index=False))
    print()
    if fails:
        print("DEMO FAILED:")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("DEMO PASSED: 恒等式（alpha = nrm*cos, R_perp = sin）と、"
          "振幅由来／方向由来の判別が両方とも期待どおり。")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delta")
    ap.add_argument("--ref")
    ap.add_argument("--genes")
    ap.add_argument("--time")
    ap.add_argument("--out", required=True)
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--n-boot", type=int, default=5000, help="状態ブートストラップ回数")
    ap.add_argument("--n-boot-gene", type=int, default=2000, help="遺伝子ブートストラップ回数")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260826)
    a = ap.parse_args()
    if a.demo:
        run_demo(a)
    else:
        missing = [k for k in ("delta", "ref", "time") if not getattr(a, k)]
        if missing:
            raise SystemExit(f"--demo でない場合は必須: {missing}")
        run_real(a)


if __name__ == "__main__":
    main()
