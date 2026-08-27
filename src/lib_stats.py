"""共通ユーティリティ：設定読み込み、BH補正、moderated t 統計量（Smyth 2004）。

n=3 の素の t 検定は分散推定が不安定なため、limma の eBayes に相当する
経験ベイズ分散縮小を実装して用いる。R に依存せず Python のみで完結させる。
"""
from __future__ import annotations
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats, special, optimize

ROOT = Path(__file__).resolve().parents[1]


def load_config(path: str | Path = None) -> dict:
    path = Path(path) if path else ROOT / "config" / "config.yaml"
    with open(path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_root"] = ROOT
    return cfg


def get_logger(name: str) -> logging.Logger:
    (ROOT / "logs").mkdir(exist_ok=True)
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    fh = logging.FileHandler(ROOT / "logs" / f"{name}.log", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg。NaN は NaN のまま返す。"""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


# --------------------------------------------------------------------------
# moderated t (empirical Bayes variance shrinkage)
# --------------------------------------------------------------------------
def _fit_fdist(s2: np.ndarray, df1: int) -> tuple[float, float]:
    """Smyth (2004) の scaled F 分布あてはめで事前分散 s0^2 と自由度 d0 を推定。

    log(s2) ~ log(s0^2) + log(F_{df1, d0}) をモーメント法で解く。
    """
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    if s2.size < 10:
        return float(np.median(s2)) if s2.size else 1.0, 0.0
    z = np.log(s2)
    ez = z.mean()
    vz = z.var(ddof=1)
    # E[log F] = digamma(df1/2) - digamma(d0/2) + log(d0/df1)
    # Var[log F] = trigamma(df1/2) + trigamma(d0/2)
    target = vz - special.polygamma(1, df1 / 2)
    if target <= 0:
        # 事前分散の自由度が無限大（全遺伝子で分散共通）に相当
        d0 = np.inf
        s0_2 = float(np.exp(ez + special.digamma(df1 / 2) - np.log(df1 / 2)))
        return s0_2, d0

    def f(x):
        return special.polygamma(1, x / 2) - target

    try:
        d0 = optimize.brentq(f, 1e-6, 1e6)
    except ValueError:
        d0 = 4.0
    s0_2 = float(np.exp(ez + special.digamma(df1 / 2) - np.log(df1 / 2)
                        - special.digamma(d0 / 2) + np.log(d0 / 2)))
    return s0_2, float(d0)


def moderated_ttest(A: pd.DataFrame, B: pd.DataFrame) -> pd.DataFrame:
    """2群比較の moderated t。A, B は行=feature, 列=sample の log スケール行列。

    返り値: lfc, t, p, q, df_total, s2, n_a, n_b
    """
    a = A.to_numpy(dtype=float)
    b = B.to_numpy(dtype=float)
    na = np.sum(~np.isnan(a), axis=1)
    nb = np.sum(~np.isnan(b), axis=1)
    ma = np.nanmean(a, axis=1)
    mb = np.nanmean(b, axis=1)
    va = np.nanvar(a, axis=1, ddof=1)
    vb = np.nanvar(b, axis=1, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        df_res = na + nb - 2
        s2 = ((na - 1) * va + (nb - 1) * vb) / df_res
        lfc = ma - mb

    usable = np.isfinite(s2) & (df_res > 0)
    df1 = int(np.nanmedian(df_res[usable])) if usable.any() else 1
    s0_2, d0 = _fit_fdist(s2[usable], df1)

    if np.isinf(d0):
        s2_post = np.full_like(s2, s0_2)
        df_total = np.full_like(df_res, np.inf, dtype=float)
    else:
        s2_post = (d0 * s0_2 + df_res * s2) / (d0 + df_res)
        df_total = df_res + d0

    with np.errstate(invalid="ignore", divide="ignore"):
        se = np.sqrt(s2_post * (1.0 / na + 1.0 / nb))
        t = lfc / se
    p = 2 * stats.t.sf(np.abs(t), df_total)

    out = pd.DataFrame(
        {
            "lfc": lfc,
            "t": t,
            "p": p,
            "s2": s2,
            "s2_post": s2_post,
            "df_total": df_total,
            "n_a": na,
            "n_b": nb,
        },
        index=A.index,
    )
    out["q"] = bh(out["p"].to_numpy())
    out.attrs["prior_df"] = d0
    out.attrs["prior_var"] = s0_2
    return out


def hedges_g(A: pd.DataFrame, B: pd.DataFrame) -> pd.Series:
    a, b = A.to_numpy(float), B.to_numpy(float)
    na, nb = a.shape[1], b.shape[1]
    sp = np.sqrt(((na - 1) * np.nanvar(a, 1, ddof=1) + (nb - 1) * np.nanvar(b, 1, ddof=1)) / (na + nb - 2))
    d = (np.nanmean(a, 1) - np.nanmean(b, 1)) / sp
    J = 1 - 3 / (4 * (na + nb) - 9)
    return pd.Series(d * J, index=A.index)


def leave_one_out_contrast(A: pd.DataFrame, B: pd.DataFrame) -> pd.DataFrame:
    """各 feature について、片群から1検体ずつ抜いたときの lfc 範囲と符号安定性。"""
    a, b = A.to_numpy(float), B.to_numpy(float)
    full = np.nanmean(a, 1) - np.nanmean(b, 1)
    vals = []
    for k in range(a.shape[1]):
        vals.append(np.nanmean(np.delete(a, k, axis=1), 1) - np.nanmean(b, 1))
    for k in range(b.shape[1]):
        vals.append(np.nanmean(a, 1) - np.nanmean(np.delete(b, k, axis=1), 1))
    V = np.vstack(vals)
    lo, hi = V.min(0), V.max(0)
    return pd.DataFrame(
        {"lfc": full, "loo_min": lo, "loo_max": hi,
         "sign_stable": np.sign(lo) == np.sign(hi)},
        index=A.index,
    )


def competitive_z(values: pd.Series, gene_set: list[str], min_n: int = 5):
    """順位ベースの競合的検定。集合の平均順位パーセンタイルを標準化して z を返す。"""
    v = values.dropna()
    v = v[np.isfinite(v)]
    g = [x for x in gene_set if x in v.index]
    if len(g) < min_n:
        return np.nan, np.nan, len(g)
    r = v.rank(pct=True)
    z = (r[g].mean() - 0.5) / np.sqrt(1 / 12 / len(g))
    return float(z), float(v[g].mean()), len(g)


def decile_matched_permutation(x: pd.Series, y: pd.Series, abundance: pd.Series,
                               n_perm: int = 1000, n_decile: int = 10, seed: int = 0):
    """発現量十分位内で y をシャッフルして Spearman rho の帰無分布を作る。

    単純な全体シャッフルは楽観的すぎるため、存在量依存の構造を保存する。
    """
    idx = x.index.intersection(y.index).intersection(abundance.index)
    x, y, ab = x.loc[idx], y.loc[idx], abundance.loc[idx]
    ok = x.notna() & y.notna() & ab.notna() & np.isfinite(x) & np.isfinite(y) & np.isfinite(ab)
    x, y, ab = x[ok], y[ok], ab[ok]
    obs = stats.spearmanr(x, y)[0]
    dec = pd.qcut(ab.rank(method="first"), n_decile, labels=False).to_numpy()
    rng = np.random.default_rng(seed)
    xv, yv = x.to_numpy(), y.to_numpy()
    null = np.empty(n_perm)
    groups = [np.where(dec == d)[0] for d in range(n_decile)]
    for i in range(n_perm):
        ys = yv.copy()
        for m in groups:
            ys[m] = rng.permutation(ys[m])
        null[i] = stats.spearmanr(xv, ys)[0]
    p = (np.sum(null >= obs) + 1) / (n_perm + 1)
    return {"rho": float(obs), "null_mean": float(null.mean()),
            "null_sd": float(null.std()), "z": float((obs - null.mean()) / null.std()),
            "p_perm": float(p), "n": int(len(x))}


def append_summary(section: str, payload: dict, cfg: dict):
    """results/summary.md に主要数値を追記（人が読む用）＋ JSON 保存（機械可読）。"""
    res = Path(cfg["_root"]) / cfg["paths"]["results"]
    res.mkdir(parents=True, exist_ok=True)
    with open(res / "summary.md", "a", encoding="utf-8") as fh:
        fh.write(f"\n## {section}\n\n")
        for k, v in payload.items():
            fh.write(f"- **{k}**: {v}\n")
    jpath = res / "summary.json"
    blob = json.loads(jpath.read_text()) if jpath.exists() else {}
    blob[section] = payload
    jpath.write_text(json.dumps(blob, ensure_ascii=False, indent=2, default=str))


def alias_lookup(target_index, aliases: str, primary: str = None):
    """データソース間でシンボルの世代がずれていても引けるようにする。

    同じ遺伝子でもファイルごとに採用シンボルが違う（FAM3D の例:
    RNA-seq FPKM は旧称 Oit1、プロテオームと UniProt は現行名 Fam3d）。
    ortholog_map.tsv の mouse_aliases 列（'|' 区切り）を順に試し、
    target_index に実在する最初の名前を返す。無ければ None。

    大文字小文字の違いも吸収する（UniProt 由来の表は大文字化されている）。
    """
    idx = {str(i).upper(): i for i in target_index}
    names = []
    if isinstance(primary, str) and primary:
        names.append(primary)
    if isinstance(aliases, str) and aliases:
        names.extend(aliases.split("|"))
    for n in names:
        hit = idx.get(str(n).upper())
        if hit is not None:
            return hit
    return None
