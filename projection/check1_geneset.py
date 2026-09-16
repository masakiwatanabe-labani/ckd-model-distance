"""CHECK 1: 旧稿との rho 不一致が「遺伝子集合の入れ替わり」で説明できるかを判定する。

問い:
    Ensembl オルソログ表で引けた遺伝子（core）に対し、大文字シンボル一致だけで
    拾われた遺伝子（added）を足すと種間 rho が動く。旧稿 0.490 と新稿の値の差は
    「added がノイズ寄りだから」で説明できるのか。

この説明を post hoc で終わらせないための手続き:
    added が core より低発現に偏っていること自体は、rho が低い理由の説明にならない。
    発現量分布を揃えた core の部分集合（added と同じ十分位構成でリサンプリング）
    を作り、その rho の分布に対して added の rho がどこに落ちるかで判定する。
    発現量マッチ後も added が低いなら「added は本当に質が低い」。
    マッチすると差が消えるなら「発現量の違いを言い換えていただけ」。
    added が低くないなら、そもそもこの説明では不一致を説明できない。

注意:
    ここでの n は遺伝子数（数千）なので、Spearman の CI は十分に狭く p 値も
    意味を持つ。状態数が一桁である Step 2/3 の時間相関とは性質が違う。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

N_DECILE = 10


def fisher_ci(rho: float, n: int, level: float = 0.95) -> tuple[float, float]:
    """Spearman rho の Fisher-z 近似 CI（Bonett-Wright の分散補正つき）。"""
    if not np.isfinite(rho) or n < 6:
        return (np.nan, np.nan)
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
    q = stats.norm.ppf(0.5 + level / 2)
    return float(np.tanh(z - q * se)), float(np.tanh(z + q * se))


def boot_ci(x: np.ndarray, y: np.ndarray, n_boot: int, rng, level: float = 0.95):
    """遺伝子をブートストラップした Spearman rho の percentile CI。"""
    n = len(x)
    if n < 6:
        return (np.nan, np.nan)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        k = rng.integers(0, n, n)
        vals[i] = stats.spearmanr(x[k], y[k])[0]
    lo, hi = np.nanpercentile(vals, [(1 - level) / 2 * 100, (0.5 + level / 2) * 100])
    return float(lo), float(hi)


def nn_match(rank_a, rank_b, caliper, rng):
    """発現量順位での 1:1 最近傍マッチング（復元なし・caliper 付き）。

    十分位リサンプリングは、片方の集合がある十分位にほとんど存在しないと
    復元抽出に落ちて推定が不安定になる（Group A は低発現側がほぼ空）。
    共通サポートの範囲だけで 1:1 対応を作り、そこで比較する。

    返り値: (a 側の採用 index, b 側の採用 index)
    """
    order = rng.permutation(len(rank_a))
    bo = np.argsort(rank_b, kind="stable")
    rb = rank_b[bo]
    used = np.zeros(len(rb), bool)
    pa, pb = [], []
    for ai in order:
        t = rank_a[ai]
        pos = int(np.searchsorted(rb, t))
        best, bestd = -1, caliper
        l, r = pos - 1, pos
        while True:
            moved = False
            if r < len(rb) and (rb[r] - t) <= bestd:
                moved = True
                if not used[r]:
                    best, bestd = r, rb[r] - t
                r += 1
            if l >= 0 and (t - rb[l]) <= bestd:
                moved = True
                if not used[l] and (t - rb[l]) < bestd:
                    best, bestd = l, t - rb[l]
                l -= 1
            if not moved:
                break
        if best >= 0:
            used[best] = True
            pa.append(ai)
            pb.append(int(bo[best]))
    return np.array(pa, int), np.array(pb, int)


def joint_expression_rank(expr: pd.DataFrame, genes: pd.Index) -> pd.Series:
    """データセットごとにパーセンタイル順位へ直し、平均して1本の発現量共変量にする。

    列は単位が異なる（ネコのシート値 / log2(FPKM+1) / log2(count+1)）ので
    絶対値は比較できない。順位に直してから平均する。
    """
    e = expr.reindex(genes)
    ranks = e.rank(pct=True)
    return ranks.mean(axis=1, skipna=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delta", required=True)
    ap.add_argument("--x", required=True, help="基準側の状態（例 cat_CKD34）")
    ap.add_argument("--y", required=True, help="比較側の状態（例 mouse_2W）")
    ap.add_argument("--core", required=True, help="Ensembl ortholog 由来遺伝子のリスト")
    ap.add_argument("--expr", required=True, help="対照群の発現量テーブル")
    ap.add_argument("--out", required=True)
    ap.add_argument("--core-label", default="core", help="core 集合の表示名（例 GroupA）")
    ap.add_argument("--added-label", default="added", help="補集合の表示名（例 GroupB）")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-resample", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--caliper", type=float, default=0.02,
                    help="1:1 マッチングの許容差（発現量パーセンタイル順位）")
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    D = pd.read_csv(a.delta, sep="\t", index_col=0)
    for s in (a.x, a.y):
        if s not in D.columns:
            raise SystemExit(f"状態 {s} が {a.delta} にありません。列: {list(D.columns)}")
    core_list = [g.strip() for g in Path(a.core).read_text().split() if g.strip()]
    expr = pd.read_csv(a.expr, sep="\t", index_col=0)

    xy = D[[a.x, a.y]].dropna()
    genes = xy.index
    is_core = genes.isin(set(core_list))
    CL, AL = a.core_label, a.added_label
    grp = pd.Series(np.where(is_core, CL, AL), index=genes, name="group")

    # -------------------------------------------------- 1. 素の rho
    rows = []
    for lab, m in [("all", np.ones(len(genes), bool)),
                   (CL, is_core), (AL, ~is_core)]:
        xv, yv = xy[a.x].to_numpy()[m], xy[a.y].to_numpy()[m]
        rho = stats.spearmanr(xv, yv)[0] if m.sum() > 5 else np.nan
        flo, fhi = fisher_ci(rho, int(m.sum()))
        blo, bhi = boot_ci(xv, yv, a.n_boot, rng)
        rows.append({"set": lab, "n": int(m.sum()), "rho": rho,
                     "fisher_lo": flo, "fisher_hi": fhi,
                     "boot_lo": blo, "boot_hi": bhi})
    summary = pd.DataFrame(rows)

    # -------------------------------------------------- 1b. 素の差の遺伝子ブートストラップ CI
    # マッチしない差にも区間を付ける。core/added それぞれの遺伝子を各リサンプルで
    # 引き直し、同じリサンプル内での rho の差を取る。
    xa, ya = xy[a.x].to_numpy()[is_core], xy[a.y].to_numpy()[is_core]
    xb, yb = xy[a.x].to_numpy()[~is_core], xy[a.y].to_numpy()[~is_core]
    # 既存の乱数列を消費しないよう専用の RNG を使う。ここで rng を進めると
    # 後段の 1:1 マッチングと十分位リサンプリングの結果が変わってしまう。
    rng_u = np.random.default_rng(a.seed + 1)
    diffs = np.empty(a.n_boot)
    for i in range(a.n_boot):
        ka = rng_u.integers(0, len(xa), len(xa))
        kb = rng_u.integers(0, len(xb), len(xb))
        diffs[i] = (stats.spearmanr(xa[ka], ya[ka])[0]
                    - stats.spearmanr(xb[kb], yb[kb])[0])
    d_lo, d_hi = np.nanpercentile(diffs, [2.5, 97.5])
    unmatched = {f"delta_{CL}_minus_{AL}_unmatched": float(np.nanmedian(diffs)),
                 "point_estimate": float(stats.spearmanr(xa, ya)[0] - stats.spearmanr(xb, yb)[0]),
                 "lo95": float(d_lo), "hi95": float(d_hi), "n_boot": a.n_boot}
    pd.Series(unmatched).to_frame("value").to_csv(out / "check1_unmatched_delta.tsv", sep="\t")

    # -------------------------------------------------- 2. 発現量分布
    jr = joint_expression_rank(expr, genes)
    have_expr = jr.notna()
    summary["median_expr_rank"] = [
        float(jr[m & have_expr.to_numpy()].median()) if (m & have_expr.to_numpy()).any() else np.nan
        for m in [np.ones(len(genes), bool), is_core, ~is_core]]
    summary["n_with_expr"] = [int((m & have_expr.to_numpy()).sum())
                              for m in [np.ones(len(genes), bool), is_core, ~is_core]]

    # 発現量が取れた遺伝子だけで十分位を切る
    g2 = genes[have_expr.to_numpy()]
    xy2, jr2, core2 = xy.loc[g2], jr.loc[g2], is_core[have_expr.to_numpy()]
    dec = pd.qcut(jr2.rank(method="first"), N_DECILE, labels=False)

    by_dec = []
    for d in range(N_DECILE):
        md = (dec == d).to_numpy()
        row = {"decile": d + 1, "n_core": int((md & core2).sum()),
               "n_added": int((md & ~core2).sum())}
        for lab, mm in [(CL, md & core2), (AL, md & ~core2)]:
            if mm.sum() > 5:
                row[f"rho_{lab}"] = float(stats.spearmanr(
                    xy2[a.x].to_numpy()[mm], xy2[a.y].to_numpy()[mm])[0])
            else:
                row[f"rho_{lab}"] = np.nan
        row[f"delta_{AL}_minus_{CL}"] = row[f"rho_{AL}"] - row[f"rho_{CL}"]
        by_dec.append(row)
    by_dec = pd.DataFrame(by_dec)

    # -------------------------------------------------- 3. 発現量マッチ・リサンプリング
    # added と同じ十分位構成の core 部分集合を繰り返し引き、rho の帰無分布を作る。
    want = by_dec.set_index("decile")["n_added"].to_dict()
    pools = {d + 1: np.where((dec == d).to_numpy() & core2)[0] for d in range(N_DECILE)}
    xv2, yv2 = xy2[a.x].to_numpy(), xy2[a.y].to_numpy()

    short = {d: want[d] - len(pools[d]) for d in want if want[d] > len(pools[d])}
    null = np.empty(a.n_resample)
    for i in range(a.n_resample):
        pick = []
        for d, k in want.items():
            pool = pools[d]
            if k == 0 or pool.size == 0:
                continue
            # core 側が足りない十分位では復元抽出に落とす（該当は下の short に記録）
            pick.append(rng.choice(pool, size=min(k, pool.size), replace=False)
                        if k <= pool.size else rng.choice(pool, size=k, replace=True))
        idx = np.concatenate(pick)
        null[i] = stats.spearmanr(xv2[idx], yv2[idx])[0]

    rho_added = float(stats.spearmanr(xv2[~core2], yv2[~core2])[0])
    nlo, nhi = np.nanpercentile(null, [2.5, 97.5])
    p_two = 2 * min((np.sum(null <= rho_added) + 1) / (a.n_resample + 1),
                    (np.sum(null >= rho_added) + 1) / (a.n_resample + 1))
    matched = {
        f"rho_{AL}": rho_added,
        f"n_{AL}": int((~core2).sum()),
        f"matched_{CL}_rho_mean": float(np.nanmean(null)),
        f"matched_{CL}_rho_lo95": float(nlo),
        f"matched_{CL}_rho_hi95": float(nhi),
        f"{AL}_minus_matched_{CL}": rho_added - float(np.nanmean(null)),
        "p_two_sided_empirical": float(min(p_two, 1.0)),
        "n_resample": a.n_resample,
        f"deciles_with_insufficient_{CL}": str(short) if short else "なし",
    }

    # -------------------------------------------------- 3b. 1:1 最近傍マッチング
    ia = np.where(core2)[0]
    ib = np.where(~core2)[0]
    ra, rb_ = jr2.to_numpy()[ia], jr2.to_numpy()[ib]
    ma, mb = nn_match(ra, rb_, a.caliper, rng)
    sel_a, sel_b = ia[ma], ib[mb]
    nn = {}
    if len(sel_a) > 50:
        r_a = float(stats.spearmanr(xv2[sel_a], yv2[sel_a])[0])
        r_b = float(stats.spearmanr(xv2[sel_b], yv2[sel_b])[0])
        # マッチ済みペアをブートストラップして差の CI を出す
        diffs = []
        for _ in range(a.n_boot):
            k = rng.integers(0, len(sel_a), len(sel_a))
            diffs.append(stats.spearmanr(xv2[sel_b[k]], yv2[sel_b[k]])[0]
                         - stats.spearmanr(xv2[sel_a[k]], yv2[sel_a[k]])[0])
        dlo, dhi = np.nanpercentile(diffs, [2.5, 97.5])
        nn = {"n_matched_pairs": len(sel_a),
              f"rho_{CL}_matched": r_a, f"rho_{AL}_matched": r_b,
              f"delta_{AL}_minus_{CL}": r_b - r_a,
              "delta_lo95": float(dlo), "delta_hi95": float(dhi),
              f"expr_rank_median_{CL}": float(np.median(jr2.to_numpy()[sel_a])),
              f"expr_rank_median_{AL}": float(np.median(jr2.to_numpy()[sel_b])),
              "caliper": a.caliper}
        pd.Series(nn).to_frame("value").to_csv(out / "check1_nn_matched.tsv", sep="\t")

    # -------------------------------------------------- 3c. 十分位層別の重み付き平均
    ok_d = by_dec[(by_dec.n_core >= 20) & (by_dec.n_added >= 20)].copy()
    strat = {}
    if len(ok_d):
        w = np.minimum(ok_d.n_core, ok_d.n_added).to_numpy(float)
        strat = {"n_deciles_used": len(ok_d), "total_weight": float(w.sum()),
                 f"rho_{CL}_stratified": float((w * ok_d[f"rho_{CL}"]).sum() / w.sum()),
                 f"rho_{AL}_stratified": float((w * ok_d[f"rho_{AL}"]).sum() / w.sum())}
        strat[f"delta_{AL}_minus_{CL}"] = strat[f"rho_{AL}_stratified"] - strat[f"rho_{CL}_stratified"]
        strat[f"deciles_{AL}_higher"] = int((ok_d[f"rho_{AL}"] > ok_d[f"rho_{CL}"]).sum())
        pd.Series(strat).to_frame("value").to_csv(out / "check1_stratified.tsv", sep="\t")

    # -------------------------------------------------- 4. 記述的 QC
    loc_added = float(np.mean([str(g).startswith("LOC") for g in g2[~core2]])) if (~core2).any() else np.nan
    loc_core = float(np.mean([str(g).startswith("LOC") for g in g2[core2]])) if core2.any() else np.nan

    # -------------------------------------------------- 5. 判定
    r_core = float(summary.loc[summary.set == CL, "rho"].iloc[0])
    r_added = float(summary.loc[summary.set == AL, "rho"].iloc[0])
    d_raw = r_added - r_core
    d_matched = rho_added - float(np.nanmean(null))
    if rho_added < nlo:
        verdict = (f"{AL} は発現量をマッチさせても {CL} より低い（素の差 {d_raw:+.4f}、"
                   f"マッチ後 {d_matched:+.4f}）。{CL} の優位は発現量だけでは説明されない。")
    elif rho_added > nhi:
        verdict = (f"{AL} は発現量をマッチさせた {CL} より **高い**（素の差 {d_raw:+.4f}、"
                   f"マッチ後 {d_matched:+.4f}）。{CL} が優れているという読みは成り立たない。")
    else:
        verdict = (f"{AL} の rho は発現量マッチ後の {CL} の分布内に収まる"
                   f"（素の差 {d_raw:+.4f} -> マッチ後 {d_matched:+.4f}）。"
                   f"{CL} と {AL} の差は発現量の違いで説明され、集合そのものの"
                   f"性質の差とは言えない。")

    summary.to_csv(out / "check1_summary.tsv", sep="\t", index=False)
    by_dec.to_csv(out / "check1_by_decile.tsv", sep="\t", index=False)
    pd.Series(matched).to_frame("value").to_csv(out / "check1_matched.tsv", sep="\t")
    pd.DataFrame({"rho_null_matched_core": null}).to_csv(
        out / "check1_matched_null.tsv", sep="\t", index=False)

    lines = [
        f"# CHECK 1: {a.x} x {a.y}", "",
        f"- 解析遺伝子: {len(genes)}（{CL} {int(is_core.sum())} / {AL} {int((~is_core).sum())}）",
        f"- 発現量が取れた遺伝子: {len(g2)}",
        f"- LOC 接頭辞の割合: {CL} {loc_core:.3%} / {AL} {loc_added:.3%}",
        "", "## rho（95% CI）", "",
        "| 集合 | n | rho | Fisher CI | bootstrap CI | 発現量順位の中央値 |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in summary.iterrows():
        lines.append(f"| {r['set']} | {r['n']} | {r['rho']:+.4f} | "
                     f"[{r['fisher_lo']:+.4f}, {r['fisher_hi']:+.4f}] | "
                     f"[{r['boot_lo']:+.4f}, {r['boot_hi']:+.4f}] | {r['median_expr_rank']:.3f} |")
    lines += [
        "", f"素の差（{CL} − {AL}）: {unmatched['point_estimate']:+.4f} "
        f"[95% {d_lo:+.4f}, {d_hi:+.4f}]（遺伝子ブートストラップ {a.n_boot} 回）",
        "", "## 発現量マッチ後の比較", "",
        f"- {AL} の rho: {rho_added:+.4f}（n={matched[f'n_{AL}']}）",
        f"- {AL} と同じ十分位構成でリサンプリングした {CL} の rho: "
        f"平均 {matched[f'matched_{CL}_rho_mean']:+.4f}, "
        f"95% 範囲 [{nlo:+.4f}, {nhi:+.4f}]（{a.n_resample} 回）",
        f"- 差（{AL} − マッチ {CL}）: {matched[f'{AL}_minus_matched_{CL}']:+.4f}",
        f"- 経験的両側 p: {matched['p_two_sided_empirical']:.4f}",
        f"- {CL} が不足した十分位: {matched[f'deciles_with_insufficient_{CL}']}",
        "", "## 1:1 最近傍マッチング（共通サポートのみ）", ""]
    if nn:
        lines += [
            f"- マッチ成立: {nn['n_matched_pairs']} 組（caliper {nn['caliper']}）",
            f"- 発現量順位の中央値: {CL} {nn[f'expr_rank_median_{CL}']:.3f} / "
            f"{AL} {nn[f'expr_rank_median_{AL}']:.3f}",
            f"- rho: {CL} {nn[f'rho_{CL}_matched']:+.4f} / {AL} {nn[f'rho_{AL}_matched']:+.4f}",
            f"- 差（{AL} − {CL}）: {nn[f'delta_{AL}_minus_{CL}']:+.4f} "
            f"[95% {nn['delta_lo95']:+.4f}, {nn['delta_hi95']:+.4f}]"]
    else:
        lines += ["- マッチが成立しなかった（共通サポートが不足）"]
    lines += ["", "## 十分位層別（重み = 各十分位の min(n) ）", ""]
    if strat:
        lines += [
            f"- 使用した十分位: {strat['n_deciles_used']} / {N_DECILE}",
            f"- 層別重み付き rho: {CL} {strat[f'rho_{CL}_stratified']:+.4f} / "
            f"{AL} {strat[f'rho_{AL}_stratified']:+.4f}",
            f"- 差（{AL} − {CL}）: {strat[f'delta_{AL}_minus_{CL}']:+.4f}",
            f"- {AL} のほうが高かった十分位: {strat[f'deciles_{AL}_higher']} / {strat['n_deciles_used']}"]
    lines += ["", "## 判定", "", verdict,
        "", f"素の rho: {CL} {r_core:+.4f} / {AL} {r_added:+.4f} / "
        f"all {float(summary.loc[summary.set == 'all', 'rho'].iloc[0]):+.4f}",
    ]
    (out / "check1_verdict.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
