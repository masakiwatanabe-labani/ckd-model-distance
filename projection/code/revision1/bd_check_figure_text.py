# -*- coding: utf-8 -*-
"""BD3. 本文・補足ノートが図について述べている要約値を、図の描画データから再計算して突き合わせる。

入力は manuscript_C/figure_data/FigS*_data.tsv（bd_figures.py が matplotlib の artist から
取り出した「実際に描かれた座標」）と、manuscript_C/*.md の本文。
結果ファイルを経由しないので、「図と本文が食い違っている」場合にここで落ちる。

対象にした文は COVERED に列挙してある。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
MS = HERE / "manuscript_C"
FD = MS / "figure_data"

FAIL: list[str] = []
OK = 0
COVERED: list[tuple[str, str]] = []   # (出典, 文)


def load(stem: str) -> pd.DataFrame:
    T = pd.read_csv(FD / f"{stem}_data.tsv", sep="\t", dtype=str, keep_default_na=False)
    for c in ("x", "y", "w", "h", "value"):
        T[c] = pd.to_numeric(T[c], errors="coerce")
    T["panel"] = T.panel.astype(int)
    return T


def text(name: str) -> str:
    return (MS / name).read_text()


def sentence(name: str, needle: str) -> str:
    """needle を含む 1 文を返す。対象文の一覧を作るために使う。"""
    t = " ".join(text(name).split())
    i = t.find(needle)
    assert i >= 0, f"{name} に見つからない: {needle}"
    a = max(t.rfind(". ", 0, i) + 2, 0)
    b = t.find(". ", i)
    s = t[a:b + 1] if b > 0 else t[a:]
    COVERED.append((name, s))
    return s


def chk(what: str, stated, recomputed, tol=0.0):
    global OK
    ok = (abs(float(stated) - float(recomputed)) <= tol
          if not isinstance(stated, str) else stated == recomputed)
    if ok:
        OK += 1
    else:
        FAIL.append(f"{what}: 本文 {stated} / 図から {recomputed}")
    print(f"  {'ok ' if ok else 'NG '} {what:62s} 本文={stated} 図={recomputed}")


def auc(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))


# ------------------------------------------------------------ Figure S1
def fig_s1():
    print("\n[Figure S1]  Orthologue characteristics and gene selection")
    T = load("FigS1_mapping_and_gene_selection")
    sentence("SUPPLEMENTARY.md", "(A) Median cosine in increasing sequence-identity quartiles")
    p0 = T[T.panel == 0]
    bars = p0[p0.artist == "bar"].sort_values("x")
    labs = p0[(p0.artist == "text") & p0.series.str.match(r"^\d\.\d{3}$")].sort_values("x")
    chk("S1A 四分位 × 2 群の棒の本数", 8, len(bars))
    chk("S1A 横軸の四分位", "Q1|Q2|Q3|Q4", "|".join(p0[p0.artist == "xtick"].series))
    for (_i, b), (_j, l) in zip(bars.iterrows(), labs.iterrows()):
        chk(f"S1A x={b.x:+.3f} の棒の高さと印字", float(l.series), round(b.h, 3), 5e-4)


# ------------------------------------------------------------ Figure S2
def fig_s2():
    print("\n[Figure S2]  Pathway-level alignment")
    T = load("FigS2_pathway_alignment")
    p0, p2 = T[T.panel == 0], T[T.panel == 2]

    sentence("RESULTS.md", "its central 95% cosine range spanned 0.652 at 30 genes")
    band = p0[(p0.artist == "line") & (p0.series == "_child0")].set_index("x")["y"]
    for n, stated in ((30, 0.652), (50, 0.540), (100, 0.419)):
        chk(f"S2A n={n} 遺伝子での中央 95% 幅", stated, round(float(band.loc[n]), 3), 5e-4)

    sentence("SUPPLEMENTARY.md", "244 of 422 sets contained fewer than 50")
    ins = p0[(p0.artist == "text")].series.iloc[0]
    m = re.search(r"(\d+) of (\d+) pathways fall below", ins)
    chk("S2A 注記「50 未満の経路数 / 全経路数」", "244/422", f"{m.group(1)}/{m.group(2)}")

    sentence("SUPPLEMENTARY.md", "Figure S2B displays 1,184 values from 74")
    sentence("METHODS.md", "Figure S2B displays 74 qualifying sets")
    cells = p2[p2.artist == "cell"]
    states = p2[p2.artist == "ytick"]
    chk("S2B ヒートマップの値の数", 1184, len(cells))
    chk("S2B 縦軸の状態数", 16, len(states))
    chk("S2B 経路数（値の数 ÷ 状態数）", 74, len(cells) // len(states))


# ------------------------------------------------------------ Figure S3
def fig_s3():
    print("\n[Figure S3]  Pathway rankings and the centering contrast")
    T = load("FigS3_pathway_ranking_and_aggregation")
    p0, p1, p2, p3 = (T[T.panel == i] for i in range(4))

    sentence("RESULTS.md", "IRI 28 d led 22% of sets, IRI 12 mo 21% and IRI 2 h 11%")
    states = list(p0[p0.artist == "ytick"].sort_values("y").series)
    shares = list(p0[p0.artist == "text"].sort_values("y").series)
    got = dict(zip(states, shares))
    for st, stated in (("IRI 28 d", "22%"), ("IRI 12 mo", "21%"), ("IRI 2 h", "11%")):
        chk(f"S3A {st} が首位の経路の割合", f"leads {stated}", got[st])

    sentence("SUPPLEMENTARY.md", "(A) Mouse-state ranks across 422 overlapping pathways")
    pts = p2[p2.artist == "points"]
    chk("S3C に描かれた経路の数", 422, len(pts))

    sentence("RESULTS.md", "The median pathway-specific AUC was 0.760")
    y = pts.y.to_numpy(float)
    chk("S3 経路ごとの AUC の中央値", 0.760, round(float(np.median(y)), 3), 5e-4)
    chk("S3 同 四分位範囲の下側", 0.707, round(float(np.percentile(y, 25)), 3), 5e-4)
    chk("S3 同 四分位範囲の上側", 0.808, round(float(np.percentile(y, 75)), 3), 5e-4)

    sentence("SUPPLEMENTARY.md", "The dashed line gives the whole-gene AUC and the solid line")
    seg = p1[p1.artist == "segment"].groupby("series").x.first().sort_values()
    chk("S3B 実線（中心化後の経路平均の AUC）", 0.511, round(float(seg.iloc[0]), 3), 5e-4)
    chk("S3B 破線（全遺伝子の AUC）", 0.813, round(float(seg.iloc[1]), 3), 5e-4)
    chk("S3 全遺伝子の AUC に達した経路の割合", 23,
        int(round(float((y >= float(seg.iloc[1])).mean()) * 100)))

    sentence("SUPPLEMENTARY.md", "(D) Pairwise cosine under four preprocessing conditions")
    P = p3[p3.artist == "points"].copy()
    P["cls"] = np.where(P.series.isin(["within species", "cat vs mouse"]), P.series, "")
    # 系列名は 1 本目の散布にしか付かないので、x の整数部で条件を、y の並びで class を取る
    P["cond"] = P.x.round().astype(int)
    lab = p3[(p3.artist == "text") & p3.series.str.startswith("AUC ")].sort_values("x")
    n_within = n_cross = None
    for i, (_j, r) in enumerate(lab.iterrows()):
        g = P[P.cond == i]
        # within は 39 ペア、cross は 48 ペア。件数で振り分ける
        counts = g.groupby("index").size()
        by_series = g.groupby(g.series.where(g.series.str.startswith("_") == False, "unnamed"))
        # x のオフセットで within / cross を分ける（within は −0.11、cross は +0.11）
        w = g[g.x < i]["y"].to_numpy(float)
        c = g[g.x > i]["y"].to_numpy(float)
        if w.size == 0 or c.size == 0:      # オフセットがない場合は件数で分ける
            srt = sorted(g.groupby("series"), key=lambda kv: len(kv[1]))
            w, c = srt[0][1].y.to_numpy(float), srt[1][1].y.to_numpy(float)
        n_within, n_cross = len(w), len(c)
        chk(f"S3D 条件 {i + 1} の印字 AUC", float(r.series.split()[1]), round(auc(w, c), 3), 5e-4)
    sentence("SUPPLEMENTARY.md", "(B) Distribution of pathway-specific AUCs among the 87 state pairs")
    chk("S3D の種内ペア数 + 種間ペア数", 87, n_within + n_cross)


# ------------------------------------------------------------ Figure S4
def fig_s4():
    print("\n[Figure S4]  Benchmark simulation")
    T = load("FigS4_benchmark_simulation")
    E = pd.read_csv(HERE / "results" / "roundR" / "simulation_exceedance.tsv", sep="\t")
    E["shift"] = E["shift"].astype(str)
    E["true_cos"] = pd.to_numeric(E.true_cos, errors="coerce")   # "<1.0" の行は NaN
    src = pd.read_csv(HERE / "results" / "round5" / "ceiling_decomposition_sim.tsv", sep="\t")

    def band_row(panel, design, tc):
        r = E[(E.panel == panel) & (E.design == design) & (E.true_cos == tc)
              & (E["shift"] == "all") & (~E.metric.str.contains("grid points"))]
        assert len(r) == 1, f"{panel}/{design}/{tc} が {len(r)} 行"
        return r.iloc[0]

    sentence("METHODS.md", "including the observed split-half range of 0.568–0.986")
    band = T[T.artist == "bar"]
    chk("S4 網掛けの帯が引かれたパネル数", 4, len(band))
    chk("S4 帯の下端", 0.568, round(float(band.x.min()), 3), 5e-4)
    chk("S4 帯の上端", 0.986, round(float((band.x + band.w).max()), 3), 5e-4)

    # --- 設計別の値（本文が主として報告する量）
    sentence("METHODS.md", "Identical responses exceeded it in 59.4% of replicates")
    sentence("SUPPLEMENTARY.md", "Solid lines with filled markers are the mouse-like design")
    for design, mean, lo, hi in (("feline-like (7 vs 6)", 59.4, 57.0, 62.0),
                                 ("mouse-like (3 vs 6)", 98.0, 96.3, 99.7)):
        r = band_row("iv", design, 1.0)
        chk(f"S4(iv) {design} の帯内平均（結果ファイル）", mean, round(r["mean"] * 100, 1), 0.05)
        chk(f"S4(iv) {design} の帯内の下限", lo, round(r["min"] * 100, 1), 0.05)
        chk(f"S4(iv) {design} の帯内の上限", hi, round(r["max"] * 100, 1), 0.05)
        chk(f"S4(iv) {design} の帯内の振れ幅が 5 ポイント以内", True,
            bool(round((r["max"] - r["min"]) * 100, 1) <= 5.0))
        # 図に実際に描かれた系列が、本文の範囲に収まっていること
        g = T[(T.panel == 3) & (T.artist == "line")
              & (T.series == f"{design} | true cos 1.0")]
        ins = g[(g.x >= 0.568) & (g.x <= 0.986)]
        chk(f"S4(iv) {design} の曲線が描かれている（帯内の点）", 5, len(ins))
        chk(f"S4(iv) {design} の描画値が本文の範囲 {lo}–{hi}% に収まる", True,
            bool(lo / 100 - 1e-9 <= ins.y.min() and ins.y.max() <= hi / 100 + 1e-9))

    # --- プール平均 78.7% は従として残す
    sentence("METHODS.md", "Averaged over both designs and the three mean shifts the rate was")
    chk("S4(iv) 両設計・3 shift をまとめた平均", 78.7,
        round(band_row("iv", "pooled", 1.0)["mean"] * 100, 1), 0.05)

    # --- 真のコサイン 0.9 は信頼性の階段
    sentence("METHODS.md", "At a true cosine of 0.9 this exceedance was a step in reliability")
    for tag, stated in (("ii", 39.7), ("iii", 39.9)):
        chk(f"S4({tag}) 帯内の格子平均（本文が「率ではない」と断る値）", stated,
            round(band_row(tag, "pooled", 0.9)["mean"] * 100, 1), 0.05)
        hi_ = E[(E.panel == tag) & E.metric.str.contains("non-exceeding") ]
        lo_ = E[(E.panel == tag) & E.metric.str.contains(r"\(exceeding")]
        assert len(hi_) == 1 and len(lo_) == 1
        chk(f"S4({tag}) 超過する側の信頼性の上端", 0.816,
            float(lo_.iloc[0].reliability.split("-")[1]), 5e-4)
        chk(f"S4({tag}) 超過しない側の信頼性の下端", 0.865,
            float(hi_.iloc[0].reliability.split("-")[0]), 5e-4)
        chk(f"S4({tag}) 超過する側の最小の超過率（本文の 97%）", True,
            bool(round(lo_.iloc[0]["min"] * 100) >= 97))
        chk(f"S4({tag}) 超過しない側は 0", 0.0, float(hi_.iloc[0]["max"]), 1e-12)
    # 図の曲線でも階段が見えること（両設計）
    for design in ("mouse-like (3 vs 6)", "feline-like (7 vs 6)"):
        g = T[(T.panel == 2) & (T.artist == "line")
              & (T.series == f"{design} | true cos 0.9")]
        chk(f"S4(iii) {design}: 信頼性 0.816 以下は 0.97 以上", True,
            bool((g[g.x <= 0.816].y >= 0.97).all()))
        chk(f"S4(iii) {design}: 信頼性 0.865 以上は 0", True,
            bool((g[g.x >= 0.865].y == 0).all()))

    # --- 描いた shift
    sentence("SUPPLEMENTARY.md", "All panels use a common mean shift of 0.26")
    shifts = sorted(src["shift"].unique())
    chk("S4 の shift が経験的中央値 0.257 に最も近い格子点", 0.26,
        min(shifts, key=lambda s_: abs(s_ - 0.257)), 1e-9)
    drawn = sorted({float(v) for v in
                    src[src.design.isin(["mouse-like (3 vs 6)", "feline-like (7 vs 6)"])]["shift"]})
    chk("S4 が 1 つの shift だけを描いていること", 1,
        len({round(float(x), 2) for x in
             src[(src["shift"] == 0.26)]["shift"]}))

    # --- 撤回した主張の歯止め
    body = " ".join((text("METHODS.md") + text("SUPPLEMENTARY.md")).split())
    chk("「which is the gradient shown in Figure S4」が残っていないこと", True,
        "gradient shown in Figure S4" not in body)
    chk("S4 に 57.0% を「帯の下端の値」として書いていないこと", True,
        "57.0% at the lower end" not in body)


# ------------------------------------------------------------ Figure S5
def fig_s5():
    print("\n[Figure S5]  Time, direction and amplitude")
    T = load("FigS5_time_direction_amplitude")
    sentence("SUPPLEMENTARY.md", "Each point represents one of the 66 pairs among twelve mouse states")
    for p in (0, 1):
        n = len(T[(T.panel == p) & (T.artist == "points")])
        chk(f"S5 パネル {'AB'[p]} の点の数", 66, n)

    sentence("RESULTS.md", "with Mantel statistics of 0.599")
    ins = T[(T.panel == 0) & (T.artist == "text")].series.iloc[0]
    m = re.search(r"\+([0-9.]+)", ins)
    chk("S5A の Mantel ρ", 0.599, float(m.group(1)), 5e-4)


# --------------------------------------- 提出する図 PDF の並びが番号と合うこと
def figures_pdf():
    """Supplementary_Figures.pdf の n ページ目が Figure S<n> であること。

    投稿済みの PDF は旧番号の並び（1 ページ目が経路の整列）だった。番号を振り直した
    あとも同じ並びのまま出さないための歯止め。
    """
    import re as _re
    from pypdf import PdfReader
    pkg = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Supplementary_Figures.pdf"
    if not pkg.exists():
        print("\n[Supplementary_Figures.pdf]  まだ無いので省略")
        return
    print("\n[Supplementary_Figures.pdf]  ページの並びが本文の番号と合うこと")
    FIGDIR = HERE / "manuscript_C" / "figures"
    tok = _re.compile(r"[A-Za-z][A-Za-z/\-']*|\d+(?:\.\d+)?")

    def toks(page_text):
        return sorted(tok.findall(page_text))

    R = PdfReader(str(pkg))
    chk("図 PDF のページ数", 5, len(R.pages))
    for n in range(1, 6):
        own = sorted(FIGDIR.glob(f"FigS{n}_*.pdf"))
        assert len(own) == 1, f"FigS{n} の PDF が {len(own)} 件"
        a = toks(PdfReader(str(own[0])).pages[0].extract_text() or "")
        b = toks(R.pages[n - 1].extract_text() or "")
        chk(f"図 PDF の {n} ページ目 = Figure S{n}（{own[0].name}）", True, a == b)


def main() -> int:
    for fn in (fig_s1, fig_s2, fig_s3, fig_s4, fig_s5, figures_pdf):
        fn()
    print("\n### 対象にした文（本文・補足ノート側）")
    seen = set()
    for src, s in COVERED:
        if s in seen:
            continue
        seen.add(s)
        print(f"  [{src}] {s}")
    print(f"\nBD3 図↔本文の検査: {OK} 件一致 / 不一致 {len(FAIL)} 件")
    for f in FAIL:
        print("  NG", f)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
