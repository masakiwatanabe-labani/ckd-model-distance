# -*- coding: utf-8 -*-
"""BD2. 補足図 S1〜S5 を現在の結果ファイルから描き直し、描画データを TSV に残す。

図は既存の描画コード（make_figures.py / figure_s1_decomposition.py / make_figures_C.py）
で描く。このスクリプトが足すのは次の 3 つ。

  1. savefig をフックして Figure オブジェクトを捕まえ、matplotlib の artist から
     「実際に描かれた座標」を取り出して TSV に落とす（結果ファイルの再計算ではなく、
     図そのものの中身）。
  2. 本文のラベル表（results/roundR/supp_labels.tsv）から番号を割り当てる。
  3. 投稿した補足図 PDF（5 ページ）と現在の PDF を文字トークンで突き合わせ、
     「同じ」「数値が変わった」「条件が変わった」を判定する。
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "code" / "revision1"))
sys.path.insert(0, str(HERE.parent / "src"))

FIG = HERE / "manuscript_C" / "figures"
DATA = HERE / "manuscript_C" / "figure_data"
SUMMARY = HERE / "results" / "roundR" / "figure_data_summary.tsv"
COMPARE = HERE / "results" / "roundR" / "figure_vs_submitted.tsv"
# 投稿済みの補足図 PDF。手元の置き場所は環境変数で与える（公開物ではない）。
import os
SUBMITTED = Path(os.environ.get("CKD_SUBMITTED_FIGURES_PDF",
                                str(Path.home() / "Downloads" /
                                    "submitted_supplementary_figures.pdf")))

# 投稿した補足図 PDF のページ順 → いまのファイル名。ページの中身で対応を取る。
SUBMITTED_PAGES = {1: "FigS2_pathway_alignment",
                   2: "FigS3_pathway_ranking_and_aggregation",
                   3: "FigS4_benchmark_simulation",
                   4: "FigS5_time_direction_amplitude",
                   5: "FigS1_mapping_and_gene_selection"}

# ラベル → 描き上がったファイルの stem。番号は supp_labels.tsv から引く。
LABEL_TO_STEM = {"figS:orthologue": "mapping_and_gene_selection",
                 "figS:pathway_alignment": "pathway_alignment",
                 "figS:pathway_rankings": "pathway_ranking_and_aggregation",
                 "figS:benchmark_sim": "benchmark_simulation",
                 "figS:time_dissimilarity": "time_direction_amplitude"}

CAPTURED: list[tuple[str, object]] = []


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------- savefig のフック
def install_hook():
    """lib_figure.savefig を包み、保存された Figure を控える。"""
    import lib_figure

    real = lib_figure.savefig

    def wrapped(fig, path, *a, **kw):
        out = real(fig, path, *a, **kw)
        CAPTURED.append((Path(path).stem, fig))
        return out

    import make_figures
    import make_figures_C
    import figure_s1_decomposition
    for mod in (lib_figure, make_figures, make_figures_C, figure_s1_decomposition):
        if getattr(mod, "savefig", None) is not None:
            setattr(mod, "savefig", wrapped)
    return wrapped


# --------------------------------------------------- artist から座標を取る
def panel_title(ax) -> str:
    for loc in ("left", "center", "right"):
        t = ax.get_title(loc=loc)
        if t:
            return t.replace("\n", " ")
    return ""


def is_reference(label: str) -> bool:
    """凡例に出ない補助線（axhline / axvline など）。"""
    return label.startswith("_")


def plot_data(fig) -> pd.DataFrame:
    """Figure に実際に描かれた座標を長い表にする。

    描画コードを読み直して再計算するのではなく、matplotlib の artist が持っている
    座標をそのまま取る。つまりこの表は図そのものの中身である。
    """
    import numpy as np
    from matplotlib.collections import LineCollection, QuadMesh
    from matplotlib.patches import PathPatch, Rectangle

    rows = []

    def add(**kw):
        base = dict(panel=0, panel_title="", artist="", series="", reference=False,
                    index=0, sub=0, x="", y="", w="", h="", value="")
        base.update(kw)
        rows.append(base)

    for pi, ax in enumerate(fig.axes):
        ttl = panel_title(ax)
        for ln in ax.lines:
            lab = ln.get_label()
            for k, (x, y) in enumerate(ln.get_xydata()):
                add(panel=pi, panel_title=ttl, artist="line", series=lab,
                    reference=is_reference(lab), index=k, x=x, y=y)
        for ci, col in enumerate(ax.collections):
            lab = col.get_label()
            if isinstance(col, QuadMesh):
                # ヒートマップ。セルの中心座標と値を出す
                A = np.ma.filled(col.get_array().astype(float), np.nan)
                co = col.get_coordinates()          # (ny+1, nx+1, 2)
                ny, nx = co.shape[0] - 1, co.shape[1] - 1
                A = np.asarray(A).reshape(ny, nx)
                for r in range(ny):
                    for c in range(nx):
                        add(panel=pi, panel_title=ttl, artist="cell", series=lab,
                            index=r, sub=c,
                            x=float(co[r:r + 2, c:c + 2, 0].mean()),
                            y=float(co[r:r + 2, c:c + 2, 1].mean()),
                            value=A[r, c])
                continue
            if isinstance(col, LineCollection):
                for si, seg in enumerate(col.get_segments()):
                    for k, (x, y) in enumerate(np.asarray(seg)):
                        add(panel=pi, panel_title=ttl, artist="segment", series=lab,
                            reference=is_reference(lab), index=si, sub=k, x=x, y=y)
                continue
            try:
                off = np.asarray(col.get_offsets(), dtype=float)
            except Exception:
                continue
            for k, pt in enumerate(off):
                if pt.size != 2:
                    continue
                add(panel=pi, panel_title=ttl, artist="points", series=lab,
                    reference=is_reference(lab), index=k, x=pt[0], y=pt[1])
        for si, pa in enumerate(ax.patches):
            lab = pa.get_label()
            if isinstance(pa, Rectangle):
                x, y = pa.get_xy()
                add(panel=pi, panel_title=ttl, artist="bar", series=lab,
                    reference=is_reference(lab), index=si, x=x, y=y,
                    w=pa.get_width(), h=pa.get_height())
            elif isinstance(pa, PathPatch):
                # バイオリンの外形。輪郭の頂点をそのまま出す
                for k, (x, y) in enumerate(np.asarray(pa.get_path().vertices, dtype=float)):
                    add(panel=pi, panel_title=ttl, artist="outline", series=lab,
                        reference=is_reference(lab), index=si, sub=k, x=x, y=y)
        for si, im in enumerate(ax.images):
            A = np.ma.filled(np.asarray(im.get_array(), dtype=float), np.nan)
            A = A.reshape(A.shape[0], -1) if A.ndim > 2 else A
            for r in range(A.shape[0]):
                for c in range(A.shape[1]):
                    add(panel=pi, panel_title=ttl, artist="cell", series=im.get_label(),
                        index=r, sub=c, value=A[r, c])
        for si, tx in enumerate(ax.texts):
            s_ = tx.get_text().replace("\n", " ")
            if not s_.strip():
                continue
            x, y = tx.get_position()
            add(panel=pi, panel_title=ttl, artist="text", series=s_, index=si, x=x, y=y)
        # 目盛りラベルも図の中身なので残す（状態名・経路名の順序が分かる）
        for axis, nm in ((ax.xaxis, "xtick"), (ax.yaxis, "ytick")):
            labs = [t.get_text().replace("\n", " ") for t in axis.get_ticklabels()]
            locs = list(axis.get_ticklocs())
            for k, (lo, la) in enumerate(zip(locs, labs)):
                if la.strip():
                    add(panel=pi, panel_title=ttl, artist=nm, series=la, index=k,
                        x=lo if nm == "xtick" else "", y="" if nm == "xtick" else lo)
    return pd.DataFrame(rows)


# ------------------------------------------------ 投稿版との突き合わせ
TOKEN = re.compile(r"[A-Za-z][A-Za-z/\-']*|\d+(?:\.\d+)?")


def tokens(text: str):
    """文字を単語と数値に割る。"""
    words, nums = [], []
    for t in TOKEN.findall(text):
        (nums if t[0].isdigit() else words).append(t)
    return words, nums


def pdf_text(path: Path, page: int | None = None) -> str:
    from pypdf import PdfReader
    r = PdfReader(str(path))
    pages = r.pages if page is None else [r.pages[page - 1]]
    return "\n".join((p.extract_text() or "") for p in pages)


def compare_to_submitted() -> pd.DataFrame:
    from collections import Counter
    rows = []
    if not SUBMITTED.exists():
        print(f"  投稿版 PDF が見つからない: {SUBMITTED}")
        return pd.DataFrame(rows)
    for page, stem in sorted(SUBMITTED_PAGES.items()):
        cur = FIG / f"{stem}.pdf"
        ow, on = tokens(pdf_text(SUBMITTED, page))
        nw, nn = tokens(pdf_text(cur))
        dw = Counter(nw) - Counter(ow), Counter(ow) - Counter(nw)
        dn = Counter(nn) - Counter(on), Counter(on) - Counter(nn)
        words_same = not (dw[0] or dw[1])
        nums_same = not (dn[0] or dn[1])
        verdict = ("same" if words_same and nums_same else
                   "numbers changed" if words_same else "conditions changed")
        rows.append(dict(submitted_page=page, figure_file=cur.name, verdict=verdict,
                         words_submitted=len(ow), words_current=len(nw),
                         numbers_submitted=len(on), numbers_current=len(nn),
                         words_added="|".join(sorted(dw[0])[:12]),
                         words_removed="|".join(sorted(dw[1])[:12]),
                         numbers_added="|".join(sorted(dn[0])[:12]),
                         numbers_removed="|".join(sorted(dn[1])[:12])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main
def main() -> int:
    L = pd.read_csv(HERE / "results" / "roundR" / "supp_labels.tsv", sep="\t")
    num = {r.label: int(r.number) for r in L.itertuples() if r.kind == "Figure"}
    want = {f"FigS{num[k]}_{v}" for k, v in LABEL_TO_STEM.items()}

    before = {p.name: sha256(p) for p in sorted(FIG.glob("Fig*"))}

    install_hook()
    import make_figures_C as MFC
    MFC.main()

    # 捕まえた Figure を C案の名前に対応づける（REUSE の付け替えをたどる）
    rename = dict(MFC.REUSE)
    figs: dict[str, object] = {}
    for stem, fig in CAPTURED:
        figs[rename.get(stem, stem)] = fig

    missing = want - set(figs)
    assert not missing, f"描画を捕まえられなかった図: {sorted(missing)}"

    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    for lab, tail in LABEL_TO_STEM.items():
        stem = f"FigS{num[lab]}_{tail}"
        png = FIG / f"{stem}.png"
        assert png.exists(), f"{png} がない"
        D = plot_data(figs[stem])
        tsv = DATA / f"{stem}_data.tsv"
        D.to_csv(tsv, sep="\t", index=False, float_format="%.6g")
        rows.append(dict(label=lab, figure=f"Figure S{num[lab]}", file=png.name,
                         data_file=tsv.name, panels=int(D.panel.nunique()),
                         data_rows=len(D),
                         sha256_before=before.get(png.name, ""),
                         sha256_after=sha256(png),
                         changed=before.get(png.name, "") != sha256(png)))
    T = pd.DataFrame(rows).sort_values("figure")
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(SUMMARY, sep="\t", index=False)

    C = compare_to_submitted()
    if len(C):
        C.to_csv(COMPARE, sep="\t", index=False)

    pd.set_option("display.width", 240)
    print("\n### 描き直した補足図と描画データ")
    print(T[["figure", "file", "data_file", "panels", "data_rows", "changed"]].to_string(index=False))
    if len(C):
        print("\n### 投稿版の補足図 PDF との突き合わせ")
        print(C[["submitted_page", "figure_file", "verdict",
                 "words_added", "words_removed", "numbers_added", "numbers_removed"]]
              .to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
