"""図 7 点の書体・級数・寸法を PDF から直接読み出して TSV にする。

pdffonts（poppler）はこの環境に入っていないので、同じ情報を PDF の
フォント辞書から読む。報告する列は pdffonts の name / type / emb / sub に対応する。
  BaseFont   -> name（先頭 6 文字 + '+' はサブセット接頭辞）
  Subtype    -> type
  FontFile2  -> emb（TrueType 実体が埋め込まれているか）
級数はページ内容ストリームの Tf 演算子から実際に使われた値を拾う。
"""
from __future__ import annotations

import re
import sys
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
FIGDIR = HERE / "manuscript" / "figures"
ORDER = ["Fig1_alpha_decomposition", "Fig2_controls_and_ceiling", "Fig3_pathway_alignment",
         "Fig4_pathway_ranking_and_aggregation", "Fig5_confounders",
         "Fig6_time_direction_amplitude", "FigS1_cosine_attenuation"]


def objects(data: bytes) -> dict[int, bytes]:
    out = {}
    for m in re.finditer(rb"(\d+)\s+0\s+obj(.*?)endobj", data, re.S):
        out[int(m.group(1))] = m.group(2)
    return out


def streams(data: bytes):
    for m in re.finditer(rb"stream\r?\n", data):
        tail = data[m.end():]
        e = tail.find(b"endstream")
        if e < 0:
            continue
        try:
            yield zlib.decompress(tail[:e])
        except Exception:
            yield tail[:e]


def fonts(data: bytes) -> list[tuple[str, str, str]]:
    objs = objects(data)
    out = []
    for body in objs.values():
        if b"/Type /Font" not in body and b"/Type/Font" not in body:
            continue
        base = re.search(rb"/BaseFont\s*/([^\s/>\]]+)", body)
        sub = re.search(rb"/Subtype\s*/([^\s/>\]]+)", body)
        if not base:
            continue
        # Type0 は合成フォントの入れ物で、実体は DescendantFonts 側にある。
        # pdffonts と同じく、実体（CIDFontType2 など）だけを数える。
        if sub and sub.group(1) == b"Type0":
            continue
        name = base.group(1).decode("latin-1")
        stype = sub.group(1).decode("latin-1") if sub else "?"
        emb = "no"
        desc = re.search(rb"/FontDescriptor\s+(\d+)\s+0\s+R", body)
        if desc:
            d = objs.get(int(desc.group(1)), b"")
            emb = "yes" if re.search(rb"/FontFile2?3?\s", d) else "no"
        out.append((name, stype, emb))
    return sorted(set(out))


def sizes(data: bytes) -> list[float]:
    found = set()
    for st in streams(data):
        for m in re.finditer(rb"/F\d+\s+([\d.]+)\s+Tf", st):
            found.add(round(float(m.group(1)), 2))
    return sorted(found)


def page_size(data: bytes) -> tuple[float, float]:
    m = re.search(rb"/MediaBox\s*\[\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)", data)
    x0, y0, x1, y1 = (float(g) for g in m.groups())
    return (x1 - x0) / 72.0, (y1 - y0) / 72.0


def images(data: bytes) -> int:
    return data.count(b"/Subtype /Image") + data.count(b"/Subtype/Image")


def main() -> int:
    rows = [("figure", "width_in", "height_in", "min_font_pt", "max_font_pt",
             "font_sizes_pt", "family_requested", "embedded_fonts", "raster_images")]
    bad = []
    for stem in ORDER:
        p = FIGDIR / f"{stem}.pdf"
        data = p.read_bytes()
        w, h = page_size(data)
        fs = sizes(data)
        fl = fonts(data)
        rows.append((stem, f"{w:.2f}", f"{h:.2f}", f"{min(fs):.1f}", f"{max(fs):.1f}",
                     " ".join(f"{v:g}" for v in fs), "Arial",
                     " ".join(f"{n}({t},emb={e})" for n, t, e in fl), str(images(data))))
        if min(fs) < 7.999:
            bad.append(f"{stem}: {min(fs)} pt の文字が残っている")
        if w > 6.301:
            bad.append(f"{stem}: 幅 {w:.2f} in が 6.3 in を超える")
        if any(e != "yes" for _n, _t, e in fl):
            bad.append(f"{stem}: 埋め込まれていないフォントがある")
        fams = {n.split("+")[-1].split("-")[0].replace("MT", "") for n, _t, _e in fl}
        if fams != {"Arial"}:
            bad.append(f"{stem}: Arial 以外の書体が入っている: {sorted(fams)}")
    out = FIGDIR.parent.parent / "results" / "revision1" / "figure_typography.tsv"
    out.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    print(out)
    for r in rows:
        print("\t".join(r))
    for b in bad:
        print("  ✗", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
