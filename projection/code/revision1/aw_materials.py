# -*- coding: utf-8 -*-
"""AW1: Group A の再計算に必要な材料の出所を確定して記録する。

Group A = ネコ皮質の検出 ∩ ネコ髄質の検出 ∩ Pod-TRECK の検出（ヒト記号空間で交差）
なので、Supplementary Data S1 だけでは Pod-TRECK 側の判定しか出ない。
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
OUT = HERE / "results" / "roundR" / "groupA_materials.tsv"
# Ensembl のリリース日（FTP の README の Last-Modified で確認）
ENSEMBL = {"115": "2025-09-03", "116": "2026-09-10"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def mtime(p: Path) -> str:
    import datetime
    return datetime.date.fromtimestamp(p.stat().st_mtime).isoformat()


def main() -> int:
    cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    cat = ROOT / cfg["paths"]["raw"] / cfg["files"]["cat_omics"]
    pod = ROOT / cfg["paths"]["raw"] / cfg["files"]["mouse_prot"]
    o1 = ROOT / "data" / "ref" / "orthologs_cat2human.tsv"
    o2 = ROOT / "data" / "ref" / "orthologs_mouse2human.tsv"
    rows = []

    def add(role, p, note=""):
        rows.append({"role": role, "file": p.name, "bytes": p.stat().st_size,
                     "downloaded": mtime(p), "sha256": sha256(p), "note": note})

    add("feline cortical and medullary proteomes",
        cat, "Supplementary Data 3 of Li et al., Commun. Biol. 2025;8:1794 "
             "(doi:10.1038/s42003-025-09164-8), sheets 'S5 Protein cortex' and "
             "'S6 Protein medulla'; publicly downloadable from the publisher")
    # 受託解析の元の表は案件番号を名前に持つので、公開する一覧では添付した
    # Supplementary Data S1 の名前とチェックサムで参照する（中身は同一）。
    att = Path.home() / "Desktop" / "CKD_xspecies_submission" / \
        "Supplementary_Data_S1_PodTRECK_proteome.xlsx"
    if att.exists():
        rows.append({"role": "Pod-TRECK proteome", "file": att.name,
                     "bytes": att.stat().st_size, "downloaded": mtime(pod),
                     "sha256": sha256(att),
                     "note": "attached as Supplementary Data S1 of this paper; "
                             "the same values as the source table, with the document "
                             "properties cleared"})
    else:
        add("Pod-TRECK proteome", pod, "attached as Supplementary Data S1 of this paper")
    # オルソログ表の取得日から、当時の現行リリースを決める
    d = mtime(o1)
    rel = max((r for r, day in ENSEMBL.items() if day <= d), key=int, default="")
    add("cat-to-human one-to-one orthologues", o1,
        f"Ensembl BioMart, www.ensembl.org; current release on the download date was {rel} "
        f"(release 116 was published {ENSEMBL['116']})")
    add("mouse-to-human one-to-one orthologues", o2,
        f"Ensembl BioMart, www.ensembl.org; current release {rel}")
    rows.append({"role": "ensembl release at download", "file": "", "bytes": "",
                 "downloaded": d, "sha256": "", "note": rel})
    T = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 46)
    print(T[["role", "file", "bytes", "downloaded", "note"]].to_string(index=False))
    print(f"\nEnsembl リリース（取得日 {d} 時点の現行版）: {rel}")
    print(f"書き出し: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
