# -*- coding: utf-8 -*-
"""AL6: 投稿ファイル一覧（ファイル名・用途・SHA-256）を作る。"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
OUT = Path(__file__).resolve().parents[2] / "results" / "roundR" / "submission_manifest.tsv"
ROLE = {
    "Main_manuscript.docx": "Manuscript, clean copy",
    "Main_manuscript_tracked.docx": "Manuscript, all changes marked",
    "Main_manuscript_blue.docx": "Manuscript, all changes accepted with the revised text in blue",
    "Response_to_reviewers.docx": "Point-by-point response",
    "Cover_letter.docx": "Cover letter to the editors",
    "Supplementary_Notes.docx": "Supplementary Notes S1-S5, Supplementary Figures S1-S5 and all captions",
    "Supplementary_Figures.pdf": "Supplementary Figures S1-S5, one page per figure",
    "Supplementary_Tables.xlsx": "Supplementary Tables S1-S12 and the Table S10 gene list",
    "Supplementary_Data_S1_PodTRECK_proteome.xlsx":
        "Supplementary Data S1: Pod-TRECK quantification table used to define Group A",
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    rows, missing, nonascii = [], [], []
    for name, role in ROLE.items():
        p = DEST / name
        if not p.exists():
            missing.append(name)
            continue
        if not name.isascii():
            nonascii.append(name)
        rows.append((name, role, f"{p.stat().st_size:,}", sha256(p)))
    extra = sorted(q.name for q in DEST.iterdir()
                   if q.is_file() and q.name not in ROLE
                   and not q.name.startswith((".", "~$")))
    w = max(len(r[0]) for r in rows)
    print(f"{'file'.ljust(w)}  {'bytes'.rjust(11)}  sha256")
    for n, role, b, h in rows:
        print(f"{n.ljust(w)}  {b.rjust(11)}  {h}")
        print(f"{''.ljust(w)}  {role}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("file\trole\tbytes\tsha256\n")
        for n, role, b, h in rows:
            fh.write(f"{n}\t{role}\t{b.replace(',', '')}\t{h}\n")
    print(f"\n書き出し: {OUT}")
    if missing:
        print(f"未作成: {missing}")
    if extra:
        print(f"一覧に無いファイル: {extra}")
    if nonascii:
        print(f"✗ ファイル名に非 ASCII: {nonascii}")
    return 1 if nonascii else 0


if __name__ == "__main__":
    sys.exit(main())
