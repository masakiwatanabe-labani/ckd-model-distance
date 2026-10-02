# -*- coding: utf-8 -*-
"""AT7: build_tracked.py が作った document.xml を docx に詰め、検査する。

容れ物は修正稿の docx（媒体・スタイル・rels がそろっている）。document.xml だけ差し替える。
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
TOOLS = Path.home() / "Desktop" / "tracked_tools"
DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
REV = DEST / "Main_manuscript.docx"
BASE = HERE / "data" / "base" / "base_manuscript.docx"
OUT = DEST / "Main_manuscript_tracked.docx"
WORK = HERE / "results" / "roundR" / "tracked_build"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(*args) -> str:
    r = subprocess.run([sys.executable, *[str(a) for a in args]], capture_output=True,
                       text=True, cwd=WORK)
    if r.returncode != 0:
        print(r.stdout + r.stderr, file=sys.stderr)
        raise SystemExit(f"失敗: {args[0]}")
    return r.stdout.strip()


def main() -> int:
    for f in ("build_tracked.py", "check_tracked.py", "blk.py"):
        if not (TOOLS / f).exists():
            print(f"道具が無い: {TOOLS / f}", file=sys.stderr)
            return 2
    WORK.mkdir(parents=True, exist_ok=True)
    for f in ("build_tracked.py", "check_tracked.py", "blk.py"):
        shutil.copy2(TOOLS / f, WORK / f)
    (WORK / "base.xml").write_bytes(zipfile.ZipFile(BASE).read("word/document.xml"))
    (WORK / "rev.xml").write_bytes(zipfile.ZipFile(REV).read("word/document.xml"))
    # 色だけの変更は文字の変更ではないので、変更履歴には記録しない。
    # build_tracked.py には unchanged として渡す（revised 側のブロックがそのまま入る）。
    import pandas as pd
    BM = pd.read_csv(HERE / "results" / "roundR" / "an_block_map.tsv", sep="\t",
                     dtype=str, keep_default_na=False)
    n_fmt = int((BM.action == "recolored (format only)").sum())
    BM.loc[BM.action == "recolored (format only)", "action"] = "unchanged"
    BM.to_csv(WORK / "an_block_map.tsv", sep="\t", index=False)
    if n_fmt:
        print(f"書式のみの変更 {n_fmt} 件は unchanged として渡した（変更履歴に残さない）")

    print(run(WORK / "build_tracked.py", "base.xml", "rev.xml", "an_block_map.tsv", "tracked.xml"))
    chk = run(WORK / "check_tracked.py", "tracked.xml", "rev.xml", "base.xml")
    print(chk)
    if "accept-all == revised: OK" not in chk or "reject-all == base   : OK" not in chk:
        print("✗ check_tracked.py が一致を報告していない", file=sys.stderr)
        return 1

    # 修正稿の docx を容れ物にして document.xml だけ差し替える
    tracked = (WORK / "tracked.xml").read_bytes()
    zin = zipfile.ZipFile(REV)
    tmp = OUT.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            data = tracked if info.filename == "word/document.xml" else zin.read(info.filename)
            zout.writestr(info, data)
    shutil.move(str(tmp), str(OUT))
    print(f"\n{OUT.name}: {OUT.stat().st_size:,} bytes")
    print(f"  sha256 {sha256(OUT)}")
    print(f"  容れ物 {REV.name} sha256 {sha256(REV)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
