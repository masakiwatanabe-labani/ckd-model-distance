# -*- coding: utf-8 -*-
"""BH: 変更履歴版から「修正箇所を青字にした版」を作り、検査して docx に詰める。

  build_blue.py <変更履歴版 document.xml> <出力 document.xml>

道具は Masaki の ~/Desktop/tracked_tools/build_blue.py を優先し、無ければ
同じ仕様で書いた code/revision1/bh_build_blue.py を使う。容れ物は修正稿の docx。
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
TRACKED = DEST / "Main_manuscript_tracked.docx"
OUT = DEST / "Main_manuscript_blue.docx"
WORK = HERE / "results" / "roundR" / "tracked_build"
OWN = HERE / "code" / "revision1" / "bh_build_blue.py"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    if not TRACKED.exists():
        print(f"変更履歴版が無い: {TRACKED}", file=sys.stderr)
        return 2
    WORK.mkdir(parents=True, exist_ok=True)
    shutil.copy2(TOOLS / "blk.py", WORK / "blk.py")

    given = TOOLS / "build_blue.py"
    if given.exists():
        shutil.copy2(given, WORK / "build_blue.py")
        tool, origin = WORK / "build_blue.py", f"{given}（Masaki）"
    else:
        shutil.copy2(OWN, WORK / "build_blue.py")
        tool, origin = WORK / "build_blue.py", f"{OWN.name}（仕様どおりの代替。Masaki の版が置かれたらそちらを使う）"
    print(f"道具: {origin}")

    (WORK / "tracked.xml").write_bytes(zipfile.ZipFile(TRACKED).read("word/document.xml"))
    r = subprocess.run([sys.executable, str(tool), "tracked.xml", "blue.xml"],
                       capture_output=True, text=True, cwd=WORK)
    if r.returncode != 0:
        print(r.stdout + r.stderr, file=sys.stderr)
        return 1
    print(" ", r.stdout.strip())

    blue = (WORK / "blue.xml").read_bytes()
    zin = zipfile.ZipFile(REV)
    tmp = OUT.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            data = blue if info.filename == "word/document.xml" else zin.read(info.filename)
            zout.writestr(info, data)
    shutil.move(str(tmp), str(OUT))
    print(f"\n{OUT.name}: {OUT.stat().st_size:,} bytes")
    print(f"  sha256 {sha256(OUT)}")
    print(f"  容れ物 {REV.name} sha256 {sha256(REV)}")
    print(f"  もと   {TRACKED.name} sha256 {sha256(TRACKED)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
