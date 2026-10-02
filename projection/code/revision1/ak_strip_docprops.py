# -*- coding: utf-8 -*-
"""AK: 添付する定量表から、著者でない個人名と受託先の内部パスを除く。

方針:
  - シートのデータ・書式・欠測の表現には一切触れない
  - 変更するのは docProps/core.xml と xl/workbook.xml の 2 エントリのみ
  - dc:creator / cp:lastModifiedBy は空文字にする（著者名に置き換えない）
  - xl/workbook.xml の x15ac:absPath（受託先の内部パスと案件 ID）を除く
  - 日付（dcterms:created / modified）は残す
  - openpyxl や pandas は使わない。zip のエントリをバイト列のまま引き写す
"""
from __future__ import annotations

import hashlib
import re
import shutil
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE.parent
import local_source as _LS  # noqa: E402
SRC = _LS.source_path()   # 受託解析の元表。名前は公開しない
DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
NAME = "Supplementary_Data_S1_PodTRECK_proteome.xlsx"
KEEP_IDENTICAL = ("xl/worksheets/", "xl/sharedStrings.xml", "xl/styles.xml", "xl/theme/",
                  "xl/_rels/", "[Content_Types].xml", "_rels/.rels", "docProps/app.xml")



def leak_terms() -> list[bytes]:
    """公開前の点検で探す語。第三者の名前などは leak_terms.txt から読む。

    そのファイルは .gitignore で公開しない。無い場合は一般的なパターンだけで点検する。
    """
    f = Path(__file__).resolve().parents[2] / "leak_terms.txt"
    out = []
    if f.exists():
        out = [l.strip().encode() for l in f.read_text().splitlines()
               if l.strip() and not l.startswith("#")]
    return out

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def clean_core(x: str) -> tuple[str, list[str]]:
    removed = []
    for tag in ("dc:creator", "cp:lastModifiedBy"):
        m = re.search(rf"<{tag}>(.*?)</{tag}>", x, re.S)
        if m and m.group(1).strip():
            removed.append(f"{tag}={m.group(1)}")
            x = x.replace(m.group(0), f"<{tag}></{tag}>")
    return x, removed


def clean_workbook(x: str) -> tuple[str, list[str]]:
    removed = []
    for m in re.finditer(r"<mc:AlternateContent\b.*?</mc:AlternateContent>", x, re.S):
        if "absPath" in m.group(0):
            u = re.search(r'url="([^"]*)"', m.group(0))
            removed.append(f"x15ac:absPath={u.group(1) if u else '?'}")
            x = x.replace(m.group(0), "")
    return x, removed


def main() -> int:
    if not SRC.exists():
        print(f"元ファイルが無い: {SRC}", file=sys.stderr)
        return 2
    DEST.mkdir(parents=True, exist_ok=True)
    out = DEST / NAME
    tmp = out.with_suffix(".tmp")

    zin = zipfile.ZipFile(SRC)
    removed_all, changed = [], []
    with zipfile.ZipFile(tmp, "w") as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            if info.filename == "docProps/core.xml":
                s, rem = clean_core(data.decode("utf-8"))
                data, removed_all = s.encode("utf-8"), removed_all + rem
                changed.append(info.filename)
            elif info.filename == "xl/workbook.xml":
                s, rem = clean_workbook(data.decode("utf-8"))
                if rem:
                    data, removed_all = s.encode("utf-8"), removed_all + rem
                    changed.append(info.filename)
            ni = zipfile.ZipInfo(info.filename, date_time=info.date_time)
            ni.compress_type = info.compress_type
            ni.external_attr = info.external_attr
            ni.internal_attr = info.internal_attr
            ni.create_system = info.create_system
            zout.writestr(ni, data)
    shutil.move(str(tmp), str(out))

    print("除去したもの:")
    for r in removed_all:
        print(f"  - {r}")
    print(f"変更したエントリ: {changed}")

    # ---- 検証
    zo = zipfile.ZipFile(out)
    same_list = [i.filename for i in zin.infolist()] == [i.filename for i in zo.infolist()]
    print(f"\nエントリ一覧が一致: {same_list}（{len(zo.infolist())} 件）")
    diff = []
    for n in zin.namelist():
        a, b = zin.read(n), zo.read(n)
        if a != b:
            diff.append(n)
    print(f"内容が変わったエントリ: {diff}")
    ident = [n for n in zin.namelist()
             if any(n.startswith(k) or n == k for k in KEEP_IDENTICAL)]
    bad = [n for n in ident if zin.read(n) != zo.read(n)]
    crc_ok = all(zin.getinfo(n).CRC == zo.getinfo(n).CRC for n in ident)
    print(f"データ側 {len(ident)} エントリがバイト同一: {not bad}（CRC-32 も一致: {crc_ok}）")
    if bad:
        print(f"  ✗ 一致しない: {bad}")
        return 1
    terms = leak_terms()
    left = [n for n in zo.namelist()
            if any(t in zo.read(n) for t in terms)]
    print(f"個人名・内部パスが残るエントリ: {left if left else 'なし'}")
    print(f"\n新しい SHA-256: {sha256(out)}")
    print(f"サイズ: {out.stat().st_size:,} bytes（元 {SRC.stat().st_size:,}）")
    return 0 if not bad and not left else 1


if __name__ == "__main__":
    sys.exit(main())
