# -*- coding: utf-8 -*-
"""AL2/AL6: 投稿する全ファイルのメタデータと残骸を点検する。

見るもの:
  docProps/core.xml, app.xml, custom.xml の全要素
  内部パス（absPath, file:///, C:\\, Z:\\, /Users/）
  全エントリへの grep（受託先・ローカルのユーザー名・生成ツール名）
  原稿の残骸（[AG:, TODO, XXX, ??, Error! Reference）
"""
from __future__ import annotations

import hashlib
import re
import sys
import zipfile
from pathlib import Path

DEST = Path.home() / "Desktop" / "CKD_xspecies_submission"
# 探す語のうち、第三者の名前など公開したくないものは leak_terms.txt から読む
# （そのファイルは公開しない）。無い場合は一般的なパターンだけで点検する。
_TERMS = Path(__file__).resolve().parents[2] / "leak_terms.txt"
LEAK = [l.strip().encode() for l in _TERMS.read_text().splitlines()
        if l.strip() and not l.startswith("#")] if _TERMS.exists() else []
LEAK += [b"absPath", b"C:\\\\", b"Z:\\\\", b"file:///",
         rb"/Users/", rb"/home/",
         b"Claude", b"Anthropic", b"ChatGPT", b"OpenAI"]
# 謝辞で生成 AI を開示しているので、本文テキストに出る分は除外して数える
TEXT_OK = {b"Claude", b"Anthropic", b"ChatGPT", b"OpenAI"}
STALE = ["[AG:", "TODO", "Error! Reference"]
# 編集部が埋める箇所・原稿 ID のプレースホルダは残骸としない
PLACEHOLDER_OK = ["Firstname Lastname", "Received: date"]
TEMPLATE_OK = ["Firstname Lastname", "Received: date", "Academic Editor"]
CORE = ["dc:creator", "cp:lastModifiedBy", "dc:title", "dc:subject", "dc:description",
        "cp:keywords", "cp:category", "cp:lastPrinted", "cp:contentStatus", "cp:revision",
        "dcterms:created", "dcterms:modified"]
APP = ["Application", "Company", "Manager", "AppVersion", "Template", "TotalTime"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def text_of(z: zipfile.ZipFile) -> str:
    out = []
    for n in z.namelist():
        if n.endswith(("document.xml", "sharedStrings.xml")) or "/worksheets/" in n:
            out.append(re.sub(r"<[^>]+>", " ", z.read(n).decode("utf-8", "replace")))
    return " ".join(out)


def audit(p: Path) -> dict:
    z = zipfile.ZipFile(p)
    rec = {"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p),
           "entries": len(z.namelist()), "core": {}, "app": {}, "custom": [],
           "leaks": {}, "stale": []}
    for part, keys, slot in (("docProps/core.xml", CORE, "core"),
                             ("docProps/app.xml", APP, "app")):
        if part in z.namelist():
            x = z.read(part).decode("utf-8", "replace")
            for k in keys:
                m = re.search(rf"<{k}[^>]*>(.*?)</{k}>", x, re.S)
                if m:
                    rec[slot][k] = m.group(1).strip() or "(空)"
    if "docProps/custom.xml" in z.namelist():
        x = z.read("docProps/custom.xml").decode("utf-8", "replace")
        rec["custom"] = re.findall(r'name="([^"]+)"[^>]*>\s*<[^>]+>(.*?)</', x, re.S)
    body = text_of(z)
    for pat in LEAK:
        hits = sorted({n for n in z.namelist() if pat in z.read(n)})
        if hits:
            if pat in TEXT_OK and all(pat.decode() in body for _ in [0]):
                rec["leaks"][pat.decode()] = hits + ["（本文テキスト。謝辞の開示）"]
            else:
                rec["leaks"][pat.decode()] = hits
    clean = body
    for okp in PLACEHOLDER_OK:
        clean = clean.replace(okp, " ")
    for s in STALE:
        if s in clean:
            rec["stale"].append(s)
    if "XXX" in clean:
        rec["stale"].append("XXX")
    rec["placeholders"] = [okp for okp in PLACEHOLDER_OK if okp in body]
    if re.search(r"\?\?(?!\?)", body):
        rec["stale"].append("??")
    return rec


def main() -> int:
    files = sorted(DEST.glob("*"))
    bad = 0
    for p in files:
        if p.suffix not in (".docx", ".xlsx") or p.name.startswith("~$"):
            continue
        r = audit(p)
        print(f"\n=== {r['file']} ({r['bytes']:,} bytes / {r['entries']} entries) ===")
        print(f"  sha256 {r['sha256']}")
        print("  core :", ", ".join(f"{k}={v}" for k, v in r["core"].items()) or "なし")
        print("  app  :", ", ".join(f"{k}={v}" for k, v in r["app"].items()) or "なし")
        print("  custom:", r["custom"] or "なし")
        if r["leaks"]:
            for k, v in r["leaks"].items():
                mark = " " if "（本文テキスト" in str(v) else "✗"
                print(f"  {mark} 検出 {k}: {v}")
                if mark == "✗":
                    bad += 1
        else:
            print("  内部パス・外部の個人名: なし")
        print("  残骸:", r["stale"] or "なし")
        if r.get("placeholders"):
            print("  想定内のプレースホルダ:", r["placeholders"])
        bad += len(r["stale"])
    print(f"\n要対応 {bad} 件")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
