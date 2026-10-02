# -*- coding: utf-8 -*-
"""BN4: 公開リポジトリだけで走る検査の件数を測り、結果ファイルに残す。

公開候補（追跡中のファイル＋.gitignore に当たらない未追跡ファイル）だけを空の
ディレクトリに写し、verify_numbers.py を 2 回走らせる。1 回目は素のクローン、
2 回目は `src/00_fetch_refs.py` 相当（data/ref/）を置いた状態。
回答書はこの表から件数を引く。
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROJ = ROOT / "projection"
OUT = PROJ / "results" / "roundR" / "verification_counts.tsv"
PY = sys.executable


def publishable() -> list[str]:
    def run(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True,
                              text=True).stdout.split("\n")
    files = [f for f in run("ls-files") if f]
    for line in run("status", "--porcelain=v1", "-uall"):
        if line.startswith("?? "):
            files.append(line[3:])
    return sorted(set(files))


def counts(where: Path) -> tuple[int, int, int]:
    r = subprocess.run([PY, "verify_numbers.py"], cwd=where / "projection",
                       capture_output=True, text=True)
    m = re.search(r"実行 (\d+) 件 / 一致 (\d+) 件 / 不一致 (\d+) 件 / 飛ばした (\d+) 組",
                  r.stdout)
    if not m:
        print(r.stdout[-1500:], r.stderr[-800:], file=sys.stderr)
        raise SystemExit("件数を読み取れない")
    return int(m.group(1)), int(m.group(3)), int(m.group(4))


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        clone = Path(td) / "clone"
        for rel in publishable():
            src = ROOT / rel
            if not src.is_file():
                continue
            dst = clone / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        n_bare, bad_bare, sk_bare = counts(clone)
        print(f"  素のクローン        : 実行 {n_bare} / 不一致 {bad_bare} / 飛ばした {sk_bare}")
        ref = ROOT / "data" / "ref"
        n_ref = bad_ref = sk_ref = None
        if ref.exists():
            shutil.copytree(ref, clone / "data" / "ref", dirs_exist_ok=True)
            n_ref, bad_ref, sk_ref = counts(clone)
            print(f"  ＋参照データ取得後  : 実行 {n_ref} / 不一致 {bad_ref} / 飛ばした {sk_ref}")

    T = pd.read_csv(OUT, sep="\t") if OUT.exists() else pd.DataFrame(
        columns=["quantity", "value"])
    T = T[~T.quantity.isin(["checks_public_repo", "checks_public_repo_with_refs"])]
    rows = [{"quantity": "checks_public_repo", "value": n_bare}]
    if n_ref is not None:
        rows.append({"quantity": "checks_public_repo_with_refs", "value": n_ref})
    pd.concat([T, pd.DataFrame(rows)], ignore_index=True).to_csv(OUT, sep="\t", index=False)
    print(f"書き出し: {OUT}")
    return 1 if (bad_bare or (bad_ref or 0)) else 0


if __name__ == "__main__":
    sys.exit(main())
