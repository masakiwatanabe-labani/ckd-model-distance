# -*- coding: utf-8 -*-
"""results/ の全ファイルについて、それを書き出すコードがツリーにあるかを判定する。

第10ラウンドで、読むコードはあるのに書くコードが無いファイルが 3 つ見つかった。
同じ欠落が他にないかを系統的に調べ、results/provenance.tsv に残す。

判定は静的な走査で行う。スクリプトから「書き出しに使われる文字列リテラル」を集め、
  1. ファイル名そのものがリテラルとして現れるか
  2. f 文字列（"auc_permutation{tag}.tsv" のような）のパターンに当たるか
で照合する。判定できたものは producer 列にスクリプト名が入り、できなかったものは
status が "no generating script found" になる。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
RES = HERE / "results"
# R 側の書き出し（write.table / writeLines / sink）も拾う
WRITE_HINT = re.compile(r"to_csv|write_text|savefig|to_markdown|open\(|write\.table|writeLines|sink\(")
ASSIGN = re.compile(r'(\w+)\s*=\s*.*["\'][^"\']*\.(?:tsv|csv|md|log|png|pdf|txt|gmt|json|xml)["\']')
LIT = re.compile(r'["\']([^"\']*\.(?:tsv|csv|md|log|png|pdf|txt|gmt|json|xml))["\']')


def script_files():
    out = sorted(HERE.glob("*.py")) + sorted((HERE / "code" / "revision1").glob("*.py"))
    out += sorted((HERE.parent / "src").glob("*.py"))
    # R のスクリプトも結果を書くので、生成元の候補に入れる
    out += sorted((HERE / "code" / "revision1").glob("*.R"))
    return [p for p in out if not p.name.startswith("apply_round")]


def index_scripts():
    """(完全一致の辞書, f 文字列パターンの一覧) を返す。"""
    exact: dict[str, set[str]] = {}
    fuzzy: list[tuple[re.Pattern, str]] = []
    for p in script_files():
        src = p.read_text(errors="replace")
        if not WRITE_HINT.search(src):
            continue
        rel = str(p.relative_to(HERE.parent))
        lines = src.splitlines()
        # 書き出し行にあるリテラルと、書き出しに渡される変数に入るリテラルだけを拾う。
        # 読み取り行（pd.read_csv(... "x.tsv")）はここで落ちる。
        write_lines = [i for i, ln in enumerate(lines) if WRITE_HINT.search(ln)]
        keep = []
        for i in write_lines:                    # 引数が次行に折り返される書き方に備える
            keep.extend(lines[i:i + 4])
        write_text = "\n".join("\n".join(lines[i:i + 4]) for i in write_lines)
        for ln in lines:
            m = ASSIGN.match(ln.strip())
            if m and re.search(rf"{re.escape(m.group(1))}\b", write_text):
                keep.append(ln)
        for m in LIT.finditer("\n".join(keep)):
            lit = m.group(1)
            if "{" in lit:                       # f 文字列。{...} を緩いパターンに開く
                pat = re.escape(lit)
                pat = re.sub(r"\\\{[^}]*\\\}", r"[A-Za-z0-9_\\-]*", pat)
                fuzzy.append((re.compile("^" + pat + "$"), rel))
            else:
                exact.setdefault(Path(lit).name, set()).add(rel)
    return exact, fuzzy


def main() -> int:
    exact, fuzzy = index_scripts()
    rows = []
    for f in sorted(RES.rglob("*")):
        # 公開しない控え（原稿の文の正本など）は由来表に載せない
        if "/local/" in str(f):
            continue
        if not f.is_file():
            continue
        rel = str(f.relative_to(HERE))
        name = f.name
        who = sorted(exact.get(name, set()))
        how = "exact filename literal" if who else ""
        if not who:
            who = sorted({s for pat, s in fuzzy if pat.match(name)})
            how = "f-string pattern" if who else ""
        # 変更履歴版の組み立てに使う中間 XML は成果物ではない
        if "tracked_build" in f.parts:
            continue
        if not who and f.suffix == ".log":
            how, who = "run log (shell redirection)", ["—"]
        if not who and rel.startswith("results/intermediate/"):
            how, who = "intermediate; see results/intermediate/README.md", ["—"]
        rows.append({"result_file": rel, "bytes": f.stat().st_size,
                     "producer": " / ".join(who) if who else "",
                     "matched_by": how,
                     "status": ("ok" if who and not how.startswith("intermediate")
                                else "intermediate, not used by the article" if who
                                else "no generating script found")})
    T = pd.DataFrame(rows).sort_values("result_file")
    out = RES / "provenance.tsv"
    T.to_csv(out, sep="\t", index=False)
    miss = T[T.status == "no generating script found"]
    inter = T[T.status.str.startswith("intermediate")]
    print(f"{len(T)} ファイル / 生成スクリプト不明 {len(miss)} 件 / "
          f"中間ファイル {len(inter)} 件 → {out}")
    for r in miss.itertuples():
        print(f"  ✗ {r.result_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
