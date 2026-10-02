"""節（## 見出し）の本文をまるごと差し替えるヘルパー。Part 4 の圧縮に使う。"""
from __future__ import annotations
import pathlib
import re
import textwrap

WIDTH = 98


def _wrap(block: str) -> str:
    out = []
    for part in re.split(r'\n\s*\n', block):
        if not part.strip():
            continue
        if part.lstrip().startswith(('|', '##', '$$', '**Table', '**Figure', '---')):
            out.append(part.strip())
        else:
            out.append('\n'.join(textwrap.wrap(' '.join(part.split()), WIDTH)))
    return '\n\n'.join(out)


def replace_section(path: str, heading: str, body: str) -> int:
    """heading で始まる節の本文を body に差し替え、削減語数を返す。"""
    p = pathlib.Path(path)
    t = p.read_text()
    i = t.index(heading)
    j = t.find('\n## ', i + len(heading))
    if j < 0:
        j = len(t)
    old = t[i + len(heading):j]
    new = '\n\n' + _wrap(body) + '\n\n'
    p.write_text(t[:i + len(heading)] + new + t[j:])
    d = len(old.split()) - len(new.split())
    print(f'  {heading[:52]}: {len(old.split())} → {len(new.split())} 語 ({-d:+d})')
    return d
