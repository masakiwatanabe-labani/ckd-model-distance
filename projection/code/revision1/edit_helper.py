"""段落単位の置換ヘルパー。

docx 由来のカーリー引用符とダッシュ（- / – / —）の違いを吸収して照合する。
置換は元の段落の該当範囲だけを対象にするので、段落内の他の文字は変わらない。
"""
from __future__ import annotations

import pathlib
import re
import textwrap

WIDTH = 98
QUOTES = {'’': "'", '‘': "'", '“': '"', '”': '"'}
DASHES = '-–—−'   # ハイフン / en / em / U+2212 マイナス


def norm(s: str) -> str:
    for a, b in QUOTES.items():
        s = s.replace(a, b)
    for d in DASHES:
        s = s.replace(d, '-')
    return ' '.join(s.split())


def wrap(s: str) -> str:
    """段落内で折り返す。空行は段落の区切りとして保つ。

    置換文字列に \n\n を含めて段落を分けられるようにするため、空行で分割してから
    それぞれを折り返す。表・見出し・数式ブロックはそのまま通す。
    """
    out = []
    for part in re.split(r'\n\s*\n', s):
        if not part.strip():
            continue
        if part.lstrip().startswith(('|', '##', '$$')):   # 表・見出し・数式はそのまま
            out.append(part.strip())
        else:
            out.append('\n'.join(textwrap.wrap(' '.join(part.split()), WIDTH)))
    return '\n\n'.join(out)


def _pattern(old: str) -> re.Pattern:
    parts = []
    for ch in old:
        if ch.isspace():
            parts.append(r'\s+')
        elif ch in DASHES:
            parts.append('[' + DASHES + r']\s*')   # 折り返しがハイフン直後に入ることがある
        elif ch in ("'", '’', '‘'):
            parts.append("['’‘]")
        elif ch in ('"', '“', '”'):
            parts.append('["“”]')
        else:
            parts.append(re.escape(ch))
    return re.compile(''.join(parts).replace(r'\s+\s+', r'\s+'))


def rp(path: str, old: str, new: str) -> None:
    """old を含む段落をちょうど1つ見つけ、その中の old だけを new に置き換える。"""
    p = pathlib.Path(path)
    paras = p.read_text().split('\n\n')
    pat = _pattern(' '.join(old.split()))
    hits = [i for i, x in enumerate(paras) if pat.search(x)]
    if len(hits) != 1:
        raise AssertionError(f'{path}: {len(hits)} hits for {old[:70]!r}')
    i = hits[0]
    is_table = paras[i].lstrip().startswith('|')
    body = pat.sub(lambda _m: new.replace('\\', r'\\'), paras[i], count=1)
    paras[i] = body if is_table else wrap(body)
    p.write_text('\n\n'.join(paras))
    print(f'  {path}: {old[:52]!r}')


def insert_before(path: str, anchor: str, block: str) -> None:
    """anchor を含む段落の直前に block（空行区切りの段落列）を差し込む。"""
    p = pathlib.Path(path)
    paras = p.read_text().split('\n\n')
    pat = _pattern(' '.join(anchor.split()))
    hits = [i for i, x in enumerate(paras) if pat.search(x)]
    if len(hits) != 1:
        raise AssertionError(f'{path}: {len(hits)} hits for anchor {anchor[:70]!r}')
    new = [b if b.lstrip().startswith(('|', '**Table', '**Figure', '##', '$$')) else wrap(b)
           for b in block.strip().split('\n\n')]
    paras[hits[0]:hits[0]] = new
    p.write_text('\n\n'.join(paras))
    print(f'  {path}: +{len(new)} para before {anchor[:40]!r}')


def append_section(path: str, block: str) -> None:
    p = pathlib.Path(path)
    p.write_text(p.read_text().rstrip() + '\n\n' + block.strip() + '\n')
    print(f'  {path}: appended')
