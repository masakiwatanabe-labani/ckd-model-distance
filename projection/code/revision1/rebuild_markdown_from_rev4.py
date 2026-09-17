"""Part A. rev4 の docx を正本として manuscript/ の Markdown を作り直す。

古い Markdown に差分を当てる方向はやらない。docx から本文・表・図キャプション・
参考文献をすべて取り出し、節の区切りで分けて書き出す。数値の当て直しは別スクリプト。
"""
from __future__ import annotations

import html
import json
import re
import textwrap
from pathlib import Path

SCR = Path('/private/tmp/claude-501/-Users-masaki/f7fba45f-f8e1-4d4a-b8f2-da5834e0fbe6/scratchpad/docx4')
OUT = Path(__file__).resolve().parents[2] / "manuscript"
WIDTH = 98


def load():
    lines = [l.rstrip() for l in (SCR / 'full.txt').read_text().split('\n')]
    tables = json.loads((SCR / 'tables.json').read_text())
    return lines, tables


def wrap(p: str) -> str:
    return '\n'.join(textwrap.wrap(p.strip(), WIDTH)) if p.strip() else ''


def md_table(rows) -> str:
    out = ['| ' + ' | '.join(c.replace('|', r'\|') for c in rows[0]) + ' |',
           '|' + '|'.join(['---'] * len(rows[0])) + '|']
    for r in rows[1:]:
        r = list(r) + [''] * (len(rows[0]) - len(r))
        out.append('| ' + ' | '.join(c.replace('|', r'\|') for c in r) + ' |')
    return '\n'.join(out)


def section(lines, a, b):
    """行 a（見出し）から b の直前まで。見出し行は除く。"""
    return [l for l in lines[a + 1:b] if l.strip()]


def is_heading(s: str) -> bool:
    return bool(re.match(r'^\d+\.\d+\.\s', s.strip()))


def _match_run(paras, i, cells):
    """paras[i:] が cells の並びと一致するか調べ、一致した長さを返す（0 なら不一致）。"""
    j, k = i, 0
    while k < len(cells):
        if j >= len(paras):
            return 0
        if not paras[j].strip():
            j += 1
            continue
        if ' '.join(paras[j].split()) != cells[k]:
            return 0
        j += 1
        k += 1
    return j - i


def drop_flat_tables(paras, tables):
    """docx から平文化して流れ込む表セルの連なりを、pipe 表に置き換える。

    rev4 のテキスト抽出では表のセルが 1 セル 1 段落として本文に混ざる。これを
    そのまま残すと、更新前の値を持つ平文コピーが pipe 表の隣に残り続ける。
    """
    flats = []
    for tb in tables:
        cells = [' '.join(c.split()) for r in tb for c in r if c.strip()]
        if len(cells) >= 4:
            flats.append((cells, tb))
    out, i, dropped = [], 0, 0
    while i < len(paras):
        hit = None
        for cells, tb in flats:
            n = _match_run(paras, i, cells)
            if n:
                hit = (n, tb)
                break
        if hit:
            n, tb = hit
            out.append(md_table(tb))
            dropped += n
            i += n
        else:
            out.append(paras[i])
            i += 1
    return out, dropped


def emit(paras, tables=(), table_for=None):
    """段落リストを Markdown に。見出しは ## に、表キャプションの直後に表を置く。"""
    paras, dropped = drop_flat_tables(list(paras), list(tables))
    out = []
    for p in paras:
        s = p.strip()
        if s.startswith('|'):
            out.append(s)
        elif is_heading(s):
            out.append('## ' + s)
        elif re.match(r'^Table \d+\.', s):
            num = int(re.match(r'^Table (\d+)\.', s).group(1))
            body = re.sub(r'^Table \d+\.\s*', '', s)
            out.append(f'**Table {num}.** *{body}*')
        else:
            out.append(wrap(s))
    emit.dropped = dropped
    return '\n\n'.join(x for x in out if x) + '\n'


def main():
    lines, tables = load()
    m = json.loads((SCR / 'marks.json').read_text())

    # ---- FRONTMATTER ----
    abstract = re.sub(r'^Abstract:\s*', '', lines[m['abstract']].strip())
    keywords = re.sub(r'^Keywords:\s*', '', lines[m['keywords']].strip())
    affil = [l.strip() for l in lines[10:m['abstract']] if l.strip()]
    fm = ['# Front matter', '',
          '**Article type.** ' + lines[0].strip(), '',
          '**Title.** ' + wrap(lines[1].strip()), '',
          '**Authors.** ' + lines[2].strip(), '',
          '**Affiliation and corresponding author.**', '']
    fm += [wrap(a) for a in affil] + ['']
    fm += ['## Abstract', '', wrap(abstract), '', '**Keywords:** ' + wrap(keywords), '']
    (OUT / 'FRONTMATTER.md').write_text('\n'.join(fm))

    # ---- INTRODUCTION ----
    intro = section(lines, m['intro'], m['results'])
    (OUT / 'INTRODUCTION.md').write_text('# 1. Introduction\n\n' + emit(intro))

    # ---- RESULTS (+ figure legends) ----
    res = section(lines, m['results'], m['disc'])
    body = emit(res, tables[1:])
    figs = section(lines, m['figleg'], m['supptab'])
    caps = '\n\n'.join(
        re.sub(r'^Figure (\d+)\.\s*', lambda mm: f'**Figure {mm.group(1)}.** ', wrap(f))
        for f in figs)
    (OUT / 'RESULTS.md').write_text('# 2. Results\n\n' + body +
                                    '\n---\n\n### Figures and tables\n\n' + caps + '\n')

    # ---- DISCUSSION ----
    disc = section(lines, m['disc'], m['methods'])
    (OUT / 'DISCUSSION.md').write_text('# 3. Discussion\n\n' + emit(disc, tables[1:]))

    # ---- METHODS ----
    met = section(lines, m['methods'], m['concl'])
    txt = emit(met, tables[1:])
    (OUT / 'METHODS.md').write_text('# 4. Materials and Methods\n\n' + txt)

    # ---- CONCLUSIONS ----
    con = section(lines, m['concl'], m['suppmat'])
    (OUT / 'CONCLUSIONS.md').write_text('# 5. Conclusions\n\n' + emit(con))

    # ---- BACK MATTER ----
    back = [l.strip() for l in lines[m['suppmat']:m['refs']] if l.strip()]
    bm = ['# Back matter', '']
    for b in back:
        head, _, rest = b.partition(':')
        bm += ['## ' + head.strip(), '', wrap(rest.strip()), '']
    (OUT / 'BACKMATTER.md').write_text('\n'.join(bm))

    # ---- REFERENCES ----
    refs = [l.strip() for l in lines[m['refs'] + 1:m['figleg']] if l.strip()]
    rf = ['# References', '']
    rf += [f'{i}. {wrap(r)}' for i, r in enumerate(refs, 1)]
    (OUT / 'REFERENCES.md').write_text('\n'.join(rf) + '\n')

    # ---- SUPPLEMENTARY ----
    sup = [l.strip() for l in lines[m['supptab'] + 1:m['disc_note']] if l.strip()]
    sp = ['# Supplementary Materials', '']
    for s in sup:
        sp.append(re.sub(r'^Table (S\d+)\.\s*', lambda mm: f'**Table {mm.group(1)}.** ', wrap(s)))
        sp.append('')
    (OUT / 'SUPPLEMENTARY.md').write_text('\n'.join(sp))

    print('rebuilt:', sorted(p.name for p in OUT.glob('*.md')))
    print('references:', len(refs))


if __name__ == "__main__":
    main()
