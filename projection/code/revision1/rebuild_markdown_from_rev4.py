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


def emit(paras, table_for=None):
    """段落リストを Markdown に。見出しは ## に、表キャプションの直後に表を挿入。"""
    out = []
    for p in paras:
        s = p.strip()
        if is_heading(s):
            out.append('## ' + s)
        elif re.match(r'^Table \d+\.', s):
            num = int(re.match(r'^Table (\d+)\.', s).group(1))
            body = re.sub(r'^Table \d+\.\s*', '', s)
            out.append(f'**Table {num}.** *{body}*')
            if table_for and num in table_for:
                out.append(md_table(table_for[num]))
        else:
            out.append(wrap(s))
    return '\n\n'.join(x for x in out if x) + '\n'


def main():
    lines, tables = load()
    m = json.loads((SCR / 'marks.json').read_text())
    tab = {1: tables[1], 2: tables[2], 3: tables[4], 4: tables[5]}
    methods_table = tables[6]

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
    body = emit(res, tab)
    figs = section(lines, m['figleg'], m['supptab'])
    caps = '\n\n'.join(
        re.sub(r'^Figure (\d+)\.\s*', lambda mm: f'**Figure {mm.group(1)}.** ', wrap(f))
        for f in figs)
    (OUT / 'RESULTS.md').write_text('# 2. Results\n\n' + body +
                                    '\n---\n\n### Figures and tables\n\n' + caps + '\n')

    # ---- DISCUSSION ----
    disc = section(lines, m['disc'], m['methods'])
    (OUT / 'DISCUSSION.md').write_text('# 3. Discussion\n\n' + emit(disc))

    # ---- METHODS ----
    met = section(lines, m['methods'], m['concl'])
    txt = emit(met)
    # 4.5 の検体数表を差し込む
    anchor = '## 4.6. Definition of Group A'
    txt = txt.replace(anchor, md_table(methods_table) + '\n\n' + anchor)
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
