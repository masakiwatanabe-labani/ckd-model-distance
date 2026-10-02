# -*- coding: utf-8 -*-
"""受託解析の元表の置き場所を、公開しない場所から引く（BN1）。

元の表のファイル名は受託先の案件番号を含むので、リポジトリには書かない。
手元では次のどちらかで与える。

  環境変数 CKD_PODTRECK_SOURCE に元表の絶対パス
  projection/results/roundR/local/attachment_source.tsv の "original file name"

公開されているのは添付（Supplementary Data S1）の名前とチェックサムだけで、
その添付は元表と同じ値を持つ。クローンした状態では元表が無いので、これを使う
スクリプトは動かない（結果ファイルは公開されている）。
"""
from __future__ import annotations

import os
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
LOCAL = HERE / "results" / "roundR" / "local" / "attachment_source.tsv"
ATTACHED = "Supplementary_Data_S1_PodTRECK_proteome.xlsx"


def source_path() -> Path:
    """元の定量表。無ければどこに置くかを告げて止まる。"""
    env = os.environ.get("CKD_PODTRECK_SOURCE")
    if env:
        return Path(env)
    if LOCAL.exists():
        for line in LOCAL.read_text().splitlines()[1:]:
            k, _, v = line.partition("\t")
            if k == "original file name" and v.strip():
                return HERE.parent / "data" / "raw" / v.strip()
    raise SystemExit(
        "元の定量表の場所が分からない。環境変数 CKD_PODTRECK_SOURCE にパスを与えるか、"
        f"{LOCAL} を用意すること。添付版（{ATTACHED}）でも同じ値が得られる。")


def sample_prefix() -> str:
    """元表のサンプル列の接頭辞。公開しないので同じ場所から引く。"""
    env = os.environ.get("CKD_PODTRECK_PREFIX")
    if env:
        return env
    if LOCAL.exists():
        for line in LOCAL.read_text().splitlines()[1:]:
            k, _, v = line.partition("\t")
            if k == "original file name" and v.strip():
                return v.split("_")[0]
    raise SystemExit("サンプル列の接頭辞が分からない。CKD_PODTRECK_PREFIX を与えること。")
