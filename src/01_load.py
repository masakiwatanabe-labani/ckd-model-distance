"""生データを読み込み、オルソログ表で ID を統一して parquet に保存する。

重要:
  - ネコのプロテオーム値は蛋白ごとに中心化されている（row mean≈0）。
    ここで検出し、`attrs['centered']=True` を立てて下流に伝える。
  - シンボル大文字一致ではなく、Ensembl オルソログ表で cat→mouse を対応させる。
    対応が取れなかった遺伝子、LOC 遺伝子の割合を QC として出す。
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger, append_summary  # noqa: E402

log = get_logger("01_load")
cfg = load_config()
ROOT = Path(cfg["_root"])
RAW = ROOT / cfg["paths"]["raw"]
REF = ROOT / cfg["paths"]["ref"]
INT = ROOT / cfg["paths"]["interim"]
INT.mkdir(parents=True, exist_ok=True)


def _rows(fname: str, sheet: str):
    wb = openpyxl.load_workbook(RAW / fname, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return rows


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


# ---------------------------------------------------------------- cat
def load_cat_rna(sheet: str):
    rows = _rows(cfg["files"]["cat_omics"], sheet)
    ids, grp = rows[1][1:], rows[2][1:]
    keep = [i for i, x in enumerate(ids) if x]
    ids = [ids[i] for i in keep]
    grp = [grp[i] for i in keep]
    data, index = [], []
    for r in rows[6:]:
        if not r[0]:
            continue
        index.append(str(r[0]).strip())
        vals = r[1:]
        data.append([_num(vals[i]) for i in keep])
    df = pd.DataFrame(data, index=index, columns=ids)
    df = df.groupby(level=0).mean()
    return df, pd.Series(grp, index=ids, name="group")


def load_cat_prot(sheet: str):
    rows = _rows(cfg["files"]["cat_omics"], sheet)
    ids, grp = rows[1][3:], rows[2][3:]
    keep = [i for i, x in enumerate(ids) if x]
    ids = [ids[i] for i in keep]
    grp = [grp[i] for i in keep]
    data, index = [], []
    n_split = 0
    for r in rows[6:]:
        if not r[1]:
            continue
        raw = str(r[1]).strip()
        # プロテオーム側は "A2M PZP" "AARS1 AARS" のように、プロテイングループの
        # 複数シンボルが空白区切りで1セルに入っていることがある（皮質2.8%/髄質2.7%）。
        # この形のままではオルソログ表にも UniProt にもヒットせず、静かに脱落する。
        # 先頭シンボル（＝マスター遺伝子）を代表として採用する。
        sym = raw.split()[0] if raw.split() else raw
        if sym != raw:
            n_split += 1
        index.append(sym)
        vals = r[3:]
        data.append([_num(vals[i]) for i in keep])
    n_rows = len(index)
    df = pd.DataFrame(data, index=index, columns=ids).groupby(level=0).mean()
    log.info("%s: 複合シンボルを分割 %d 行 / 全%d行、統合後 %d 蛋白（重複統合 %d 行）",
             sheet, n_split, n_rows, len(df), n_rows - len(df))
    rowmean = df.mean(axis=1)
    centered = bool(np.nanmedian(np.abs(rowmean)) < 0.1)
    df.attrs["centered"] = centered
    if centered:
        log.warning("%s: 行中心化された値です。log2FCの絶対値は種間比較に使えません", sheet)
    return df, pd.Series(grp, index=ids, name="group")


# ---------------------------------------------------------------- mouse
def load_mouse_rna():
    rows = _rows(cfg["files"]["mouse_rna"], "Sheet1")
    samples = rows[0][7:]
    acc = {}
    for r in rows[1:]:
        if not r[1]:
            continue
        g = str(r[1]).strip()
        vals = [_num(x) for x in r[7:]]
        acc[g] = list(np.nansum([acc[g], vals], axis=0)) if g in acc else vals
    df = pd.DataFrame.from_dict(acc, orient="index", columns=samples)
    groups = pd.Series({c: c.rsplit("-", 1)[0] for c in df.columns}, name="group")
    return df, groups


def load_mouse_prot():
    rows = _rows(cfg["files"]["mouse_prot"], "解析結果")
    samples = rows[0][10:19]
    acc = {}
    for r in rows[1:]:
        if not r[5]:
            continue
        acc.setdefault(str(r[5]).strip(), []).append([_num(x) for x in r[10:19]])
    df = pd.DataFrame.from_dict(
        {g: list(np.nanmax(np.array(v), axis=0)) for g, v in acc.items()},
        orient="index", columns=samples)
    df = df.replace(0, np.nan)
    logdf = np.log2(df)
    logdf = logdf - logdf.median() + logdf.median().mean()   # 列中央値正規化
    groups = pd.Series({c: c.split("_")[1].rsplit("-", 1)[0] for c in logdf.columns}, name="group")
    return logdf, groups


# ---------------------------------------------------------------- orthology
def mouse_rna_ensembl_symbol_map() -> dict:
    """マウス RNA-seq raw ファイル（GeneID, GeneName 列）から

        Ensembl gene ID -> このファイルで実際に使われているシンボル

    の対応を作る。BioMart の現行シンボル（external_gene_name）は、
    旧いゲノムアノテーションでアラインされた raw データのシンボルと
    食い違うことがある（例: ENSMUSG00000021749 は BioMart 現行では
    "Fam3d" だが、この FPKM ファイルでは旧称 "Oit1" のまま）。
    Ensembl ID を橋渡しにすれば、この種の命名規則の更新ズレを
    Ensembl ID ベースで正しく解決できる。
    """
    rows = _rows(cfg["files"]["mouse_rna"], "Sheet1")
    m = {}
    for r in rows[1:]:
        if not r[0] or not r[1]:
            continue
        gid = str(r[0]).strip().split(".")[0]
        m.setdefault(gid, str(r[1]).strip())
    return m


def mouse_synonym_map() -> dict:
    """Ensembl mouse gene ID -> {現行シンボル, 別名...} の集合。"""
    path = REF / "synonyms_mouse.tsv"
    if not path.exists():
        log.warning("synonyms_mouse.tsv がありません。00_fetch_refs.py を先に実行してください。"
                    "synonym 経由のシンボル解決はスキップされます。")
        return {}
    df = pd.read_csv(path, sep="\t")
    m: dict[str, set] = {}
    for gid, name, syn in df.itertuples(index=False):
        s = m.setdefault(str(gid), set())
        if isinstance(name, str):
            s.add(name)
        if isinstance(syn, str):
            s.add(syn)
    return m


def build_ortholog_map(mouse_rna_index) -> pd.DataFrame:
    """cat_symbol -> mouse_symbol の一意な対応表を作る。

    キーは tgt_symbol（BioMart の現行シンボル文字列）ではなく
    tgt_ensembl（Ensembl gene ID）。理由: BioMart の現行シンボルは raw データの
    アノテーション時点のシンボルと食い違うことがあり（FAM3D→Fam3d(BioMart) vs
    Oit1(raw FPKM)）、シンボル文字列だけで結合すると「シンボル大文字一致で
    見落とす」のと同じ失敗を、名称が更新された形で再生産してしまう。

    解決の優先順位:
      1. mouse_rna raw データの GeneID(Ensembl) -> 実際のシンボル表で
         Ensembl ID が直接ヒットすればそれを採用（最も確実）。
      2. 1で見つからない場合、BioMart の synonym 表からその Ensembl ID の
         現行名+別名候補を取り、mouse_rna.index に実在するものがあれば採用。
      3. どちらも失敗したら BioMart の現行シンボルをそのまま使うが、
         「未解決」としてカウント・報告する（qc_report 経由）。
    """
    path = REF / "orthologs_cat2mouse.tsv"
    if not path.exists():
        log.error("オルソログ表がありません。00_fetch_refs.py を先に実行してください。"
                  "フォールバックとして大文字一致を使いますが、FAM3D/Oit1 のような対応は失われます。")
        return pd.DataFrame(columns=["cat_symbol", "mouse_symbol", "orthology_type", "resolved_by"])
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"src_symbol": "cat_symbol", "tgt_symbol": "mouse_symbol_biomart",
                             "tgt_ensembl": "mouse_ensembl"})
    df = df.dropna(subset=["cat_symbol", "mouse_symbol_biomart", "mouse_ensembl"])
    keep = cfg["orthology"]["keep_types"]
    o2o = df[df.orthology_type.isin(keep)].copy()

    # 同じ cat_symbol が、猫ゲノム上の複数の異なる Ensembl 遺伝子座（重複領域/
    # アノテーションの重複）から one2one 判定されるケースが 453 シンボルで発生する
    # （felCat のアノテーション精度の問題、CLAUDE.md の既知の落とし穴と同根）。
    # 例: FAM3D は真のオルソログ Fam3d(=Oit1, perc_id 50.3) の他に、無関係な
    # Pmch という遺伝子への one2one 判定（perc_id 58.7）を持つ。perc_id だけで
    # 選ぶと偽陽性の方が勝ってしまうため、まず「シンボル名の大文字小文字無視の
    # 一致」を最優先のタイブレークにする（種を跨いだ命名委員会がオルソログの
    # 名前を極力揃えている、という一致性を独立した証拠として使う）。
    # それでも複数残る場合のみ confidence → perc_id で決める。
    o2o["name_match"] = (o2o["cat_symbol"].str.upper() == o2o["mouse_symbol_biomart"].str.upper())
    o2o["collision"] = o2o.groupby("cat_symbol")["src_ensembl"].transform("nunique") > 1
    o2o["rank"] = ((~o2o["name_match"]).astype(int) * 1_000_000
                   - o2o.get("confidence", 0).fillna(0) * 1000
                   - o2o.get("perc_id", 0).fillna(0))
    n_collision_symbols = int(o2o.loc[o2o.collision, "cat_symbol"].nunique())
    o2o = o2o.sort_values("rank").drop_duplicates("cat_symbol")

    id_map = mouse_rna_ensembl_symbol_map()
    syn_map = mouse_synonym_map()
    rna_index = set(mouse_rna_index)

    def resolve(row):
        gid = str(row["mouse_ensembl"]).split(".")[0]
        if gid in id_map:
            return id_map[gid], "ensembl_id"
        candidates = syn_map.get(gid, set())
        hit = sorted(candidates & rna_index)
        if hit:
            return hit[0], "synonym"
        return row["mouse_symbol_biomart"], "unresolved"

    resolved = o2o.apply(resolve, axis=1, result_type="expand")
    o2o["mouse_symbol"], o2o["resolved_by"] = resolved[0], resolved[1]
    o2o["unresolved"] = ~o2o["mouse_symbol"].isin(rna_index)

    # 単一の「正しいシンボル」は存在しない。同じ遺伝子でもデータソースごとに
    # 採用しているシンボルの世代が違う（FAM3D の例: RNA-seq FPKM は旧称 Oit1、
    # マウスプロテオームと UniProt は現行名 Fam3d）。どれか1つに正規化すると
    # 必ず他方の突合が静かに落ちるので、既知の別名を全部持ち回る。
    # 下流は lib_stats.alias_lookup() でこの列を使って引く。
    def aliases(row):
        gid = str(row["mouse_ensembl"]).split(".")[0]
        names = set(syn_map.get(gid, set()))
        names.add(row["mouse_symbol"])
        names.add(row["mouse_symbol_biomart"])
        if gid in id_map:
            names.add(id_map[gid])
        return "|".join(sorted(str(n) for n in names if isinstance(n, str) and n))

    o2o["mouse_aliases"] = o2o.apply(aliases, axis=1)

    n_unresolved = int(o2o["unresolved"].sum())
    n_via_synonym = int((o2o["resolved_by"] == "synonym").sum())
    log.info("orthologs: %d cat genes mapped (%s)", len(o2o), keep)
    log.info("シンボル衝突（同一cat_symbolが複数の異なるEnsembl遺伝子座から"
              "one2one判定）: %d シンボル。name_match優先タイブレークで解消。",
              n_collision_symbols)
    log.info("シンボル解決: ensembl_id直接 %d 件 / synonym経由 %d 件 / 未解決 %d 件",
              int((o2o.resolved_by == "ensembl_id").sum()), n_via_synonym, n_unresolved)
    if n_unresolved:
        log.warning("未解決 %d 件（先頭の例）: %s",
                    n_unresolved,
                    o2o.loc[o2o["unresolved"], "cat_symbol"].head(10).tolist())
    o2o.attrs["n_collision_symbols"] = n_collision_symbols
    return o2o[["cat_symbol", "mouse_symbol", "orthology_type", "resolved_by", "mouse_ensembl",
                "mouse_symbol_biomart", "mouse_aliases", "name_match", "collision"]]


def qc_report(cat_rna, mouse_rna, omap):
    cat_genes = set(cat_rna.index)
    loc_frac = np.mean([g.startswith("LOC") for g in cat_genes])
    mapped = set(omap.cat_symbol) & cat_genes
    naive = {g for g in cat_genes if g.upper() in {m.upper() for m in mouse_rna.index}}
    gained = mapped - naive
    lost = naive - mapped
    payload = {
        "ネコ遺伝子数": len(cat_genes),
        "LOC遺伝子の割合": f"{loc_frac:.1%}",
        "オルソログで対応が取れた数": len(mapped),
        "大文字一致のみで対応した数": len(naive),
        "オルソログ表で新たに拾えた数": len(gained),
        "大文字一致では拾えたがオルソログ表にない数": len(lost),
        "新たに拾えた例": sorted(gained)[:20],
    }
    for probe in ["FAM3D", "C3", "LCN2", "SLC14A2", "NR3C2"]:
        rows = omap.loc[omap.cat_symbol == probe, ["mouse_symbol", "resolved_by"]]
        payload[f"probe:{probe}"] = rows.to_dict("records") or "未対応"

    mouse_rna_index = set(mouse_rna.index)
    unresolved = omap[~omap.mouse_symbol.isin(mouse_rna_index)] if "resolved_by" in omap else omap.iloc[0:0]
    n_by_ensembl = int((omap.get("resolved_by") == "ensembl_id").sum()) if "resolved_by" in omap else None
    n_by_synonym = int((omap.get("resolved_by") == "synonym").sum()) if "resolved_by" in omap else None
    payload["シンボル解決: ensembl_id直接"] = n_by_ensembl
    payload["シンボル解決: synonym経由"] = n_by_synonym
    payload["シンボル衝突（複数のcat遺伝子座が同名）で採用した件数"] = (
        int(omap.get("collision", pd.Series(dtype=bool)).sum())
    )
    payload["シンボル未解決（mouse_rnaに実在しない）件数"] = len(unresolved)
    payload["シンボル未解決の例（最大100件）"] = (
        unresolved[["cat_symbol", "mouse_ensembl", "mouse_symbol_biomart"]].head(100).to_dict("records")
    )
    append_summary("01_load / ID対応QC", payload, cfg)
    log.info("QC: %s", payload)


if __name__ == "__main__":
    store = {}
    store["cat_rna_ctx"], store["cat_rna_ctx_grp"] = load_cat_rna(cfg["cat_sheets"]["rna_cortex"])
    store["cat_rna_med"], store["cat_rna_med_grp"] = load_cat_rna(cfg["cat_sheets"]["rna_medulla"])
    store["cat_prot_ctx"], store["cat_prot_ctx_grp"] = load_cat_prot(cfg["cat_sheets"]["prot_cortex"])
    store["cat_prot_med"], store["cat_prot_med_grp"] = load_cat_prot(cfg["cat_sheets"]["prot_medulla"])
    store["mouse_rna"], store["mouse_rna_grp"] = load_mouse_rna()
    store["mouse_prot"], store["mouse_prot_grp"] = load_mouse_prot()

    omap = build_ortholog_map(store["mouse_rna"].index)
    omap.to_csv(INT / "ortholog_map.tsv", sep="\t", index=False)

    # 皮質と髄質で共有される個体（対応ありデザイン用）
    shared = sorted(set(store["cat_rna_ctx"].columns) & set(store["cat_rna_med"].columns))
    pd.Series(shared, name="cat_id").to_csv(INT / "paired_cats.tsv", sep="\t", index=False)
    log.info("皮質・髄質で共通の個体: %d 頭", len(shared))

    for k, v in store.items():
        (v.to_frame() if isinstance(v, pd.Series) else v).to_parquet(INT / f"{k}.parquet")
        log.info("%s: %s", k, getattr(v, "shape", None))

    qc_report(store["cat_rna_ctx"], store["mouse_rna"], omap)
    append_summary("01_load / データ形状", {
        "ネコ皮質RNA": str(store["cat_rna_ctx"].shape),
        "ネコ髄質RNA": str(store["cat_rna_med"].shape),
        "ネコ皮質プロテオーム": str(store["cat_prot_ctx"].shape),
        "ネコプロテオームは行中心化": store["cat_prot_ctx"].attrs.get("centered"),
        "マウスRNA": str(store["mouse_rna"].shape),
        "マウスプロテオーム": str(store["mouse_prot"].shape),
        "対応のあるネコ個体数": len(shared),
    }, cfg)
