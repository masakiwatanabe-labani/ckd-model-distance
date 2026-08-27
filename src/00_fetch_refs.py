"""公的DBからリファレンスを取得する。

1. Ensembl BioMart: ネコ⇄マウス、ネコ⇄ヒト、マウス⇄ヒトのオルソログ表
   → シンボル大文字一致では FAM3D/Oit1、C5/Hc のようなペアを取り逃す。
2. UniProt: マウスとヒトの「シグナルペプチド保有／分泌／細胞外」注釈
   → 尿バイオマーカー候補の分泌型フィルタに使う。
3. MSigDB (gseapy 経由): Hallmark, GO BP, KEGG
   → 手作りモジュールとは別に、独立した遺伝子セットで再現性を確認する。

ネットワークが使えない環境では data/ref/ に既存ファイルがあればそれを使う。
"""
from __future__ import annotations
import io
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger  # noqa: E402

log = get_logger("00_fetch_refs")
cfg = load_config()
REF = Path(cfg["_root"]) / cfg["paths"]["ref"]
REF.mkdir(parents=True, exist_ok=True)

# www.ensembl.org の martservice は現在 POST を 405 で拒否する
# （GET のみ許可、`allow: GET` ヘッダーで確認済み）。
# asia.ensembl.org ミラーは同じクエリで POST を受理する（2026-08-26 時点で確認済み）。
# useast.ensembl.org は接続はできるが応答が極端に遅く（60秒超えても応答なしを確認）、
# 全体をブロックするだけなので除外。
# また、requests のデフォルト User-Agent (python-requests/x.y) だと asia ミラーの WAF に
# 403 Forbidden で弾かれることを確認済み。ブラウザ相当の UA を明示する。
BIOMART_MIRRORS = [
    "https://www.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice",
]
BIOMART_HEADERS = {
    "User-Agent": ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                   "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"),
}


def _biomart_post(query: str) -> str:
    """共通のミラーフォールバック付き BioMart POST。"""
    last_exc = None
    for mirror in BIOMART_MIRRORS:
        try:
            r = requests.post(mirror, data={"query": query}, headers=BIOMART_HEADERS, timeout=90)
            r.raise_for_status()
            if mirror != BIOMART_MIRRORS[0]:
                log.info("  mirror %s で成功", mirror)
            return r.text
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            log.warning("  mirror %s 失敗: %s", mirror, exc)
    raise last_exc


def biomart_orthologs(src: str, tgt: str) -> pd.DataFrame:
    """src 種の遺伝子と tgt 種オルソログの対応表を BioMart から取得。"""
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
  <Dataset name="{src}_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="{tgt}_homolog_ensembl_gene"/>
    <Attribute name="{tgt}_homolog_associated_gene_name"/>
    <Attribute name="{tgt}_homolog_orthology_type"/>
    <Attribute name="{tgt}_homolog_orthology_confidence"/>
    <Attribute name="{tgt}_homolog_perc_id"/>
  </Dataset>
</Query>"""
    log.info("BioMart %s -> %s", src, tgt)
    text = _biomart_post(query)
    df = pd.read_csv(io.StringIO(text), sep="\t")
    df.columns = ["src_ensembl", "src_symbol", "tgt_ensembl", "tgt_symbol",
                  "orthology_type", "confidence", "perc_id"]
    df = df.dropna(subset=["tgt_symbol"])
    return df


def biomart_gene_synonyms(species: str) -> pd.DataFrame:
    """species の ensembl_gene_id / external_gene_name / external_synonym 対応表。

    BioMart の現行シンボル（external_gene_name）が、旧アノテーションでアラインされた
    raw データのシンボル（例: マウス FAM3D の raw データ表記 "Oit1"）と食い違うことが
    ある。Ensembl ID を主キーにしても raw 側の ID が古い/欠落しているケースに備え、
    シンボルの別名（synonym）経由でも解決できるようにするための表。
    1 遺伝子に複数 synonym があれば複数行になる。
    """
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
  <Dataset name="{species}_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="external_synonym"/>
  </Dataset>
</Query>"""
    log.info("BioMart synonyms: %s", species)
    text = _biomart_post(query)
    df = pd.read_csv(io.StringIO(text), sep="\t")
    df.columns = ["ensembl_gene_id", "external_gene_name", "external_synonym"]
    df = df.dropna(subset=["ensembl_gene_id"])
    return df


def fetch_orthologs():
    pairs = [
        (cfg["orthology"]["cat_species"], cfg["orthology"]["mouse_species"], "cat2mouse"),
        (cfg["orthology"]["cat_species"], cfg["orthology"]["human_species"], "cat2human"),
        (cfg["orthology"]["mouse_species"], cfg["orthology"]["human_species"], "mouse2human"),
    ]
    for src, tgt, name in pairs:
        out = REF / f"orthologs_{name}.tsv"
        if out.exists():
            log.info("skip (exists): %s", out.name)
            continue
        for attempt in range(3):
            try:
                df = biomart_orthologs(src, tgt)
                df.to_csv(out, sep="\t", index=False)
                log.info("%s: %d rows, one2one %d", name, len(df),
                         (df.orthology_type == "ortholog_one2one").sum())
                break
            except Exception as exc:  # noqa: BLE001
                log.warning("attempt %d failed: %s", attempt + 1, exc)
                time.sleep(10)
        else:
            log.error("FAILED %s — 手動で BioMart から取得して %s に置いてください", name, out)


def fetch_synonyms():
    """マウス（01_load.py の build_ortholog_map で ID→シンボル解決に使う）と
    ネコ・ヒト（将来同種の問題が出た場合に備え）の synonym 表を取得。"""
    species_map = {
        cfg["orthology"]["mouse_species"]: "mouse",
        cfg["orthology"]["cat_species"]: "cat",
        cfg["orthology"]["human_species"]: "human",
    }
    for species, name in species_map.items():
        out = REF / f"synonyms_{name}.tsv"
        if out.exists():
            log.info("skip (exists): %s", out.name)
            continue
        for attempt in range(3):
            try:
                df = biomart_gene_synonyms(species)
                df.to_csv(out, sep="\t", index=False)
                log.info("%s: %d rows, %d unique genes", name, len(df), df.ensembl_gene_id.nunique())
                break
            except Exception as exc:  # noqa: BLE001
                log.warning("attempt %d failed: %s", attempt + 1, exc)
                time.sleep(10)
        else:
            log.error("FAILED synonyms_%s — 手動で BioMart から取得して %s に置いてください", name, out)


def fetch_entrez_symbols():
    """ヒト Entrez Gene ID -> シンボルの対応表。

    ERCB (GSE104954) は Brainarray ENTREZG カスタム CDF で、
    プローブセット ID が "<EntrezID>_at" 形式（1プローブセット=1遺伝子）。
    シンボルに直すためにこの表が要る。
    """
    out = REF / "entrez2symbol_human.tsv"
    if out.exists():
        log.info("skip (exists): %s", out.name)
        return
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
  <Dataset name="{cfg['orthology']['human_species']}_gene_ensembl" interface="default">
    <Attribute name="entrezgene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="ensembl_gene_id"/>
  </Dataset>
</Query>"""
    log.info("BioMart entrez->symbol (human)")
    for attempt in range(3):
        try:
            text = _biomart_post(query)
            df = pd.read_csv(io.StringIO(text), sep="\t")
            df.columns = ["entrez_id", "symbol", "ensembl_gene_id"]
            df = df.dropna(subset=["entrez_id", "symbol"])
            df["entrez_id"] = df.entrez_id.astype(float).astype(int)
            df = df.drop_duplicates("entrez_id")
            df.to_csv(out, sep="\t", index=False)
            log.info("entrez2symbol: %d 遺伝子", len(df))
            return
        except Exception as exc:  # noqa: BLE001
            log.warning("attempt %d failed: %s", attempt + 1, exc)
            time.sleep(10)
    log.error("FAILED entrez2symbol — 手動で BioMart から取得して %s に置いてください", out)


def fetch_uniprot_secretome(taxon: int, name: str):
    """シグナルペプチド保有 or 細胞外/分泌に注釈された遺伝子を取得。"""
    out = REF / f"secretome_{name}.tsv"
    if out.exists():
        log.info("skip (exists): %s", out.name)
        return
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": f"(taxonomy_id:{taxon}) AND (reviewed:true)",
        "fields": "accession,gene_primary,protein_name,ft_signal,cc_subcellular_location,ft_transmem,length,mass",
        "format": "tsv",
    }
    log.info("UniProt %s", name)
    r = requests.get(url, params=params, timeout=900)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    df.columns = [c.strip() for c in df.columns]
    gene = [c for c in df.columns if "Gene" in c][0]
    sig = [c for c in df.columns if "Signal" in c]
    loc = [c for c in df.columns if "Subcellular" in c]
    tm = [c for c in df.columns if "Transmembrane" in c]
    df["symbol"] = df[gene].astype(str).str.split().str[0].str.upper()
    df["has_signal_peptide"] = df[sig[0]].notna() if sig else False
    locstr = df[loc[0]].fillna("") if loc else ""
    df["is_secreted"] = locstr.str.contains("Secreted", case=False, na=False)
    df["is_extracellular"] = locstr.str.contains("extracellular", case=False, na=False)
    df["n_transmembrane"] = (df[tm[0]].fillna("").str.count("TRANSMEM") if tm else 0)
    keep = ["symbol", "has_signal_peptide", "is_secreted", "is_extracellular",
            "n_transmembrane", "Length", "Mass"]
    keep = [c for c in keep if c in df.columns]
    df[keep].to_csv(out, sep="\t", index=False)
    log.info("%s: %d entries, secreted %d", name, len(df), int(df.is_secreted.sum()))


def fetch_msigdb():
    """gseapy 経由で Hallmark / GO BP / KEGG の GMT をローカルに保存。"""
    try:
        import gseapy as gp
    except ImportError:
        log.error("gseapy が無いため MSigDB は取得しません（pip install gseapy）")
        return
    sets = {
        "hallmark": "MSigDB_Hallmark_2020",
        "gobp": "GO_Biological_Process_2023",
        "kegg": "KEGG_2021_Human",
    }
    for name, lib in sets.items():
        out = REF / f"geneset_{name}.gmt"
        if out.exists():
            log.info("skip (exists): %s", out.name)
            continue
        try:
            d = gp.get_library(name=lib, organism="Human")
            with open(out, "w", encoding="utf-8") as fh:
                for term, genes in d.items():
                    fh.write("\t".join([term, ""] + list(genes)) + "\n")
            log.info("%s: %d terms", name, len(d))
        except Exception as exc:  # noqa: BLE001
            log.warning("%s 取得失敗: %s", name, exc)


if __name__ == "__main__":
    fetch_orthologs()
    fetch_synonyms()
    fetch_entrez_symbols()
    fetch_uniprot_secretome(10090, "mouse")
    fetch_uniprot_secretome(9606, "human")
    fetch_uniprot_secretome(9685, "cat")
    fetch_msigdb()
    log.info("done. data/ref/ を確認してください")
