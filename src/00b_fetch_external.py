"""外部データセットの取得：マウスIRI（GSE98622）と KPMP（ヒト）。

設計上の役割:
  - GSE98622 = 尿細管入口のマウス。Pod-TRECK（糸球体入口）と対にすることで
    「同一種内で入口を変えたときの距離」が測れる。これが無いと、ネコ×マウスの
    種間距離に病因差が交絡したままになる。
  - KPMP = ヒト。区画分離（糸球体 G / 尿細管間質 TI）されているのが利点。

KPMP の注意:
  regional transcriptomics で公開されているのは「区画特異性スコア」であって
  CKD vs 対照の直接比較ではない。CKD と HRT のスコア差を取る必要があり、
  通常の fold change とは別物。トランスクリプトームで正面から比較するなら
  snRNA-seq の疑似バルクを使うほうが素直。proteomics は pre-computed の
  CKD vs HRT コントラストがそのまま使える。
"""
from __future__ import annotations
import gzip
import io
import shutil
import sys
import tarfile
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import load_config, get_logger  # noqa: E402

log = get_logger("00b_fetch_external")
cfg = load_config()
ROOT = Path(cfg["_root"])
EXT = ROOT / "data" / "external"
EXT.mkdir(parents=True, exist_ok=True)

GEO_SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/suppl/"


def fetch_geo_supplementary(acc: str, keep: list[str] | None = None) -> Path:
    """GEO series の supplementary ファイル一覧を取得してダウンロードする。

    keep を渡すと、その文字列を名前に含むファイルだけを取得する。
    CEL の塊である *_RAW.tar のような巨大ファイルを避けるために使う。
    """
    out = EXT / acc
    out.mkdir(exist_ok=True)
    stub = acc[:-3]
    base = GEO_SUPP.format(stub=stub, acc=acc)
    log.info("GEO supplementary: %s", base)
    try:
        idx = requests.get(base, timeout=120).text
    except Exception as exc:  # noqa: BLE001
        log.error("一覧取得に失敗: %s。ブラウザで %s を開いて手動DLし %s に置いてください",
                  exc, base, out)
        return out
    names = sorted({n.split('"')[0] for n in idx.split('href="')[1:] if not n.startswith("?")})
    for n in names:
        if n in ("/", "..") or n.endswith("/") or n.startswith(("http://", "https://")):
            continue
        if keep and not any(k in n for k in keep):
            log.info("skip (filtered): %s", n)
            continue
        dest = out / n
        if dest.exists():
            log.info("skip (exists): %s", n)
            continue
        log.info("download %s", n)
        try:
            r = requests.get(base + n, timeout=1800, stream=True)
            r.raise_for_status()
            with open(dest, "wb") as fh:
                shutil.copyfileobj(r.raw, fh)
        except Exception as exc:  # noqa: BLE001
            log.warning("%s の取得に失敗: %s", n, exc)
    return out


def load_gse98622(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """カウント/FPKM 行列とサンプル群ラベルを返す。

    ファイル構成は GEO 側の都合で変わるため、見つかった表形式ファイルを
    総当たりで読み、遺伝子×サンプルらしい行列を採用する。
    サンプル名から時点（2h..12mo）と処置（IRI/sham/normal）を推定する。
    """
    text_cands = [p for p in path.iterdir()
                  if p.suffix in (".txt", ".csv", ".tsv") or p.name.endswith((".txt.gz", ".csv.gz", ".tsv.gz"))]
    xlsx_cands = [p for p in path.iterdir() if p.suffix == ".xlsx"]
    best = None  # (path, df, symbol_col_name_or_None)
    for p in text_cands:
        try:
            opener = gzip.open if p.name.endswith(".gz") else open
            with opener(p, "rt", encoding="utf-8", errors="ignore") as fh:
                head = fh.readline()
            sep = "\t" if head.count("\t") >= head.count(",") else ","
            df = pd.read_csv(p, sep=sep, index_col=0, low_memory=False)
            df_num = df.select_dtypes("number")
            if df_num.shape[0] > 5000 and df_num.shape[1] >= 20:
                if best is None or df_num.shape[1] > best[1].shape[1]:
                    best = (p, df_num, None)
        except Exception:  # noqa: BLE001
            continue
    for p in xlsx_cands:
        try:
            xl = pd.ExcelFile(p)
            for sheet in xl.sheet_names:
                df = xl.parse(sheet)
                sym_col = next((c for c in df.columns
                                if str(c).strip().lower() in ("symbol", "gene", "gene_symbol", "genesymbol")), None)
                df_num = df.select_dtypes("number")
                if df_num.shape[0] > 5000 and df_num.shape[1] >= 20:
                    if best is None or df_num.shape[1] > best[1].shape[1]:
                        best = (p, df_num, df[sym_col] if sym_col else None)
        except Exception:  # noqa: BLE001
            continue
    if best is None:
        raise FileNotFoundError(
            f"{path} に遺伝子×サンプル行列が見つかりません。"
            "GEO の supplementary を手動で確認してください。")
    src_path, mat, sym_series = best
    if sym_series is not None:
        # 遺伝子シンボル列をインデックスに採用（Ensembl ID より下流の解析で扱いやすい）
        mat = mat.set_axis(sym_series.values, axis=0)
    log.info("GSE98622 行列: %s %s", src_path.name, mat.shape)

    def label(col: str) -> str:
        """列名の replicate 接尾辞 (-R\\d+) を除いた部分をそのまま群名にする。

        以前はタイムポイントのタグ文字列一致（'2h' 等）を使っていたが、
        'IRI72h' が '2h' に、'IRI28d' が unassigned になるなど誤判定を生んでいた
        （GSE98622 実データの列名: SHAM4h/24h/12m, NORM3m/9m/15m,
        IRI2h/4h/24h/48h/72h/7d/14d/28d/6mN/12m）。
        列名そのものを群キーとして使う方が確実。
        """
        import re
        c = re.sub(r"-R\d+$", "", str(col))
        return c

    groups = pd.Series({c: label(c) for c in mat.columns}, name="group")
    unassigned = int((groups == "unassigned").sum())
    if unassigned:
        log.warning("群ラベルを推定できない列が %d 本あります。"
                    "data/external/GSE98622/ の series matrix を見て "
                    "config で明示指定してください", unassigned)
    return mat, groups


KPMP_URLS = {
    # 論文で使用したもの。URL が変わることがあるので必ずログを確認すること。
    "regional_proteomics": "https://doi.org/10.48698/mg7h-bc51",
    "regional_transcriptomics": "https://doi.org/10.48698/t9fh-qn48",
}


def kpmp_instructions():
    log.info(
        "KPMP は Atlas Explorer から手動取得してください（クリックスルー同意が必要）。\n"
        "  regional proteomics : %s  → data/external/KPMP/regional_proteomics.csv\n"
        "  regional transcriptomics : %s → data/external/KPMP/regional_transcriptomics.csv\n"
        "  proteomics は CKD.vs.HRT.in.TI / CKD.vs.HRT.in.G の pre-computed コントラストを"
        "そのまま使う。transcriptomics は区画特異性スコアなので CKD-HRT の差を取る。",
        KPMP_URLS["regional_proteomics"], KPMP_URLS["regional_transcriptomics"])
    (EXT / "KPMP").mkdir(exist_ok=True)


GEO_MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/matrix/"


def fetch_geo_matrix(acc: str, keep: list[str] | None = None) -> Path:
    """GEO series matrix ファイルを取得して gunzip する。

    keep を渡すとその名前を含むファイルだけを落とす（巨大ファイル回避）。
    """
    import gzip as _gz
    out = EXT / acc
    out.mkdir(parents=True, exist_ok=True)
    base = GEO_MATRIX.format(stub=acc[:-3], acc=acc)
    log.info("GEO series matrix: %s", base)
    try:
        idx = requests.get(base, timeout=120).text
    except Exception as exc:  # noqa: BLE001
        log.error("一覧取得に失敗: %s", exc)
        return out
    names = sorted({n.split('"')[0] for n in idx.split('href="')[1:]})
    for n in names:
        if not n.endswith(".txt.gz") or n.startswith(("http://", "https://")):
            continue
        if keep and not any(k in n for k in keep):
            continue
        dest = out / n[:-3]
        if dest.exists():
            log.info("skip (exists): %s", dest.name)
            continue
        log.info("download %s", n)
        try:
            r = requests.get(base + n, timeout=1800)
            r.raise_for_status()
            dest.write_bytes(_gz.decompress(r.content))
        except Exception as exc:  # noqa: BLE001
            log.warning("%s の取得に失敗: %s", n, exc)
    return out


def fetch_gse79443():
    """UUO（GSE79443）。正規化カウント行列と series matrix を取得する。"""
    # *_RAW.tar（CEL の塊）は本解析で使わないので取得しない
    p = fetch_geo_supplementary("GSE79443", keep=["norm_counts"])
    for f in list(p.glob("*.gz")):
        import gzip as _gz
        dest = p / f.name[:-3]
        if not dest.exists():
            dest.write_bytes(_gz.decompress(f.read_bytes()))
            log.info("gunzip %s", dest.name)
        f.unlink()
    fetch_geo_matrix("GSE79443")
    return p


def fetch_gse104954():
    """ERCB 尿細管間質（GSE104954）。series matrix のみでよい。

    supplementary の RAW.tar は CEL ファイルの塊で数百MBあり、本解析では
    使わないので取得しない。
    """
    return fetch_geo_matrix("GSE104954")


if __name__ == "__main__":
    p = fetch_geo_supplementary("GSE98622")
    try:
        mat, groups = load_gse98622(p)
        mat.to_parquet(ROOT / cfg["paths"]["interim"] / "mouse_iri_matrix.parquet")
        groups.to_frame().to_parquet(ROOT / cfg["paths"]["interim"] / "mouse_iri_grp.parquet")
        log.info("群構成: %s", groups.value_counts().to_dict())
    except Exception as exc:  # noqa: BLE001
        log.error("GSE98622 の読み込みに失敗: %s", exc)
    fetch_gse79443()
    fetch_gse104954()
    kpmp_instructions()
