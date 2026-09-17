# ENVIRONMENT — 環境の再現とデータの再取得

作成: 2026-09-15。`HANDOVER.md` の補足。

---

## 1. 元の環境

| 項目 | 値 |
|---|---|
| OS | macOS 26.6（build 25G72）、Apple Silicon（arm64）、8 コア、メモリ 8 GB |
| リポジトリ | `~/Downloads/ckd_xspecies`（GitHub: masakiwatanabe-labani/ckd-model-distance、branch main） |
| Python | 3.9.6（`/usr/bin/python3` から作った `.venv`） |
| 主要パッケージ | pandas 2.3.3、numpy 2.0.2、scipy 1.13.1、matplotlib 3.9.4、openpyxl 3.1.5。一式は `requirements-lock.txt` |
| BLAST+ | 2.17.0（Homebrew、`/opt/homebrew/bin/blastp`。**`brew` と `/opt/homebrew/bin` は PATH に入っていない**） |
| DIAMOND | 入っていない（`src/31` は blastp で実行した） |
| curl | 8.7.1 |
| git | 2.39.5 |
| シェル | zsh |

---

## 2. Python 環境

```bash
cd ckd_xspecies
python3 -m venv .venv
.venv/bin/pip install -r requirements-lock.txt    # 完全再現（2026-09-15 の pip freeze）
# または
.venv/bin/pip install -r requirements.txt         # 下限だけ指定した版
```

- `requirements-lock.txt` は 2026-09-15 に `pip freeze` で更新した。前の版からの差分は `pypdf==6.16.2` と `typing_extensions==4.16.0` の追加だけ。
- 元の環境の Python は 3.9 で、`urllib3` が LibreSSL について警告を出すが、動作には影響しない。
- **新しい環境で Python 3.10 以降を使う場合**: lock の一部（numpy 2.0.2、scipy 1.13.1 など）は 3.9 向けの版なので、`requirements.txt` から入れて、主要な数値（`results/summary.md`、`projection/manuscript/NUMBERS.md` の `verify_numbers.py`）が再現するか確かめること。

---

## 3. 外部ツール

| ツール | 使う場所 | 元の環境 | インストール例 |
|---|---|---|---|
| BLAST+ | `src/31_urine_reannotation.py` | 2.17.0 | macOS: `brew install blast`。Linux: `conda install -c bioconda blast` か NCBI の配布バイナリ（`https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/`） |
| DIAMOND（任意） | `src/31` を DIAMOND に置き換える場合 | なし | `conda install -c bioconda diamond`、または GitHub（bbuchfink/diamond）の Linux 用バイナリ。**blastp の結果と閾値を揃えて比べてから切り替えること** |
| ProteoWizard msconvert | 段階2 | なし | Windows 版、または Docker イメージ（`PHASE2_PLAN.md` 2節） |
| FragPipe（MSFragger、Philosopher、IonQuant） | 段階2 | なし | GitHub（Nesvilab/FragPipe）。Java が必要。MSFragger は学術用ライセンスへの同意が必要 |
| MaxQuant 2.6.4.0 | 段階2（論文の再現） | なし | maxquant.org から登録してダウンロード。.NET が必要 |
| Skyline | 段階2（標的抽出） | なし | Windows |
| NCBI datasets CLI | 段階2（RefSeq の取得） | なし | `conda install -c conda-forge ncbi-datasets-cli` |
| bcftools ≥ 1.21、samtools、SpliceAI | コットンラット（サーバー側） | このマシンには無い | — |

`src/31_urine_reannotation.py` は、環境変数で BLAST の場所とスレッド数を変えられる。

```bash
BLAST_BIN=/usr/bin BLAST_THREADS=16 .venv/bin/python src/31_urine_reannotation.py
```

---

## 4. データの再取得

**完全な一覧とチェックサム**は `DATA_CHECKSUMS.md5`（`data/raw`、`data/external`、`data/ref` の57ファイル）。確認するときは次のコマンドを使う。

```bash
# macOS
md5 -r $(awk '{print $2}' DATA_CHECKSUMS.md5) | diff - DATA_CHECKSUMS.md5
# Linux
md5sum -c <(awk '{print $1"  "$2}' DATA_CHECKSUMS.md5)
```

`data/interim/` は再生成できるので対象外（`src/01_load.py` と `src/31` の BLAST）。

### 4.1 自動で取れるもの

```bash
.venv/bin/python src/00_fetch_refs.py        # data/ref/ の orthologs_*, secretome_*, synonyms_*, geneset_*.gmt, entrez2symbol_human.tsv
.venv/bin/python src/00b_fetch_external.py   # data/external/GSE98622, GSE79443, GSE104954
```

- Ensembl BioMart、MSigDB、GEO のバージョンが変わると、チェックサムは一致しなくなる。
- 旧稿の数値は `data/README.md` に書かれた版で計算されている。

### 4.2 手で取るもの

| ファイル | 取得方法 | 旧稿で使うか |
|---|---|---|
| `data/raw/42003_2025_9164_MOESM3_ESM.xlsx` | Li et al. Commun Biol 2025 の Supplementary Data 3 | 使う |
| `data/raw/all_fpkm_TRECK.xlsx` | GSE299326 | 使う |
| `data/raw/PK25055_解析結果.xlsx`、`PK25055-1〜3_全タンパク質.xlsx` | 研究室内（非公開） | 使う |
| `data/external/KPMP/DataLake_DEPs.txt` | KPMP Atlas Explorer（doi:10.48698/mg7h-bc51）、利用規約に同意 | 使う |
| `data/external/plasma_RamirezMedina2023/12014_2023_9405_MOESM2_ESM.csv` | Ramírez Medina et al. Clin Proteomics 2023 の補足 MOESM2 | 30 |
| `data/external/HPA/hpa_secretome_blood.tsv` | 下のコマンド | 30 |

### 4.3 尿プロテオームと参照配列（コマンド）

```bash
# PXD076696 の処理済みファイルとメタデータ（RAW は PHASE2_PLAN.md）
B=https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/08/PXD076696
mkdir -p data/external/urine_PXD076696/pride && cd data/external/urine_PXD076696/pride
for f in MSvsStudyGroup_IDmask_20260828.xlsx MSvsStudyGroup_IDmask.xlsx 20260829_Notes.txt catUniprot_withADH1.fasta \
         proteinGroups_20251010_IDmask.xlsx proteinGroups_20251016_IDmask.xlsx proteinGroups_20251017_IDmask.xlsx \
         proteinGroups_20251024_IDmask.xlsx proteinGroups_20251103_IDmask.xlsx proteinGroups_20251107_IDmask.xlsx \
         proteinGroups_20251117_IDmask.xlsx proteinGroups_20260301_IDmask.xlsx; do curl -sS -O "$B/$f"; done
cd -

# 論文本文と補足（Europe PMC）
mkdir -p data/external/urine_PXD076696/paper
curl -sS "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13498864/fullTextXML" -o data/external/urine_PXD076696/paper/PMC13498864.xml
curl -sS -L "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13498864/supplementaryFiles" -o /tmp/suppl.zip   # docx を取り出す
# appendixS1_table_S*.csv と fulltext.txt は、2026-09-14 のセッションで docx / XML から抽出して作った（抽出コードは保存していない）

# UniProt reviewed（HTTP/2 だと途中で切れることがあるので --http1.1）
for sp in human:9606 mouse:10090; do n=${sp%%:*}; t=${sp##*:}
  curl --http1.1 -sS -D data/ref/uniprot_reviewed_${n}.headers \
    "https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A${t}%29&format=fasta&compressed=true" \
    -o data/ref/uniprot_reviewed_${n}.fasta.gz && gunzip -f data/ref/uniprot_reviewed_${n}.fasta.gz
done
grep -i x-uniprot-release data/ref/uniprot_reviewed_*.headers    # 元の環境は 2026_03

# Human Protein Atlas（2026-09-15 に取得した版とは中身が変わりうる）
mkdir -p data/external/HPA
curl -sS "https://www.proteinatlas.org/api/search_download.php?search=&format=tsv&compress=no&columns=g,eg,up,pc,scl,sec,secl,sl,sa,sa_location,secretome,secloc,blconcms,blconcib" \
  -o data/external/HPA/hpa_secretome_blood.tsv
```

**自動化されていないもの**（再取得するならスクリプト化を勧める）

| ファイル | 元の作り方 |
|---|---|
| `data/external/plasma_RamirezMedina2023/panel_uniprot_gene.tsv` | 血漿パネル623 accession を UniProt REST（`rest.uniprot.org/uniprotkb/search`、`fields=accession,gene_primary,reviewed,organism_id`）に90件ずつ問い合わせた |
| `data/external/urine_PXD076696/paper/appendixS1_table_S*.csv` | 補足 docx の表をそのまま CSV にした |
| `data/external/urine_PXD076696/pride/pride_filelist_20260914.json`、`raw_d_manifest.tsv` | PRIDE API のファイル一覧（2ページ分）を保存し、`.d.tar.gz` だけを抜き出した |

---

## 5. Claude Code で作業するときの注意（元の環境で踏んだもの）

`HANDOVER.md` 4.4節と同じ内容。

- zsh では `grep --include='*.csv'` のようにグロブを引用する。
- 前景の `sleep` は使えない。長い処理は `run_in_background` で走らせる。
- PMC 本体、Wiley、figshare はボット対策で取れない。回避を試みない。
- 結果ファイルを上書きしない（ユーザーの基準）。訂正は追記か別ファイルで行う。
