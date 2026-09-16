"""Step 0: Δ行列（遺伝子 × 状態）を組み立てる。

出力:
    delta_matrix.tsv            行=ヒトシンボル（大文字）, 列=状態, 値=Δlog2FC
    time_map.tsv                state, hours（ネコ2状態は NA）
    groupA_intersection.txt     主解析の遺伝子集合（cortex ∩ medulla ∩ マウス蛋白）
    groupA_union.txt            感度分析の遺伝子集合（(cortex ∪ medulla) ∩ マウス蛋白）
    ensembl_ortholog_genes.txt  CHECK 1 の core 集合（全データセットで ortholog 由来）
    control_log2cpm.tsv         各データセットの対照群の平均発現（発現量マッチ用）
    gene_provenance.tsv         遺伝子ごとの写像由来（ortholog / fallback）
    results/build_report.md     実際に得られた n と、想定値とのズレ

設計上の判断（結果を左右するので明示する）:

  1. 層は RNA に統一する。IRI と Pod-TRECK Day5 には対応するプロテオームが無く、
     さらに CLAUDE.md の落とし穴2のとおりネコのプロテオーム値は蛋白ごとに
     行中心化されているため |v| が定義できない。alpha / nrm は絶対値を使う
     指標なので、順位でしか読めない層は使えない。
     Group A は「プロテオームで検出されたか」という遺伝子集合の定義にのみ用い、
     Δ の値は RNA から取る。

  2. 共通空間はヒトシンボル。写像は Ensembl one2one オルソログ表を第一とし、
     引けない場合のみ大文字化にフォールバックする（src/09_dataset_matrix.py の
     to_human_space と完全に同じ挙動）。フォールバック分は
     gene_provenance.tsv で "fallback" と記録され、CHECK 1 の added 集合になる。
     主解析（Step 2）は Group A に絞られ、Group A は蛋白検出とオルソログの
     両方を通った遺伝子なので、実質 core のみで走る。

  3. IRI の対照は時点ごとに変える。src/14_time_axis_full.py と同じ:
     2h〜28d は若齢 sham（SHAM4h+SHAM24h, n=6）、12mo は同週齢 SHAM12m（n=3）。
     全 sham をプールすると IRI12m だけが「加齢+傷害 vs 若齢」になり、
     加齢シグナルが時間差として混入して時間効果を過大評価する。
     GPL19057 の IRI6mN は別プラットフォームなので除外し、9時点すべてを
     GPL13112 内に揃える。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))
from lib_stats import load_config, get_logger, moderated_ttest  # noqa: E402

log = get_logger("projection_build")
cfg = load_config(ROOT / "config" / "config.yaml")
INT = ROOT / cfg["paths"]["interim"]
REF = ROOT / cfg["paths"]["ref"]
TH = cfg["thresholds"]

# ---------------------------------------------------------------- 状態の定義
# IRI: GPL13112 の9時点のみ。(ラベル, 経過時間[hours], 対照ラベル群)
YOUNG_SHAM = ["SHAM4h", "SHAM24h"]
IRI_STATES = [
    ("IRI_2h",   "IRI2h",  2.0,    YOUNG_SHAM),
    ("IRI_4h",   "IRI4h",  4.0,    YOUNG_SHAM),
    ("IRI_24h",  "IRI24h", 24.0,   YOUNG_SHAM),
    ("IRI_48h",  "IRI48h", 48.0,   YOUNG_SHAM),
    ("IRI_72h",  "IRI72h", 72.0,   YOUNG_SHAM),
    ("IRI_7d",   "IRI7d",  168.0,  YOUNG_SHAM),
    ("IRI_14d",  "IRI14d", 336.0,  YOUNG_SHAM),
    ("IRI_28d",  "IRI28d", 672.0,  YOUNG_SHAM),
    ("IRI_12mo", "IRI12m", 8760.0, ["SHAM12m"]),
]
# Pod-TRECK: (状態名, 群ラベル, 経過時間[hours])
PODTRECK_STATES = [("mouse_5D", "DT100ng-5D", 120.0),
                   ("mouse_2W", "DT100ng-2W", 336.0),
                   ("mouse_3W", "DT100ng-3W", 504.0)]
# ネコ: 自然発症のため経過時間は定義できない -> hours = NA
# (状態名, 群ラベル, 区画)
CAT_STATES = [("cat_CKD12", "CKD1/2", "ctx"), ("cat_CKD34", "CKD3/4", "ctx"),
              ("cat_med_CKD12", "CKD1/2", "med"), ("cat_med_CKD34", "CKD3/4", "med")]
# 補助ベクトル: 進行軸（CKD3/4 vs CKD1/2）。対照群比ではないので Δ行列には入れず
# delta_matrix_aux.tsv に分ける。旧稿の「皮質と髄質の進行軸が直交」(rho 0.204) は
# この軸の話であり、CKD3/4 vs Control 同士の相関 (rho 0.763) とは別物。
CAT_AUX = [("cat_ctx_prog", "CKD3/4", "CKD1/2", "ctx"),
           ("cat_med_prog", "CKD3/4", "CKD1/2", "med")]


# ---------------------------------------------------------------- オルソログ
def one2one(fname: str) -> dict:
    """BioMart オルソログ表から src_symbol -> tgt_symbol の一意な辞書。

    src/09_dataset_matrix.py の _one2one_map と同一のタイブレーク
    （confidence 降順 -> perc_id 降順）を使う。ここを変えると旧稿の
    rho 0.490 が再現しなくなり、CHECK 1 の比較対象が消える。
    """
    o = pd.read_csv(REF / fname, sep="\t")
    o = o[o.orthology_type.isin(cfg["orthology"]["keep_types"])]
    o = o.dropna(subset=["src_symbol", "tgt_symbol"]).copy()
    o["rank"] = (-o.confidence.fillna(0)) * 1000 - o.perc_id.fillna(0)
    o = o.sort_values("rank").drop_duplicates("src_symbol")
    return dict(zip(o.src_symbol, o.tgt_symbol))


CAT2HUMAN = one2one("orthologs_cat2human.tsv")
MOUSE2HUMAN = one2one("orthologs_mouse2human.tsv")
MAPS = {"cat": CAT2HUMAN, "mouse": MOUSE2HUMAN}


def to_human(s: pd.Series, species: str) -> tuple[pd.Series, pd.Series]:
    """系列をヒトシンボル空間へ。値と、遺伝子ごとの写像由来を返す。

    src/09_dataset_matrix.py の to_human_space と同じ順序で重複を落とす
    （元の index 順で最初に現れたものを残す）。
    """
    mapping = MAPS[species]
    src = [str(g) for g in s.index]
    human = [str(mapping[g]).upper() if g in mapping else str(g).upper() for g in src]
    prov = ["ortholog" if g in mapping else "fallback" for g in src]
    out = pd.DataFrame({"value": s.to_numpy(), "prov": prov, "src": src}, index=human)
    out = out[~out.index.duplicated(keep="first")]
    return out["value"], out["prov"]


# ---------------------------------------------------------------- Δ の計算
def read_int(name: str):
    df = pd.read_parquet(INT / f"{name}.parquet")
    return df.iloc[:, 0] if name.endswith("grp") else df


def _cat_mat(tissue: str):
    mat = read_int(f"cat_rna_{tissue}")
    grp = read_int(f"cat_rna_{tissue}_grp").squeeze()
    return mat[mat.std(axis=1) > 0], grp     # 02_de.py と同じフィルタ


def cat_deltas() -> dict[str, pd.Series]:
    """ネコ皮質・髄質 RNA。シート値はすでに log2 スケール（範囲 -5〜16.6）。"""
    out = {}
    for state, label, tissue in CAT_STATES:
        mat, grp = _cat_mat(tissue)
        case = mat[grp[grp == label].index]
        ctrl = mat[grp[grp == "Control"].index]
        res = moderated_ttest(case, ctrl)
        out[state] = res["lfc"]
        log.info("%s: case n=%d vs Control n=%d, 遺伝子 %d",
                 state, case.shape[1], ctrl.shape[1], len(res))
    return out


def cat_aux_deltas() -> dict[str, pd.Series]:
    """進行軸（CKD3/4 vs CKD1/2）。対照解析の下振れ対照に使う。"""
    out = {}
    for state, a, b, tissue in CAT_AUX:
        mat, grp = _cat_mat(tissue)
        res = moderated_ttest(mat[grp[grp == a].index], mat[grp[grp == b].index])
        out[state] = res["lfc"]
        log.info("%s: %s n=%d vs %s n=%d, 遺伝子 %d", state, a,
                 int((grp == a).sum()), b, int((grp == b).sum()), len(res))
    return out


def podtreck_deltas() -> dict[str, pd.Series]:
    """Pod-TRECK RNA。FPKM -> 群平均>=1 のフィルタ後に log2(FPKM+1)（02_de.py と同じ）。"""
    mat = read_int("mouse_rna")
    grp = read_int("mouse_rna_grp").squeeze()
    gm = mat.T.groupby(grp).mean().T
    keep = gm.max(axis=1) >= TH["fpkm_min"]
    log.info("Pod-TRECK 発現フィルタ: %d / %d 遺伝子を保持", int(keep.sum()), len(mat))
    m = np.log2(mat[keep] + 1)
    ctrl = m[grp[grp == "Ctrl"].index]
    out = {}
    for state, label, _h in PODTRECK_STATES:
        case = m[grp[grp == label].index]
        res = moderated_ttest(case, ctrl)
        out[state] = res["lfc"]
        log.info("%s: case n=%d vs Ctrl n=%d, 遺伝子 %d",
                 state, case.shape[1], ctrl.shape[1], len(res))
    return out


def iri_deltas() -> dict[str, pd.Series]:
    """GSE98622 の9時点。時点ごとに対照を変える（14_time_axis_full.py と同じ）。"""
    mat = read_int("mouse_iri_matrix")
    grp = read_int("mouse_iri_grp").squeeze()
    out = {}
    for state, label, _h, ctrl_labels in IRI_STATES:
        A = mat[grp[grp == label].index]
        B = mat[grp[grp.isin(ctrl_labels)].index]
        A, B = np.log2(A + 1), np.log2(B + 1)
        keep = pd.concat([A, B], axis=1).max(axis=1) > 1
        res = moderated_ttest(A[keep], B[keep])
        out[state] = res["lfc"]
        log.info("%s: case n=%d vs %s n=%d, 遺伝子 %d",
                 state, A.shape[1], ctrl_labels, B.shape[1], len(res))
    return out


# ---------------------------------------------------------------- 対照発現量
def control_expression() -> pd.DataFrame:
    """各データセットの対照群平均発現（ヒト空間）。CHECK 1 の発現量マッチに使う。

    単位はデータセットごとに異なる（ネコ=シートの log2 正規化値、
    マウス=log2(FPKM+1)、IRI=log2(正規化カウント+1)）。CHECK 1 では
    データセット内のパーセンタイル順位に変換してから使うので、単位の
    違いは問題にならない。列間で絶対値を比較してはいけない。
    """
    cols = {}
    mat = read_int("cat_rna_ctx"); grp = read_int("cat_rna_ctx_grp").squeeze()
    v, _ = to_human(mat[grp[grp == "Control"].index].mean(axis=1), "cat")
    cols["cat_ctx_control"] = v

    mat = read_int("mouse_rna"); grp = read_int("mouse_rna_grp").squeeze()
    gm = mat.T.groupby(grp).mean().T
    m = np.log2(mat[gm.max(axis=1) >= TH["fpkm_min"]] + 1)
    v, _ = to_human(m[grp[grp == "Ctrl"].index].mean(axis=1), "mouse")
    cols["podtreck_control"] = v

    mat = read_int("mouse_iri_matrix"); grp = read_int("mouse_iri_grp").squeeze()
    B = np.log2(mat[grp[grp.isin(YOUNG_SHAM)].index] + 1)
    v, _ = to_human(B[B.max(axis=1) > 1].mean(axis=1), "mouse")
    cols["iri_young_sham"] = v
    return pd.DataFrame(cols)


# ---------------------------------------------------------------- Group A
def group_a_sets() -> tuple[list[str], list[str], dict]:
    """プロテオームで検出された遺伝子集合をヒト空間で作る。

    検出基準: そのデータセットのサンプルの 50% 以上で非欠測。
    """
    frac = 0.5

    def detected(name: str, species: str) -> set:
        d = read_int(name)
        s = pd.Series(1.0, index=d.index[d.notna().mean(axis=1) >= frac])
        h, _ = to_human(s, species)
        return set(h.index)

    ctx = detected("cat_prot_ctx", "cat")
    med = detected("cat_prot_med", "cat")
    mo = detected("mouse_prot", "mouse")
    inter = sorted((ctx & med) & mo)
    union = sorted((ctx | med) & mo)
    stats = {"cat_cortex_detected": len(ctx), "cat_medulla_detected": len(med),
             "podtreck_prot_detected": len(mo),
             "cat_ctx_and_med": len(ctx & med), "cat_ctx_or_med": len(ctx | med),
             "groupA_intersection": len(inter), "groupA_union": len(union)}
    log.info("Group A: %s", stats)
    return inter, union, stats


# ---------------------------------------------------------------- 組み立て
def main():
    series, prov = {}, {}
    for state, s in cat_deltas().items():
        series[state], prov[state] = to_human(s[np.isfinite(s)], "cat")
    for state, s in podtreck_deltas().items():
        series[state], prov[state] = to_human(s[np.isfinite(s)], "mouse")
    for state, s in iri_deltas().items():
        series[state], prov[state] = to_human(s[np.isfinite(s)], "mouse")

    delta = pd.DataFrame(series)
    delta.index.name = "gene"
    delta = delta.sort_index()

    # 写像由来: あるデータセットに現れる遺伝子がすべて ortholog 由来なら core。
    P = pd.DataFrame(prov).reindex(delta.index)
    is_ortho = P.eq("ortholog") | P.isna()
    core = sorted(delta.index[is_ortho.all(axis=1) & P.notna().any(axis=1)])
    P.index.name = "gene"

    inter, union, ga_stats = group_a_sets()
    inter_in = [g for g in inter if g in delta.index]
    union_in = [g for g in union if g in delta.index]

    # 全状態で有限な遺伝子（Step 2 はここに絞られる。|v| を状態間で比べるため）
    complete = delta.index[delta.notna().all(axis=1)]

    aux_series = {}
    for state, sv in cat_aux_deltas().items():
        aux_series[state], _ = to_human(sv[np.isfinite(sv)], "cat")
    pd.DataFrame(aux_series).rename_axis("gene").sort_index().to_csv(
        HERE / "delta_matrix_aux.tsv", sep="\t", na_rep="NA")

    tm = pd.DataFrame(
        [{"state": s, "hours": np.nan} for s, _l, _t in CAT_STATES]
        + [{"state": s, "hours": h} for s, _l, h in PODTRECK_STATES]
        + [{"state": s, "hours": h} for s, _l, h, _c in IRI_STATES])

    delta.to_csv(HERE / "delta_matrix.tsv", sep="\t", na_rep="NA")
    tm.to_csv(HERE / "time_map.tsv", sep="\t", index=False, na_rep="NA")
    P.to_csv(HERE / "gene_provenance.tsv", sep="\t", na_rep="NA")
    control_expression().to_csv(HERE / "control_log2cpm.tsv", sep="\t", na_rep="NA")
    for fn, genes in [("groupA_intersection.txt", inter_in),
                      ("groupA_union.txt", union_in),
                      ("ensembl_ortholog_genes.txt", core)]:
        (HERE / fn).write_text("\n".join(genes) + "\n")

    n_fb = int((~is_ortho).any(axis=1).sum())
    rep = [
        "# Step 0 build report",
        "",
        f"- Δ行列: {delta.shape[0]} 遺伝子 x {delta.shape[1]} 状態",
        f"- 全状態で有限な遺伝子（complete case）: {len(complete)}",
        f"- core（全データセットで Ensembl ortholog 由来）: {len(core)}",
        f"- added（どこかで大文字フォールバック）: {n_fb}",
        "",
        "## Group A",
        "",
        "| 集合 | n（Δ行列内） | n（写像後の素の集合） |",
        "|---|---|---|",
        f"| intersection: (cortex ∩ medulla) ∩ mouse prot | {len(inter_in)} | {ga_stats['groupA_intersection']} |",
        f"| union: (cortex ∪ medulla) ∩ mouse prot | {len(union_in)} | {ga_stats['groupA_union']} |",
        "",
        "中間量: "
        + ", ".join(f"{k}={v}" for k, v in ga_stats.items()
                    if k not in ("groupA_intersection", "groupA_union")),
        "",
        "## 状態あたりの有限値数",
        "",
        "| 状態 | 有限値数 |",
        "|---|---|",
        *[f"| {k} | {v} |" for k, v in delta.notna().sum().items()],
    ]
    (HERE / "results" / "build_report.md").write_text("\n".join(rep) + "\n")
    log.info("書き出し完了: %s", HERE)
    print("\n".join(rep))


if __name__ == "__main__":
    main()
