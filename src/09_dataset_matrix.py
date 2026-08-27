"""種 × 発症区画（入口）の 2×2 設計による分散分解。

これが本パイプラインの主結果になる想定。

設計:
    |        | 糸球体入口        | 尿細管入口          |
    | ヒト   | KPMP (DKD/HTN主体) | —                  |
    | ネコ   | —                 | 自然発症CKD         |
    | マウス | Pod-TRECK         | IRI (GSE98622 後期) |

マウスだけが両方の入口を持つ。ここから

    d_entry   = マウス内で Pod-TRECK vs IRI の距離        （病因の効果）
    d_species = 入口を揃えた ネコ自然発症 vs マウスIRI     （種の効果）

を別々に測り、大小を比較する。従来の「ネコ(尿細管) × マウス(糸球体)」だけの
比較は両者が交絡しており、この分解なしには解釈できない。

必ず先に走らせるチェック（優先度最高）:
    GSE98622 の後期時点で適応免疫マーカー（Cd3e, Cd8a, Ms4a1, Cd79a, Mzb1）が
    動いているか。動いていれば「マウスに適応免疫軸が欠落」という主張は
    「Pod-TRECK 3週が急性すぎるだけ」に修正される。→ check_adaptive_immunity()
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import (load_config, get_logger, moderated_ttest,  # noqa: E402
                       append_summary)

log = get_logger("09_dataset_matrix")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
EXT = ROOT / "data" / "external"
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]

DE = pd.read_parquet(INT / "de_all.parquet")
OMAP = pd.read_csv(INT / "ortholog_map.tsv", sep="\t")
CAT2MOUSE = dict(zip(OMAP.cat_symbol, OMAP.mouse_symbol))


def _one2one_map(fname: str) -> dict:
    """BioMart オルソログ表から src_symbol -> tgt_symbol の一意な辞書を作る。"""
    p = Path(cfg["_root"]) / cfg["paths"]["ref"] / fname
    if not p.exists():
        log.warning("%s がありません。大文字化のみで写像します", fname)
        return {}
    o = pd.read_csv(p, sep="\t")
    o = o[o.orthology_type.isin(cfg["orthology"]["keep_types"])]
    o = o.dropna(subset=["src_symbol", "tgt_symbol"]).copy()
    o["rank"] = (-o.confidence.fillna(0)) * 1000 - o.perc_id.fillna(0)
    o = o.sort_values("rank").drop_duplicates("src_symbol")
    return dict(zip(o.src_symbol, o.tgt_symbol))


CAT2HUMAN = _one2one_map("orthologs_cat2human.tsv")
MOUSE2HUMAN = _one2one_map("orthologs_mouse2human.tsv")

# 種 × 入口のメタ情報。共通のヒトシンボル空間に写像して比較する。
# 層をまたぐ比較（RNA×蛋白）は距離行列に入れない。測定量が違うものを
# 同じ距離尺度に混ぜると、種・入口の効果と測定層の効果が交絡するため。
# 層ごとに独立した行列を作り、それぞれ答えられる問いだけを答える。
DESIGN_RNA = [
    # key,                    species, entry,        source
    ("cat_natural_ctx_rna",   "cat",   "tubular",    "internal:cat_rna_ctx_late"),
    ("cat_natural_med_rna",   "cat",   "tubular",    "internal:cat_rna_med_late"),
    ("mouse_podtreck_rna",    "mouse", "glomerular", "internal:m_rna_2w"),
    ("mouse_iri_late_rna",    "mouse", "tubular",    "external:GSE98622"),
]
# 蛋白層にはマウスIRIの対応データが存在しないため 2x2 は完成しない。
# したがって種×入口の分解は RNA 層でのみ行い、蛋白層はヒトへの距離に使う。
DESIGN_PROT = [
    ("cat_natural_ctx_prot",  "cat",   "tubular",    "internal:cat_prot_ctx_late"),
    ("cat_natural_med_prot",  "cat",   "tubular",    "internal:cat_prot_med_late"),
    ("mouse_podtreck_prot",   "mouse", "glomerular", "internal:m_prot_d21"),
    ("human_kpmp_ti",         "human", "glomerular", "external:KPMP_TI"),
    ("human_kpmp_g",          "human", "glomerular", "external:KPMP_G"),
]
DESIGN = DESIGN_RNA + DESIGN_PROT

# GSE98622 のラベルは config から取る（ハードコードのタグ一致は誤判定の温床だった）。
# series matrix で確認した実際の構成:
#   GPL13112: SHAM4h/24h/12m, IRI2h..IRI28d, IRI12m
#   GPL19057: NORM3m/9m/15m, IRI6mN
# → IRI6mN はプラットフォームが違うので、対照をまたいでプールしてはいけない。
IRI_LATE = cfg["external"]["iri_late_labels"]
IRI_CTRL = cfg["external"]["iri_control_labels"]
IRI_LATE_SENS = cfg["external"].get("iri_late_labels_sens", [])
IRI_CTRL_SENS = cfg["external"].get("iri_control_labels_sens", [])


# ------------------------------------------------------------------ helpers
def to_human_space(s: pd.Series, species: str) -> pd.Series:
    """ネコ／マウス／ヒトの系列を共通の *ヒト* シンボル空間（大文字）へ。

    KPMP を同じ行列に載せる以上、共通空間はヒトに揃える必要がある。
    写像は「オルソログ表を引き、無ければ大文字化にフォールバック」。
    片方だけでは取りこぼすため（実測、cat_rna_ctx_late × m_rna_2w）:
        大文字化のみ            n=7356  rho=0.2683
        mouse2human 写像のみ    n=7155  rho=0.2697   ← 335遺伝子を失う
        写像＋fallback（採用）  n=7470  rho=0.2706   ← 最大かつ最良
    ヒト synonym による正規化も試したが n・rho とも改善しなかったため不採用。
    """
    mapping = {"cat": CAT2HUMAN, "mouse": MOUSE2HUMAN}.get(species)
    if mapping:
        s = s.copy()
        s.index = [(mapping.get(g) or str(g)) for g in s.index]
    s = s.copy()
    s = s[pd.notna(s.index)]
    s.index = [str(i).upper() for i in s.index]
    return s[~s.index.duplicated()]


def internal_lfc(name: str) -> pd.Series:
    s = DE.loc[name, "lfc"]
    return s[np.isfinite(s)]


def mouse_iri_lfc(late_labels=None, ctrl_labels=None, tag: str = "late") -> pd.Series | None:
    """GSE98622 の後期IRI vs 対照の log2FC。

    プラットフォーム交絡を避けるため、late/ctrl は config で明示指定する
    （既定は IRI12m vs SHAM12m = 同一プラットフォーム・同週齢）。
    """
    late_labels = late_labels or IRI_LATE
    ctrl_labels = ctrl_labels or IRI_CTRL
    p = INT / "mouse_iri_matrix.parquet"
    if not p.exists():
        log.warning("マウスIRI行列がありません。00b_fetch_external.py を実行してください")
        return None
    mat = pd.read_parquet(p)
    grp = pd.read_parquet(INT / "mouse_iri_grp.parquet").squeeze()
    late = [g for g in grp.unique() if g in late_labels]
    ctrl = [g for g in grp.unique() if g in ctrl_labels]
    if not late or not ctrl:
        log.error("ラベルが見つかりません（late=%s ctrl=%s）。実際のラベル: %s",
                  late_labels, ctrl_labels, sorted(grp.unique()))
        return None
    A = mat[grp[grp.isin(late)].index]
    B = mat[grp[grp.isin(ctrl)].index]
    log.info("IRI[%s] 後期 %s (n=%d) vs 対照 %s (n=%d)",
             tag, late, A.shape[1], ctrl, B.shape[1])
    A, B = np.log2(A + 1), np.log2(B + 1)
    keep = (pd.concat([A, B], axis=1).max(axis=1) > 1)
    res = moderated_ttest(A[keep], B[keep])
    res.to_csv(RES / f"de_mouse_iri_{tag}.csv")
    return res["lfc"]


def kpmp_lfc(compartment: str) -> pd.Series | None:
    """KPMP regional proteomics の CKD vs HRT（pre-computed）を返す。

    DataLake_DEPs.txt は long 形式（Comparison 列でコントラストを持つ）、
    タブ区切り・CRLF 改行。ここで使うのは
    'CKD.vs.HRT.in.TI'（尿細管間質）と 'CKD.vs.HRT.in.G'（糸球体）。

    同一 Gene_name が複数の Accession（アイソフォーム/プロテオグループ）で
    出るため、Adj_pvalue 最小 → 同値なら |LogFC| 最大 の順で1行に代表させる。
    """
    p = EXT / "KPMP" / "DataLake_DEPs.txt"
    if not p.exists():
        log.warning("KPMP proteomics がありません（%s）。手動取得してください", p)
        return None
    df = pd.read_csv(p, sep="\t")
    need = {"Gene_name", "Comparison", "Adj_pvalue", "LogFC"}
    if not need.issubset(df.columns):
        log.error("KPMP の列名が想定と異なります: %s", list(df.columns))
        return None
    target = f"CKD.vs.HRT.in.{compartment}"
    d = df[df.Comparison == target].copy()
    if d.empty:
        log.error("コントラスト %s が見つかりません。存在するのは: %s",
                  target, sorted(df.Comparison.unique()))
        return None
    n_raw = len(d)

    gn = d.Gene_name.astype(str).str.strip()
    blank = d.Gene_name.isna() | (gn == "") | (gn.str.lower() == "nan")
    d, gn = d[~blank], gn[~blank]
    n_blank = int(blank.sum())

    novalue = d.LogFC.isna()
    d, gn = d[~novalue], gn[~novalue]
    n_novalue = int(novalue.sum())

    d = d.assign(_sym=gn.str.upper(), _abs=d.LogFC.abs())
    # Adj_pvalue 昇順 → |LogFC| 降順 で並べ、各シンボルの先頭を代表にする
    d = d.sort_values(["_sym", "Adj_pvalue", "_abs"], ascending=[True, True, False])
    n_before_dedup = len(d)
    d = d.drop_duplicates("_sym", keep="first")
    n_dup_dropped = n_before_dedup - len(d)

    log.info("KPMP %s: %d行 → 有効 %d蛋白（Gene_name欠損 %d行, LogFC欠損 %d行, "
             "重複シンボル統合で %d行を除外）",
             target, n_raw, len(d), n_blank, n_novalue, n_dup_dropped)
    return d.set_index("_sym")["LogFC"]


def collect(design) -> dict[str, pd.Series]:
    out = {}
    for key, species, _entry, src in design:
        if src.startswith("internal:"):
            out[key] = to_human_space(internal_lfc(src.split(":")[1]), species)
        elif src == "external:GSE98622":
            v = mouse_iri_lfc()
            if v is not None:
                out[key] = to_human_space(v, "mouse")
        elif src.startswith("external:KPMP"):
            v = kpmp_lfc(src.split("_")[-1])
            if v is not None:
                out[key] = v
    log.info("利用可能なデータセット: %s", list(out))
    return out


# ------------------------------------------------------------------ core
def distance_matrix(series: dict[str, pd.Series]) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = list(series)
    rho = pd.DataFrame(np.nan, index=keys, columns=keys)
    n = pd.DataFrame(0, index=keys, columns=keys)
    for a, b in combinations(keys, 2):
        x, y = series[a], series[b]
        idx = x.index.intersection(y.index)
        if len(idx) < 500:
            log.warning("%s x %s: 共通遺伝子が %d しかありません", a, b, len(idx))
        r = stats.spearmanr(x.loc[idx], y.loc[idx])[0] if len(idx) > 50 else np.nan
        rho.loc[a, b] = rho.loc[b, a] = r
        n.loc[a, b] = n.loc[b, a] = len(idx)
    np.fill_diagonal(rho.values, 1.0)
    return rho, n


def decompose_rna(rho: pd.DataFrame) -> dict:
    """RNA層：入口の効果と種の効果を、同じ距離尺度（1 - rho）で並べて比較する。

    この2x2はマウスだけが両方の入口（Pod-TRECK=糸球体, IRI=尿細管）を
    持つことで成立する。蛋白層にはマウスIRIが無いので分解できない。
    """
    d = 1 - rho

    def get(a, b):
        return float(d.loc[a, b]) if a in d.index and b in d.columns else np.nan

    out = {}
    # 種を固定して入口を変える
    out["d_entry_mouse (PodTRECK vs IRI)"] = get("mouse_podtreck_rna", "mouse_iri_late_rna")
    # 入口を固定して種を変える
    out["d_species_tubular (cat natural vs mouse IRI)"] = get("cat_natural_ctx_rna", "mouse_iri_late_rna")
    # 交絡したままの従来比較（参考）
    out["d_confounded (cat natural vs PodTRECK)"] = get("cat_natural_ctx_rna", "mouse_podtreck_rna")

    e = out["d_entry_mouse (PodTRECK vs IRI)"]
    sp = out["d_species_tubular (cat natural vs mouse IRI)"]
    if np.isfinite(e) and np.isfinite(sp):
        out["結論"] = ("種差 > 病因差（入口を揃えても種が効く）" if sp > e
                       else "病因差 > 種差（入口を揃えれば種は問題にならない）")
        out["比 (種/病因)"] = round(sp / e, 3)
    return out


def decompose_prot(rho: pd.DataFrame) -> dict:
    """蛋白層：ヒト（KPMP）への距離。

    絶対値ではなく順位で読むこと。KPMP はレーザーマイクロダイセクションで
    TI/G を分離しているのに対し、マウスは全腎、ネコは皮質/髄質の分離であり、
    区画解像度が非対称だから。またネコのプロテオームは蛋白ごとに行中心化
    されているため、順位ベース（Spearman）以外は使えない。
    """
    d = 1 - rho

    def get(a, b):
        return float(d.loc[a, b]) if a in d.index and b in d.columns else np.nan

    out = {}
    out["d_cat_to_human_TI"] = get("cat_natural_ctx_prot", "human_kpmp_ti")
    out["d_catmed_to_human_TI"] = get("cat_natural_med_prot", "human_kpmp_ti")
    out["d_podtreck_to_human_TI"] = get("mouse_podtreck_prot", "human_kpmp_ti")
    out["d_cat_to_human_G"] = get("cat_natural_ctx_prot", "human_kpmp_g")
    out["d_catmed_to_human_G"] = get("cat_natural_med_prot", "human_kpmp_g")
    out["d_podtreck_to_human_G"] = get("mouse_podtreck_prot", "human_kpmp_g")
    out["d_human_TI_vs_G"] = get("human_kpmp_ti", "human_kpmp_g")

    nice = {"cat": "ネコ皮質", "catmed": "ネコ髄質", "podtreck": "Pod-TRECK"}
    for comp in ["TI", "G"]:
        cand = {k: out.get(f"d_{k}_to_human_{comp}") for k in nice}
        cand = {k: v for k, v in cand.items() if v is not None and np.isfinite(v)}
        if cand:
            order = sorted(cand, key=cand.get)
            out[f"ヒト{comp}への近さ順位"] = " < ".join(
                f"{nice[k]}({cand[k]:.3f})" for k in order)
    out["注意"] = ("区画解像度が非対称（KPMP=LMDでTI/G分離, マウス=全腎, "
                   "ネコ=皮質/髄質）。距離の絶対値ではなく順位で解釈すること。")
    return out


def permanova(rho: pd.DataFrame, n_perm: int = 9999, seed: int = 0) -> pd.DataFrame:
    """距離行列に対する種／入口の効果の並べ替え検定（PERMANOVA 相当）。

    n が 6 程度と小さいので検出力は無い。効果量（疑似 F）の記述として扱い、
    p 値は参考値にとどめること。
    """
    d = (1 - rho).to_numpy()
    keys = list(rho.index)
    meta = {k: (sp, en) for k, sp, en, _ in DESIGN}
    rows = []
    rng = np.random.default_rng(seed)
    for factor, pos in [("species", 0), ("entry", 1)]:
        labels = np.array([meta[k][pos] for k in keys])

        def pseudo_f(lab):
            within, between = [], []
            for i, j in combinations(range(len(keys)), 2):
                (within if lab[i] == lab[j] else between).append(d[i, j])
            if not within or not between:
                return np.nan
            return float(np.mean(between) / np.mean(within))

        obs = pseudo_f(labels)
        null = np.array([pseudo_f(rng.permutation(labels)) for _ in range(n_perm)])
        null = null[np.isfinite(null)]
        p = (np.sum(null >= obs) + 1) / (null.size + 1) if null.size else np.nan
        rows.append({"factor": factor, "pseudo_F": obs, "p_perm": p, "n_datasets": len(keys)})
    return pd.DataFrame(rows)


def check_adaptive_immunity() -> pd.DataFrame:
    """最優先チェック：マウスIRI後期で適応免疫が動くか。

    動いていれば「マウスに適応免疫が欠落」という主張は
    「Pod-TRECK 3週が急性すぎる（時間スケールの差）」に修正が必要。
    """
    p = INT / "mouse_iri_matrix.parquet"
    if not p.exists():
        log.warning("IRI 行列がないためスキップ")
        return pd.DataFrame()
    mat = pd.read_parquet(p)
    grp = pd.read_parquet(INT / "mouse_iri_grp.parquet").squeeze()
    mat.index = [str(i).upper() for i in mat.index]
    mat = mat[~mat.index.duplicated()]
    genes = ["CD3D", "CD3E", "CD8A", "CD2", "LCK", "IL7R",
             "MS4A1", "CD79A", "CD79B", "IGHM", "MZB1", "JCHAIN", "CCL19", "PTPRC"]
    rows = []
    for g in genes:
        if g not in mat.index:
            rows.append({"gene": g, "note": "未検出"})
            continue
        rec = {"gene": g}
        for lab in sorted(grp.unique()):
            rec[lab] = round(float(mat.loc[g, grp[grp == lab].index].mean()), 2)
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "mouse_iri_adaptive_immunity.csv", index=False)
    log.warning("適応免疫チェック（この結果次第で主張が変わる）:\n%s", df.to_string(index=False))
    return df


def main():
    # ---------------- RNA 層（2x2 が完成する。種と病因の分解はここ） ------------
    rna = collect(DESIGN_RNA)
    rho_r, n_r = distance_matrix(rna)
    rho_r.to_csv(RES / "design_rho_matrix_rna.csv")
    n_r.to_csv(RES / "design_overlap_matrix_rna.csv")
    log.info("RNA層 rho 行列:\n%s", rho_r.round(3).to_string())
    log.info("RNA層 overlap:\n%s", n_r.to_string())

    dec_r = decompose_rna(rho_r)
    pd.Series(dec_r).to_csv(RES / "design_decomposition_rna.csv")
    log.info("RNA層 分解:\n%s", pd.Series(dec_r).to_string())

    pf = permanova(rho_r, seed=SEED)
    pf.to_csv(RES / "design_permanova_rna.csv", index=False)

    # 感度解析: 後期IRI を 6か月（別プラットフォーム）に差し替えても
    # 分解の大小関係が変わらないかを確認する。
    sens = {}
    if IRI_LATE_SENS and IRI_CTRL_SENS:
        v = mouse_iri_lfc(IRI_LATE_SENS, IRI_CTRL_SENS, tag="late_sens")
        if v is not None:
            s2 = dict(rna)
            s2["mouse_iri_late_rna"] = to_human_space(v, "mouse")
            rho2, _ = distance_matrix(s2)
            rho2.to_csv(RES / "design_rho_matrix_rna_sens.csv")
            sens = decompose_rna(rho2)
            pd.Series(sens).to_csv(RES / "design_decomposition_rna_sens.csv")
            log.info("RNA層 感度解析（IRI 6か月）:\n%s", pd.Series(sens).to_string())

    # ---------------- 蛋白層（ヒトへの距離はここ） --------------------------
    prot = collect(DESIGN_PROT)
    rho_p, n_p = distance_matrix(prot)
    rho_p.to_csv(RES / "design_rho_matrix_prot.csv")
    n_p.to_csv(RES / "design_overlap_matrix_prot.csv")
    log.info("蛋白層 rho 行列:\n%s", rho_p.round(3).to_string())
    log.info("蛋白層 overlap:\n%s", n_p.to_string())

    dec_p = decompose_prot(rho_p)
    pd.Series(dec_p).to_csv(RES / "design_decomposition_prot.csv")
    log.info("蛋白層 ヒトへの距離:\n%s", pd.Series(dec_p).to_string())

    # overlap が薄いペアを明示的に拾って警告する
    thin = []
    for label, nmat in [("RNA", n_r), ("蛋白", n_p)]:
        for a, b in combinations(list(nmat.index), 2):
            if nmat.loc[a, b] < 1000:
                thin.append(f"{label}: {a} x {b} = {int(nmat.loc[a, b])}")
    if thin:
        log.warning("共通遺伝子が少ないペア（<1000）: %s", thin)

    ai = check_adaptive_immunity()

    append_summary("09_dataset_matrix / 種×入口の分解", {
        "【RNA層】データセット": list(rna),
        "IRI後期の定義（主）": f"{IRI_LATE} vs {IRI_CTRL}（同一プラットフォーム GPL13112・同週齢）",
        "IRI後期の定義（感度）": f"{IRI_LATE_SENS} vs {IRI_CTRL_SENS}（同一プラットフォーム GPL19057）",
        "RNA層 rho行列": rho_r.round(3).to_dict(),
        "RNA層 overlap": n_r.to_dict(),
        "RNA層 分解": dec_r,
        "RNA層 分解（感度: IRI 6か月）": sens,
        "PERMANOVA(参考値)": pf.round(3).to_dict("records"),
        "【蛋白層】データセット": list(prot),
        "蛋白層 rho行列": rho_p.round(3).to_dict(),
        "蛋白層 overlap": n_p.to_dict(),
        "蛋白層 ヒトへの距離": dec_p,
        "共通遺伝子が少ないペア(<1000)": thin or "なし",
        "適応免疫チェック": ("実施済み。results/mouse_iri_adaptive_immunity.csv 参照"
                             if len(ai) else "IRIデータ未取得のため未実施"),
        "制約1（層の分離）": ("RNA層と蛋白層は独立した距離行列として扱い、層をまたぐ比較は"
                              "行列に含めていない。測定量が異なるため、種・入口の効果と"
                              "測定層の効果が交絡するため。"),
        "制約2（蛋白層の2x2は不完全）": ("マウスIRIに対応する腎プロテオームが存在しないため、"
                                          "蛋白層では種×入口の2x2が完成しない。したがって "
                                          "d_entry / d_species の分解は RNA層でのみ実施し、"
                                          "蛋白層はヒト（KPMP）への距離の評価に限定している。"),
        "制約3（区画解像度の非対称）": ("KPMP はレーザーマイクロダイセクションで TI/G を分離、"
                                        "マウスは全腎、ネコは皮質/髄質の分離。区画の粒度が"
                                        "揃っていないため、ヒトへの距離は絶対値ではなく"
                                        "順位でのみ解釈すること。"),
        "制約4（ネコ蛋白の行中心化）": ("ネコのプロテオーム値は蛋白ごとに行中心化されている"
                                        "（row mean≈0）。log2FC の絶対値は種間比較に使えないため、"
                                        "蛋白層も Spearman（順位ベース）のみで扱っている。"),
        "警告": ("データセット数が少ないので PERMANOVA の p 値は参考値。"
                 "主結果は d_entry と d_species の直接比較で述べること。"),
    }, cfg)


if __name__ == "__main__":
    main()
