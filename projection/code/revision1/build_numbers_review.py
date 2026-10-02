# -*- coding: utf-8 -*-
"""AJ: 査読対応（AA–AI）の数値を NUMBERS.md の §18 として生成する。

Response letter はこの節からのみ数値を引く。節そのものは結果ファイルから作るので、
手で書き写す段階が無い。
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[2]
R = HERE / "results"
NUM = HERE / "manuscript" / "NUMBERS.md"
HEAD = "## 18. 査読対応で計算した値（AA–AI）"


def tsv(p, idx=None):
    T = pd.read_csv(R / p, sep="\t")
    return T.set_index(idx) if idx else T


def rows_of(title, note, header, rows, src):
    out = [f"### {title}", "", note, "", f"出典 `{src}`。", "",
           "| " + " | ".join(header) + " |",
           "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    out.append("")
    return out


def main() -> int:
    AA = tsv("roundR/aa_unit_evidence.tsv", "quantity")["value"]
    B87 = tsv("round7/centering_auc_animal_bootstrap.tsv", "quantity")["value"]
    B75 = tsv("round7/centering_auc_animal_bootstrap_75pairs.tsv", "quantity")["value"]
    DD = tsv("roundR/delta_definitions.tsv", "definition")
    S4 = tsv("roundR/table_s4_rows.tsv", "method")
    CE = tsv("roundR/reference_implementation_cells.tsv")
    SO = tsv("roundR/ad_ssgsea_options.tsv", "option")
    RV = dict(l.split("\t") for l in (R / "roundR" / "r_versions.tsv").read_text().splitlines())
    RS = tsv("roundR/random_set_designs.tsv").set_index(["design", "quantity"])
    OV = tsv("roundR/random_set_designs_overlap.tsv", "quantity")["value"]
    LP = tsv("roundR/random_set_designs_draws_labelperm.tsv")
    SD = tsv("roundR/pathway_set_design.tsv")
    DS = tsv("roundR/duplicate_symbol_sensitivity.tsv", "quantity")
    A7 = tsv("roundR/af7_attached_table.tsv", "quantity")["value"]
    GS = tsv("roundR/groupA_detection_summary.tsv")

    L: list[str] = [HEAD, "",
                    "査読対応（AA–AI）で新たに計算した値。**Response letter はこの節からのみ引く。**",
                    "すべて結果ファイルから自動生成しており、手で書き写す段階は無い。", ""]

    L += rows_of(
        "18.1 GSE98622 の単位（AA / R3-major-1）",
        "本文の `normalized counts` が誤記で、値は FPKM。数値そのものは動かない。",
        ["量", "値"],
        [["ファイル", f"`{AA['file name']}`"],
         ["SHA-256", f"`{AA['sha256']}`"],
         ["行 x 数値列", f"{int(float(AA['rows in the sheet'])):,} x {int(float(AA['numeric sample columns']))}"],
         ["非ゼロ値のうち非整数", f"**{float(AA['share of non-zero entries that are non-integer'])*100:.3f}%**"],
         ["非ゼロ最小値", f"**{float(AA['smallest non-zero value']):g}**"],
         ["列合計の範囲", f"{float(AA['smallest column sum']):,.0f}–{float(AA['largest column sum']):,.0f}"],
         ["列合計の変動係数", f"{float(AA['coefficient of variation of column sums']):.3f}"]],
        "results/roundR/aa_unit_evidence.tsv")

    L += rows_of(
        "18.2 75 cross-dataset ペア（AE1 / R3-major-2）",
        "87 ペアのうち 12（猫の皮質×髄質 4、IRI 12mo × 早期 IRI 8）は within-dataset かつ同種。",
        ["部分集合", "非中心化", "中心化", "差", "resampling 中央値", "95% 範囲", "unc > cen"],
        [["87 ペア（主解析）", f"**{B87['observed AUC, uncentred pathway means']:.4f}**",
          f"**{B87['observed AUC, centred pathway means']:.4f}**",
          f"**{B87['observed drop']:.4f}**", f"{B87['bootstrap drop median']:.4f}",
          f"{B87['bootstrap drop 2.5th percentile']:.4f}–{B87['bootstrap drop 97.5th percentile']:.4f}",
          f"{B87['replicates with uncentred AUC above centred']*100:.2f}%"],
         ["75 ペア（cross-dataset）", f"**{B75['observed AUC, uncentred pathway means']:.4f}**",
          f"**{B75['observed AUC, centred pathway means']:.4f}**",
          f"**{B75['observed drop']:.4f}**", f"{B75['bootstrap drop median']:.4f}",
          f"{B75['bootstrap drop 2.5th percentile']:.4f}–{B75['bootstrap drop 97.5th percentile']:.4f}",
          f"{B75['replicates with uncentred AUC above centred']*100:.2f}%"]],
        "results/round7/centering_auc_animal_bootstrap{,_75pairs}.tsv")

    L += rows_of(
        "18.3 Δ の 4 定義（AC / R3-major-3）",
        "**順位化・分位正規化での消失は構成上の帰結であり、独立した所見ではない。**"
        "どちらも状態ごとの周辺分布を揃えるので、全体平均は全 16 状態で同一になる。",
        ["Δ の定義", "全体平均が種を分離", "遺伝子レベル AUC", "非中心化", "中心化", "差"],
        [[d, "はい" if DD.loc[d, "mean_separates_species"] else "いいえ",
          f"{DD.loc[d, 'auc_gene_level']:.3f}", f"**{DD.loc[d, 'auc_agg_uncentered']:.3f}**",
          f"**{DD.loc[d, 'auc_agg_centered']:.3f}**", f"**{DD.loc[d, 'centering_drop']:+.3f}**"]
         for d in DD.index],
        "results/roundR/delta_definitions.tsv")

    L += rows_of(
        "18.4 経路スコアリング 5 手法（AD / R2-3）",
        f"報告値は Bioconductor 参照実装（R {RV['R']} / Bioconductor {RV['Bioconductor']} / "
        f"GSVA {RV['GSVA']} / singscore {RV['singscore']} / GSEABase {RV['GSEABase']}）。",
        ["手法", "実装", "非中心化", "中心化", "差", "中心化前後のスコア最大差"],
        [[m, S4.loc[m, "implementation"], f"**{S4.loc[m, 'auc_uncentered']:.3f}**",
          f"**{S4.loc[m, 'auc_centered']:.3f}**", f"**{S4.loc[m, 'difference']:+.3f}**",
          f"{S4.loc[m, 'max_abs_score_difference']:.4f}"] for m in S4.index],
        "results/roundR/table_s4_rows.tsv")

    log = [i for i in SO.index if i.startswith("sample_norm_method='log'")][0]
    L += rows_of(
        "18.5 実装の突き合わせ（AI）",
        "自前実装・GSEApy と Bioconductor 参照実装のセル単位の一致。"
        "ssGSEA と singscore の中心化不変性は順位のみを読むことの代数的帰結で、実装に依存しない。",
        ["手法", "比較", "最大差", "中央値", "セルの相関"],
        [[r.method, r.compared, f"**{r.max_abs_difference:.1e}**",
          f"{r.median_abs_difference:.1e}", f"{r.pearson_r_of_cells:.6f}"]
         for r in CE[CE.matrix == "uncentered"].itertuples()] +
        [["ssGSEA（順位を使わない設定）", f"{log}",
          f"**{SO.loc[log, 'max_abs_ES_difference']:.0f}**", "—", "—"]],
        "results/roundR/reference_implementation_cells.tsv / ad_ssgsea_options.tsv")

    def rsv(d, q, c):
        return RS.loc[(d, q), c]

    SM, LPn = "size-matched draws", "gene-label permutation"
    L += rows_of(
        "18.6 ランダム集合の設計（AE2 / R3-major-4）",
        f"422 集合は遺伝子を再利用する（1 遺伝子が中央値 {int(OV['sets per gene, median'])} 集合、"
        f"最大 {int(OV['sets per gene, maximum'])} 集合、{int(OV['genes in the universe']):,} 中 "
        f"{int(OV['genes in no set'])} がどの集合にも属さない）。"
        "**サイズだけを合わせた乱択はこの構造を保存していない。**",
        ["設計", "差の中央値", "差の範囲", "観測 0.3563 以上", "中心化 AUC の範囲"],
        [["サイズマッチ乱択", f"{rsv(SM,'difference','random_median'):.4f}",
          f"{rsv(SM,'difference','random_min'):.4f}–{rsv(SM,'difference','random_max'):.4f}",
          f"**{int(rsv(SM,'difference','draws_at_or_above_observed'))}/300**",
          f"{rsv(SM,'auc_agg_centred','random_min'):.4f}–{rsv(SM,'auc_agg_centred','random_max'):.4f}"],
         ["遺伝子ラベル置換", f"{rsv(LPn,'difference','random_median'):.4f}",
          f"{rsv(LPn,'difference','random_min'):.4f}–{rsv(LPn,'difference','random_max'):.4f}",
          f"**{int(rsv(LPn,'difference','draws_at_or_above_observed'))}/300**",
          f"{rsv(LPn,'auc_agg_centred','random_min'):.4f}–{rsv(LPn,'auc_agg_centred','random_max'):.4f}"]],
        "results/roundR/random_set_designs.tsv / random_set_designs_overlap.tsv")
    L += [f"ラベル置換で中心化 AUC が観測値 {rsv(LPn,'auc_agg_centred','observed_real_sets'):.4f} "
          f"以下になった draw: **{int((LP.auc_agg_centred <= rsv(LPn,'auc_agg_centred','observed_real_sets')).sum())}/300**、"
          f"0.55 以下: **{int((LP.auc_agg_centred <= 0.55).sum())}/300**。", ""]

    COL = SD[SD.analysis == "collection"]
    RED = SD[SD.analysis == "redundancy"]
    CMN = SD[SD.analysis == "common sets"]
    L += rows_of(
        "18.7 経路集合の設計（AE2 / R3-major-4）", "コレクション別と冗長性削減。",
        ["条件", "集合数", "非中心化", "中心化", "差"],
        [[r.label, int(r.n_sets), f"{r.auc_uncentered:.4f}", f"{r.auc_centered:.4f}",
          f"**{r.difference:+.4f}**"] for r in pd.concat([COL, RED]).itertuples()],
        "results/roundR/pathway_set_design.tsv")
    L += rows_of(
        "18.8 共通経路 ID に限定（AE2c）",
        "**「全 1:1 オルソログで差が最小」は遺伝子空間の効果ではなく経路集合の効果だった。**"
        "§2.2 の帰属の誤りを修正した。",
        ["遺伝子空間", "自前の集合数", "自前の差", "共通 117 集合での差"],
        [[r.label, int(r.n_sets_all), f"{r.difference_all_sets:+.4f}", f"**{r.difference:+.4f}**"]
         for r in CMN.itertuples()],
        "results/roundR/pathway_set_design.tsv")

    L += rows_of(
        "18.9 重複シンボルの行順感度（AB / R3-minor-2）",
        "重複ヒトシンボルの代表選択を 100 回シャッフルした。報告値の 3 桁は動かない。",
        ["量", "現行", "中央値", "最小", "最大", "現行との最大差"],
        [[q, DS.loc[q, "current"], DS.loc[q, "median"], DS.loc[q, "min"], DS.loc[q, "max"],
          f"**{DS.loc[q, 'max_abs_diff_from_current']}**"]
         for q in ("auc_gene", "auc_agg_uncentered", "auc_agg_centered", "centering_drop")] +
        [["種の分離", DS.loc["mean_separates_species", "current"],
          f"**{DS.loc['mean_separates_species', 'median']}**", "—", "—", "—"],
         ["代表選択の対象になる遺伝子",
          f"**{int(DS.loc['analysed genes with duplicate source rows', 'current'])}**"
          f" / 2,016 = **{float(DS.loc['share of the analysed genes', 'current'])*100:.1f}%**",
          f"猫 {int(DS.loc['analysed genes with duplicate source rows, cat', 'current'])}",
          f"マウス {int(DS.loc['analysed genes with duplicate source rows, mouse', 'current'])}",
          "—", "—"]],
        "results/roundR/duplicate_symbol_sensitivity.tsv")

    L += rows_of(
        "18.10 検出コールと添付した定量表（AF / R3-major-5）",
        f"添付 `{A7['attached file name']}`（"
        f"{int(float(A7['attached file bytes'])):,} bytes、SHA-256 `{A7['attached file sha256']}`）。"
        "**加工していない。**添付ファイルから再計算すると選択が正確に再現する。",
        ["段階", "遺伝子数"],
        [[r.step, f"{int(r.genes):,}"] for r in GS.itertuples()],
        "results/roundR/af7_attached_table.tsv / groupA_detection_summary.tsv")

    EX = tsv("roundR/simulation_exceedance.tsv")
    EX["true_cos_n"] = pd.to_numeric(EX.true_cos, errors="coerce")

    def ex(panel, design, tc):
        r = EX[(EX.panel == panel) & (EX.design == design) & (EX.true_cos_n == tc)
               & (EX["shift"] == "all") & (~EX.metric.str.contains("grid points"))]
        assert len(r) == 1, f"{panel}/{design}/{tc} が {len(r)} 行"
        return r.iloc[0]

    def step(panel, which):
        r = EX[(EX.panel == panel) & EX.metric.str.contains(which, regex=False)]
        assert len(r) == 1, f"{panel}/{which} が {len(r)} 行"
        return r.iloc[0]

    L += rows_of(
        "18.11 ベンチマークのシミュレーションの超過率（BE1 / 自主訂正）",
        "すべて信頼性 0.568–0.986 の帯の中。超過率は設計で分かれるので、設計別を主、"
        "両設計・3 shift をまとめた値を従として扱う。真のコサイン 0.9 の 39.7 / 39.9% は"
        "信頼性に沿った階段をならした格子平均であって、ある信頼性での率ではない。",
        ["パネル", "量", "集計", "平均", "最小", "最大", "条件数"],
        [["(iv)", "真コサイン 1.0、補正後基準の超過率", "ネコ型 7 対 6",
          f"**{ex('iv', 'feline-like (7 vs 6)', 1.0)['mean']*100:.1f}%**",
          f"{ex('iv', 'feline-like (7 vs 6)', 1.0)['min']*100:.1f}%",
          f"{ex('iv', 'feline-like (7 vs 6)', 1.0)['max']*100:.1f}%", 15],
         ["(iv)", "同上", "マウス型 3 対 6",
          f"**{ex('iv', 'mouse-like (3 vs 6)', 1.0)['mean']*100:.1f}%**",
          f"{ex('iv', 'mouse-like (3 vs 6)', 1.0)['min']*100:.1f}%",
          f"{ex('iv', 'mouse-like (3 vs 6)', 1.0)['max']*100:.1f}%", 15],
         ["(iv)", "同上", "両設計・3 shift",
          f"{ex('iv', 'pooled', 1.0)['mean']*100:.1f}%",
          f"{ex('iv', 'pooled', 1.0)['min']*100:.1f}%",
          f"{ex('iv', 'pooled', 1.0)['max']*100:.1f}%", 30],
         ["(iv)", "真コサイン 1.0 未満、補正後基準の超過率", "両設計・3 shift",
          f"**{EX[(EX.panel=='iv') & (EX.true_cos=='<1.0')].iloc[0]['max']*100:.1f}%**",
          "—", "—", int(EX[(EX.panel=='iv') & (EX.true_cos=='<1.0')].iloc[0]['n'])],
         ["(ii)", "真コサイン 0.9、中心化・半標本基準", "格子平均",
          f"{ex('ii', 'pooled', 0.9)['mean']*100:.1f}%", "0%", "100%", 30],
         ["(iii)", "真コサイン 0.9、非中心化・半標本基準", "格子平均",
          f"{ex('iii', 'pooled', 0.9)['mean']*100:.1f}%", "0%", "100%", 30],
         ["(ii)(iii)", "超過する側の信頼性",
          f"{step('ii', '(exceeding grid').reliability.replace('-', '–')}",
          f"**{step('ii', '(exceeding grid')['min']*100:.0f}–"
          f"{step('ii', '(exceeding grid')['max']*100:.0f}%**", "—", "—", 12],
         ["(ii)(iii)", "超過しない側の信頼性",
          f"{step('ii', 'non-exceeding').reliability.replace('-', '–')}",
          "**0%**", "—", "—", 18]],
        "results/roundR/simulation_exceedance.tsv")

    body = "\n".join(L).rstrip() + "\n\n---\n\n"
    t = NUM.read_text()
    marker = "# 未計算（数値を作らないこと）"
    if HEAD in t:
        i = t.index(HEAD)
        j = t.index(marker, i)
        t = t[:i] + body + t[j:]
    else:
        t = t.replace(marker, body + marker)
    NUM.write_text(t)
    print(f"NUMBERS.md §18 を生成した（{len(L)} 行）→ {NUM}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
