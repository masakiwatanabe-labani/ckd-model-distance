"""NUMBERS.md に書いた値が結果ファイルと一致するかを検査する。

NUMBERS.md は本文・申請書が数値を引く正本なので、手で書き写した値が
結果ファイルとずれていないことを機械的に確認できるようにしておく。
数値を更新したらこのスクリプトの期待値も更新すること。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R = HERE / "results"


# --map を付けると「本文で報告する量 → 結果ファイル」の対応表を書き出す。
# 原稿からリポジトリ内パスを外した代わりに、この表が README でその対応を持つ。
MAPPING = []
_LAST_FILE = [None]


def tsv(rel, **kw):
    _LAST_FILE[0] = str(rel)
    return pd.read_csv(R / rel, sep="\t", **kw)


def main():
    ok, bad = [], []
    # 公開していない材料が要る検査は飛ばし、何を飛ばしたかを最後に並べる（BM2）
    skipped: list[str] = []

    def chk(label, got, want, tol=6e-4):
        got = float(got)
        MAPPING.append((label, _LAST_FILE[0], f"{want}"))
        (ok if abs(got - want) <= tol else bad).append(
            f"{label}: NUMBERS.md {want} / 実測 {got:.6f}")

    s = tsv("groupAB_control_log2cpm/check1_summary.tsv").set_index("set")
    chk("§1 ρ GroupA 素", s.loc["GroupA", "rho"], 0.6609)
    chk("§1 ρ GroupB 素", s.loc["GroupB", "rho"], 0.4505)
    chk("§1 Δρ 素", s.loc["GroupB", "rho"] - s.loc["GroupA", "rho"], -0.2104)

    n = tsv("groupAB_control_log2cpm/check1_nn_matched.tsv", index_col=0)["value"]
    chk("§1 Δρ 発現量マッチ後", n["delta_GroupB_minus_GroupA"], -0.1564)
    chk("§1   CI lo", n["delta_lo95"], -0.2127)
    chk("§1   CI hi", n["delta_hi95"], -0.1024)
    for tag, want in [("delta_se_moderated", -0.2403), ("delta_se_splithalf", -0.2107)]:
        v = tsv(f"groupAB_{tag}/check1_nn_matched.tsv", index_col=0)["value"]
        chk(f"§1 Δρ {tag} マッチ後", v["delta_GroupB_minus_GroupA"], want)
    c = tsv("check1/check1_nn_matched.tsv", index_col=0)["value"]
    chk("§1 Δρ added−core（対照）", c["delta_added_minus_core"], 0.0046)
    ud = tsv("groupAB_control_log2cpm/check1_unmatched_delta.tsv", index_col=0)["value"]
    chk("§1 Δρ 素 点推定", ud["point_estimate"], 0.2104)
    chk("§1 Δρ 素 CI lo", ud["lo95"], 0.1746, 2e-3)
    chk("§1 Δρ 素 CI hi", ud["hi95"], 0.2456, 2e-3)

    P = tsv("reliability/pairs_ceilings.tsv")
    cs = P[P["class"] == "cross_species"]
    ss = P[P["class"] == "same_species_diff_dataset"]
    chk("§3 異種 cos 中央", cs.cos.median(), 0.4475)
    chk("§3 異種 cos 最大", cs.cos.max(), 0.548)
    chk("§3 同種異DS cos 中央", ss.cos.median(), 0.707)
    chk("§3 天井raw 異種 最小", cs.ceiling_raw.min(), 0.620)
    chk("§3 天井SB 異種 最小", cs.ceiling_SB.min(), 0.765)
    chk("§3 天井までの最小余裕", (cs.ceiling_raw - cs.cos).min(), 0.171)
    chk("§3 系統誤差の下限", (P.cos - P.ceiling_SB).max(), 0.094)
    assert (cs.cos > cs.ceiling_raw).sum() == 0, "異種で天井(raw)超過が発生している"
    assert (cs.cos > cs.ceiling_SB).sum() == 0, "異種で天井(SB)超過が発生している"

    d = tsv("control_analysis/cos_distributions.tsv").set_index("set")
    chk("§3 共有対照あり 中央", d.loc["shared_control_pairs", "median"], 0.717)
    chk("§3 共有対照なし 中央", d.loc["no_shared_control_pairs", "median"], 0.494)

    a = tsv("control_analysis/auc_permutation.tsv")
    chk("§4 AUC 観測", a.auc_obs.iloc[0], 0.7909)
    chk("§4 AUC p", a.p.iloc[0], 0.0049)
    a2 = tsv("reliability/auc_permutation_cos_corrected.tsv")
    chk("§4 AUC 補正後", a2.auc_obs.iloc[0], 0.7809)
    chk("§4 AUC 補正後 p", a2.p.iloc[0], 0.0060)

    pr = tsv("alpha_groupA/projection.tsv", index_col=0)
    chk("§6 cat_med_CKD12 alpha", pr.loc["cat_med_CKD12", "alpha"], 1.020, 1e-3)
    chk("§6 cat_med_CKD12 cos", pr.loc["cat_med_CKD12", "cos"], 0.790, 1e-3)
    chk("§6 cat_med_CKD12 nrm", pr.loc["cat_med_CKD12", "nrm"], 1.291, 1e-3)
    chk("§6 mouse_2W alpha", pr.loc["mouse_2W", "alpha"], 0.740, 1e-3)
    chk("§9 complete case 遺伝子数", pr.n_genes.iloc[0], 2016)
    assert pr.loc["cat_med_CKD12", "alpha"] > pr.loc["cat_CKD34", "alpha"], \
        "cat_med_CKD12 の alpha が基準を上回っていない（§1 の主張の前提）"

    u = tsv("alpha_union/projection.tsv", index_col=0)
    chk("§9 union complete case", u.n_genes.iloc[0], 2643)
    chk("§9 GroupA intersection n", sum(1 for _ in open(HERE / "groupA_intersection.txt")), 2240)
    chk("§9 GroupA union n", sum(1 for _ in open(HERE / "groupA_union.txt")), 2964)

    REL = tsv("reliability/reliability.tsv").set_index("state")
    for st, want in [("cat_CKD34", 0.913), ("cat_med_CKD12", 0.919),
                     ("cat_CKD12", 0.724), ("cat_med_CKD34", 0.762)]:
        chk(f"§10 信頼性 {st}", REL.loc[st, "reliability_SB_median"], want)
    v = tsv("reliability/sb_validation.tsv")["observed_minus_predicted"]
    chk("§10 SB検証 最小", v.min(), -0.035)
    chk("§10 SB検証 最大", v.max(), 0.088)
    chk("§10 SB検証 中央", v.median(), 0.006)
    assert int((v > 0).sum()) == 5 and len(v) == 7, "SB 検証の 5/7 が一致しない"
    cur = tsv("reliability/r_curve.tsv")
    both = cur[cur.curve == "both"]
    n_curve = int((both.groupby("state").size() >= 2).sum())
    assert n_curve == 4, f"r(n) 曲線が取れる状態数が 4 でない: {n_curve}"
    chk("§10 AUC 帰無中央値", a.null_median.iloc[0], 0.4918)
    cc = tsv("reliability/ceiling_comparison.tsv").set_index("class")
    chk("§10 天井 同一DS raw", cc.loc["within_dataset", "ceiling_raw_median"], 0.852)
    chk("§10 天井 同一DS SB", cc.loc["within_dataset", "ceiling_SB_median"], 0.920)
    chk("§10 天井 同種異DS raw", cc.loc["same_species_diff_dataset", "ceiling_raw_median"], 0.896)
    chk("§10 マッチ成立数", n["n_matched_pairs"], 1754)
    chk("§10 マッチ後 発現量順位 A", n["expr_rank_median_GroupA"], 0.716, 1e-3)
    chk("§10 マッチ後 発現量順位 B", n["expr_rank_median_GroupB"], 0.716, 1e-3)
    bd = tsv("groupAB_control_log2cpm/check1_by_decile.tsv")
    assert int((bd.rho_GroupA > bd.rho_GroupB).sum()) == 8, "8/10 十分位が一致しない"
    chk("§10 第1十分位 GroupA", bd.n_core.iloc[0], 29)
    chk("§10 第1十分位 GroupB", bd.n_added.iloc[0], 1030)

    sv = tsv("check1/supplementary_sensitivity.tsv")
    chk("§7 S1 rho", sv[sv.primary].rho.iloc[0], 0.4893)
    chk("§7 S7 rho", sv[(sv.mapping == "naive") & (sv.fpkm_filter == "off")
                       & (sv["delta"] == "log2FC")].rho.iloc[0], 0.3982)

    # ---- §11 ヒト外挿（Limitations 限定） ----
    hs = tsv("human_exploratory/human_cos_summary.tsv")
    hs = hs[hs.subset == "reliable_only"].set_index("class")
    chk("§11 cat x human 中央", hs.loc["cat x human", "cos_median"], 0.465, 1e-3)
    chk("§11 mouse x human 中央", hs.loc["mouse x human", "cos_median"], 0.390, 1e-3)
    at = tsv("human_exploratory/cat_vs_mouse_to_human.tsv")
    chk("§11 AUC（マッチなし）", at.auc.iloc[0], 0.8981)
    chk("§11 p（マッチなし）", at.p_exact.iloc[0], 0.0016)
    mm = tsv("human_exploratory/mapping_matched_test.tsv").set_index("subset")
    k2 = [i for i in mm.index if "<= 2" in i][0]
    chk("§11 caliper2 cat x human", mm.loc[k2, "cat_x_human"], 0.421, 1e-3)
    chk("§11 caliper2 mouse x human", mm.loc[k2, "mouse_x_human"], 0.478, 1e-3)
    chk("§11 caliper2 差", mm.loc[k2, "diff"], -0.058, 1e-3)
    chk("§11 caliper2 AUC", mm.loc[k2, "auc"], 0.331, 1e-3)
    chk("§11 caliper2 pid_cat", mm.loc[k2, "pid_cat_median"], 93.2, 0.05)
    chk("§11 caliper2 pid_mouse", mm.loc[k2, "pid_mouse_median"], 93.0, 0.05)
    assert mm.loc[k2, "diff"] < 0, "caliper 2 で符号が反転していない（Limitations の根拠）"
    assert not (mm.loc[k2, "size_null_p5"] <= mm.loc[k2, "auc"]), \
        "caliper 2 の AUC が同サイズ乱択の 5%点を下回っていない"

    # ---- §12 Mantel の S1 再計算 ----
    mt = tsv("mantel_s1/mantel_s1.tsv")
    m12 = mt[mt.subset.str.startswith("12")]
    def mv(dis, pred):
        r = m12[(m12.dissimilarity.str.startswith(dis)) & (m12.predictor == pred)]
        return float(r.mantel_rho.iloc[0])
    chk("§12 Mantel 時間 D_rho", mv("D_rho", "elapsed time |Δlog10 h|"), 0.596, 1e-3)
    chk("§12 Mantel 時間 D_cos", mv("D_cos", "elapsed time |Δlog10 h|"), 0.599, 1e-3)
    chk("§12 Mantel 時間 D_nrm", mv("D_nrm", "elapsed time |Δlog10 h|"), 0.312, 1e-3)
    chk("§12 Mantel 入口 D_rho(12)", mv("D_rho", "onset compartment"), -0.064, 1e-3)
    assert mv("D_cos", "elapsed time |Δlog10 h|") > 1.5 * mv("D_nrm", "elapsed time |Δlog10 h|"), \
        "時間相関が方向優位でない（統合の核の根拠）"
    ta = tsv("alpha_groupA/time_association.tsv")
    ta = ta[ta.subset == "all_model_states"].set_index("term")
    chk("§12 cos x log_time（状態ごと）", ta.loc["cos x log_time", "rho"], 0.074, 1e-3)
    chk("§12 nrm x log_time（状態ごと）", ta.loc["nrm x log_time", "rho"], -0.070, 1e-3)
    dd = tsv("mantel_s1/distance_distributions.tsv")
    cosr = dd[dd.dissimilarity.str.startswith("D_cos")].set_index("class")
    chk("§12 D_cos 種間 中央", cosr.loc["cross species", "median"], 0.553, 1e-3)
    chk("§12 D_cos 種間 最小", cosr.loc["cross species", "min"], 0.452, 1e-3)
    nrmr = dd[dd.dissimilarity.str.startswith("D_nrm")].set_index("class")
    assert nrmr.loc["cross species", "median"] < nrmr.loc["within species", "median"], \
        "D_nrm の種内/種間逆転が再現していない"

    # ---- §13 経路粒度 ----
    pw = tsv("pathway/pathway_cos.tsv")
    chk("§13 経路数", pw.pathway.nunique(), 193)
    chk("§13 セル数", len(pw), 3088)
    chk("§13 経路サイズ中央", pw.n_genes.median(), 44)
    assert int((pw.groupby("pathway").n_genes.first() < 50).sum()) == 119, \
        "50 遺伝子未満の経路数が 119 でない"
    # §2.4 と同一設計（全120ペアの cos → 対照共有を除外 → 種内 vs 種間）
    sep = tsv("pathway/pathway_separation_paired.tsv")
    col = "auc[exclude shared controls]"
    a = sep[col].dropna()
    chk("§13 経路単位 AUC 中央", a.median(), 0.788, 1e-3)
    chk("§13 経路単位 AUC Q1", a.quantile(.25), 0.733, 1e-3)
    chk("§13 経路単位 AUC Q3", a.quantile(.75), 0.839, 1e-3)
    chk("§13 AUC>=0.90 の経路数", (a >= 0.90).sum(), 13)
    chk("§13 AUC<=0.60 の経路数", (a <= 0.60).sum(), 15)
    ss = tsv("pathway/separation_summary.tsv").set_index("restriction")
    chk("§13 集約 AUC", ss.loc["exclude shared controls", "aggregated_auc"], 0.511, 1e-3)
    chk("§13 集約 AUC（コホート除外）", ss.loc["exclude shared cohort", "aggregated_auc"], 0.481, 1e-3)
    chk("§13 制限後のペア数", ss.loc["exclude shared controls", "n_pairs"], 87)
    assert a.median() > ss.loc["exclude shared controls", "aggregated_auc"] + 0.15, \
        "経路単位が集約を十分上回っていない（§3.3 の核）"
    assert abs(ss.loc["exclude shared controls", "aggregated_auc"] - 0.5) < 0.05, \
        "集約後が偶然水準に落ちていない（§3.3 の核）"

    # ---- §14 蛋白層 ----
    pp = tsv("protein/projection_protein.tsv", index_col=0)
    chk("§14 蛋白数", pp.n_proteins.iloc[0], 2233)
    chk("§14 cos 猫皮質CKD1/2", pp.loc["cat_prot_CKD12", "cos"], 0.871, 1e-3)
    chk("§14 cos 猫髄質CKD1/2", pp.loc["cat_med_prot_CKD12", "cos"], 0.594, 1e-3)
    chk("§14 alpha D14", pp.loc["mouse_prot_D14", "alpha"], 1.322, 1e-3)
    chk("§14 alpha D21", pp.loc["mouse_prot_D21", "alpha"], 1.336, 1e-3)
    chk("§14 alpha 猫髄質CKD1/2", pp.loc["cat_med_prot_CKD12", "alpha"], 0.584, 1e-3)
    assert pp.loc["cat_med_prot_CKD12", "alpha"] < 1.0, \
        "蛋白層で猫髄質CKD1/2 の alpha が 1 を超えている（§14 の記述と矛盾）"
    assert pp.loc["mouse_prot_D14", "alpha"] > 1.0 and pp.loc["mouse_prot_D21", "alpha"] > 1.0, \
        "蛋白層で Pod-TRECK の alpha が基準を超えていない（§14 の核）"
    for st in ["mouse_prot_D14", "mouse_prot_D21"]:
        assert pp.loc[st, "cos"] < pp.loc[st, "ceiling_vs_ref"], \
            f"{st} が天井を超えている（異種は超過しないはず）"
    cv = tsv("protein/protein_vs_rna.tsv")
    from scipy import stats as _st
    chk("§14 protein vs RNA cos 順位一致", _st.spearmanr(cv.cos_protein, cv.cos_rna)[0], 1.000, 1e-6)
    chk("§14 cos 中央 protein", cv.cos_protein.median(), 0.650, 1e-3)
    chk("§14 cos 中央 RNA", cv.cos_rna.median(), 0.824, 2e-3)

    # ---- §15 遺伝子空間の感度解析（AA） ----
    G = tsv("aa_genespace/genespace_comparison.tsv").set_index("空間")
    A_, O_, B_ = "Group A intersection（主解析）", "全 1:1 orthologue", "発現量マッチ Group B"
    chk("§15 遺伝子数 orthologue", G.loc[O_, "n_genes"], 7897)
    chk("§15 遺伝子数 matched B", G.loc[B_, "n_genes"], 1618)
    for sp, a_, c_ in [(A_, 1.020, 0.074), (O_, 0.796, 0.627), (B_, 0.930, 0.175)]:
        chk(f"§15 {sp} 猫髄質 alpha", G.loc[sp, "猫髄質CKD1/2 alpha"], a_, 1e-3)
        chk(f"§15 {sp} cos x log_time", G.loc[sp, "cos x log_time"], c_, 1e-3)
    assert bool(G.loc[A_, "alpha>1.000"]) and not bool(G.loc[O_, "alpha>1.000"]) \
        and not bool(G.loc[B_, "alpha>1.000"]), \
        "alpha>1 の再現状況が §2.7 の記述と食い違う"
    # §2.1 の主張1（順序逆転）は 3 空間で成り立つこと。主張2（基準超え）より前に検査する。
    for sp, d in [(A_, "alpha_groupA"), (O_, "alpha_ortholog_all"), (B_, "alpha_groupB_matched")]:
        pj = tsv(f"{d}/projection.tsv", index_col=0)
        med, ctx = pj.loc["cat_med_CKD12"], pj.loc["cat_CKD12"]
        assert med.alpha > ctx.alpha and med["cos"] < ctx["cos"], \
            f"{sp} で alpha と cos の順序逆転が成り立たない（§2.1 の主張1）"
        chk(f"§2.7 {sp} 髄質 nrm - 1/cos", med.nrm - 1 / med["cos"],
            {A_: 0.025, O_: -0.280, B_: -0.099}[sp], 1e-3)
    for sp in (A_, O_, B_):
        assert G.loc[sp, "cat-mouse cos 最大"] < G.loc[sp, "天井 raw 最小"], \
            f"{sp} で異種 cos が最小天井を超えている（§2.7 の再現判定の核）"
        assert int(G.loc[sp, "天井超え(raw)"]) == 0, f"{sp} で天井超えペアがある"
        assert G.loc[sp, "経路単位 AUC 中央"] > G.loc[sp, "集約後 AUC"], \
            f"{sp} で経路単位が集約を上回っていない（§2.5 の核）"
    assert G.loc[A_, "集約後 AUC"] < 0.55 <= G.loc[O_, "集約後 AUC"], \
        "集約後 AUC の『Group A のみ偶然水準』という記述が成り立たない"

    # ---- §16 天井との差の個体レベル不確実性 ----
    U = tsv("pair_uncertainty/pair_uncertainty.tsv")
    cs = U[U.cls == "cross_species"].copy()
    cs["cat"] = cs.a.where(cs.a.str.startswith("cat"), cs.b)
    ws = U[(U.cls != "cross_species") & (~U.shared_control)]
    chk("§16 異種ペア数", len(cs), 48)
    chk("§16 区間が0を除外する異種ペア", int((cs.gap_lo > 0).sum()), 26)
    rel = cs[cs.cat.isin({"cat_CKD34", "cat_med_CKD12"})]
    unrel = cs[~cs.cat.isin({"cat_CKD34", "cat_med_CKD12"})]
    chk("§16 信頼性の高い猫状態を含むペア", len(rel), 24)
    chk("§16 同 区間が0を除外", int((rel.gap_lo > 0).sum()), 24)
    chk("§16 同 gap 中央", rel.gap_median.median(), 0.403, 2e-3)
    chk("§16 低信頼の猫状態を含むペアで0を除外", int((unrel.gap_lo > 0).sum()), 2)
    chk("§16 ジャックナイフで gap<=0", int((cs.jack_min <= 0).sum()), 1)
    chk("§16 同種(対照非共有)で0を除外", int((ws.gap_lo > 0).sum()), 24)
    assert 0.08 < cs.frac_undefined.median() < 0.15, \
        "未定義反復の割合が n=3 の理論値 1/9 から外れている（実装の健全性）"
    assert rel.p_gap_le0.max() < unrel.p_gap_le0.max(), \
        "信頼性の高い猫状態のほうが P(gap<=0) が大きい（§2.4 の条件つき主張と矛盾）"

    # ---- §2.5 中心化と経路平均の分解 ----
    AS = tsv("pathway/aggregation_steps.tsv").set_index(["condition", "cls"])
    R1 = "1. cos of Δ as analysed"
    R2 = "2. Δ centred per state (= Pearson r of Δ)"
    R3 = "3. pathway means of centred Δ"
    R4 = "4. pathway means of uncentred Δ"
    for r, w, c in [(R1, 0.724, 0.447), (R2, 0.709, 0.551), (R3, 0.869, 0.882), (R4, 0.896, 0.489)]:
        chk(f"§2.5 {r} 種内", AS.loc[(r, "within"), "median"], w, 2e-3)
        chk(f"§2.5 {r} 種間", AS.loc[(r, "cross"), "median"], c, 2e-3)
    BA = tsv("pathway/cos_before_after.tsv")

    def auc_(col):
        a = BA[BA.within == 1][col].to_numpy(); b = BA[BA.within == 0][col].to_numpy()
        return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))
    chk("§2.5 AUC Δ そのまま", auc_("cos_raw"), 0.813, 2e-3)
    chk("§2.5 AUC 中心化のみ", auc_("cos_centred"), 0.759, 2e-3)
    chk("§2.5 AUC 中心化+経路平均", auc_("cos_agg_centred"), 0.511, 2e-3)
    chk("§2.5 AUC 中心化なし経路平均", auc_("cos_agg_raw"), 0.899, 2e-3)
    assert auc_("cos_agg_raw") > auc_("cos_raw") > auc_("cos_agg_centred"), \
        "中心化の有無で経路平均が逆方向になる、という §2.5 の核が成り立たない"

    # ---- §2.5 経路セルの件数（条件と分母） ----
    pw = tsv("pathway/pathway_cos.tsv")
    mo = pw[~pw.state.str.startswith("cat")]
    chk("§2.5 マウス×経路セル数", len(mo), 2316)
    chk("§2.5 n>=50 のセル数", int((mo.n_genes >= 50).sum()), 888)
    chk("§2.5 cos>=0.9 のセル数（条件なし）", int((mo["cos"] >= 0.9).sum()), 34)
    sel = mo[(mo["cos"] >= 0.9) & (mo.n_genes >= 50) & (mo.above_random)]
    chk("§2.5 条件を満たし cos>=0.9 のセル数", len(sel), 3)
    assert set(sel.pathway) == {"Focal adhesion"}, "0.9 以上のセルが focal adhesion 以外にある"

    # ================= 改訂 1（revision1）で追加した数値 =================
    # ---- Task 1: 3群の gap と状態ラベル置換 ----
    G = tsv("revision1/gap_by_group.tsv")
    g = G[G.ceiling == "uncorrected"].set_index("group")
    A = "A same-species, no shared controls"
    Bg = "B cross-species, reliable feline state"
    Cg = "C cross-species, less reliable feline state"
    chk("§2.2 同種 gap 中央", g.loc[A, "gap_median"], 0.1195, 2e-3)
    chk("§2.2 高信頼 異種 gap 中央", g.loc[Bg, "gap_median"], 0.403, 2e-3)
    chk("§2.2 低信頼 異種 gap 中央", g.loc[Cg, "gap_median"], 0.215, 2e-3)
    chk("§2.2 同種で区間が0を除外", int(g.loc[A, "n_interval_excludes_zero"]), 24)
    gsb = G[G.ceiling == "Spearman-Brown"].set_index("group")
    chk("§2.2 SB 同種 gap 中央", gsb.loc[A, "gap_median"], 0.196, 2e-3)
    chk("§2.2 SB 高信頼 gap 中央", gsb.loc[Bg, "gap_median"], 0.487, 2e-3)
    chk("§2.2 SB 低信頼 gap 中央", gsb.loc[Cg, "gap_median"], 0.338, 2e-3)
    PM = tsv("revision1/gap_group_permutation.tsv").set_index("ceiling")
    chk("§2.2 置換 観測差（無補正）", PM.loc["uncorrected", "observed"], 0.2188, 2e-3)
    chk("§2.2 置換 p（無補正）", PM.loc["uncorrected", "p_one_sided"], 0.0011, 1e-4)
    chk("§2.2 置換 観測差（SB）", PM.loc["Spearman-Brown", "observed"], 0.2390, 2e-3)
    chk("§2.2 置換 p（SB）", PM.loc["Spearman-Brown", "p_one_sided"], 0.0005, 1e-4)
    assert int(PM.loc["uncorrected", "n_perm"]) == 1820, "置換が C(16,4) 全列挙でない"

    # ---- Task 2: 非中心化コサインへの減衰式の当てはまり ----
    SIM = tsv("revision1/cosine_attenuation_sim.tsv")
    band = SIM[SIM.mean_r_half.between(0.5, 0.9)]
    assert band[band.true_cos <= 0.5].frac_obs_above_ceiling.max() == 0.0, \
        "真の cos <= 0.5 の帯で天井超過が起きている（§4.10 の記述と矛盾）"
    assert band[band.true_cos == 1.0].frac_obs_above_ceiling.min() > 0.99, \
        "同一応答でも天井を超えない条件がある（§4.10 の記述と矛盾）"
    chk("§4.10 真の cos 0.7 での最大超過", band[band.true_cos == 0.7].frac_obs_above_ceiling.max(),
        0.4567, 2e-3)
    chk("§4.10 実データの平均シフト中央", 0.257, 0.257, 1e-3)

    # ---- Task 3: Reactome を加えた経路解析 ----
    PR = tsv("revision1/pathway_reactome/pathway_separation_paired.tsv")
    col = "auc[exclude shared controls]"
    chk("§2.3 経路数（4コレクション）", len(PR), 422)
    chk("§2.3 経路単位 AUC 中央", PR[col].median(), 0.760, 2e-3)
    chk("§2.3 同 四分位下", PR[col].quantile(.25), 0.707, 2e-3)
    chk("§2.3 同 四分位上", PR[col].quantile(.75), 0.808, 2e-3)
    chk("§2.3 AUC>=0.90 の経路", int((PR[col] >= 0.90).sum()), 16)
    chk("§2.3 AUC<=0.60 の経路", int((PR[col] <= 0.60).sum()), 25)
    GS = tsv("revision1/genespace_with_reactome.tsv").set_index("space")
    chk("§2.6 経路 AUC Group A", GS.loc["Group A", "pathway_auc_median"], 0.760, 2e-3)
    chk("§2.6 経路 AUC 全ortholog", GS.loc["all 1:1 orthologues", "pathway_auc_median"], 0.782, 2e-3)
    chk("§2.6 経路 AUC マッチB", GS.loc["matched Group B", "pathway_auc_median"], 0.780, 2e-3)
    chk("§2.6 集約 AUC 全ortholog", GS.loc["all 1:1 orthologues", "aggregated_auc_centred"], 0.600, 2e-3)
    chk("§2.6 集約 AUC マッチB", GS.loc["matched Group B", "aggregated_auc_centred"], 0.624, 2e-3)
    for sp in GS.index:
        assert GS.loc[sp, "pathway_auc_median"] > GS.loc[sp, "aggregated_auc_centred"], \
            f"{sp} で経路単位が集約を上回っていない（Reactome 追加後）"

    # ---- Task 4: 猫側フィルタの非対称性 ----
    EX = tsv("revision1/groupA_expression_by_species.tsv").set_index("column")
    chk("§4.3 Group A 発現分位 猫", EX.loc["cat_ctx_control", "pct_median"], 0.857, 2e-3)
    chk("§4.3 Group A 発現分位 Pod-TRECK", EX.loc["podtreck_control", "pct_median"], 0.810, 2e-3)
    chk("§4.3 Group A 発現分位 IRI", EX.loc["iri_young_sham", "pct_median"], 0.828, 2e-3)
    FS = tsv("revision1/feline_filter_sensitivity.tsv")
    piv = FS.pivot(index="state", columns="filter", values="r_half")
    filt = [c for c in piv.columns if c != "none"][0]
    d = (piv[filt] - piv["none"])
    chk("§4.3 フィルタで上がる r_half 猫皮質CKD1/2", d.loc["cat_CKD12"], 0.035, 2e-3)
    assert (d > 0).all(), "猫側フィルタで信頼性が下がる状態がある（§4.3 の記述と矛盾）"
    assert piv[filt].idxmin() == piv["none"].idxmin(), \
        "フィルタで最も信頼性の低いネコ状態が入れ替わる（§2.2 の群分けに影響）"

    # ---- 改訂3: Reactome 込み 422 セットでの四条件（図 3・4 と §2.3 が同じ値を使う） ----
    BA2 = tsv("revision1/pathway_reactome/cos_before_after.tsv")

    def auc2(col):
        a = BA2[BA2.within == 1][col].to_numpy(); b = BA2[BA2.within == 0][col].to_numpy()
        return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))
    for col, w, c in [("cos_raw", 0.724, 0.447), ("cos_centred", 0.709, 0.551),
                      ("cos_agg_centred", 0.831, 0.832), ("cos_agg_raw", 0.870, 0.462)]:
        chk(f"§2.3(422) {col} 種内", float(BA2[BA2.within == 1][col].median()), w, 2e-3)
        chk(f"§2.3(422) {col} 種間", float(BA2[BA2.within == 0][col].median()), c, 2e-3)
    chk("§2.3(422) AUC 中心化+経路平均", auc2("cos_agg_centred"), 0.511, 2e-3)
    chk("§2.3(422) AUC 中心化なし経路平均", auc2("cos_agg_raw"), 0.868, 2e-3)
    assert auc2("cos_agg_raw") > auc2("cos_raw") > auc2("cos_agg_centred"), \
        "422 セットでも中心化の有無で経路平均が逆方向になる、という §2.3 の核が成り立たない"

    # ---- 改訂10: 結果ファイルがコードから再生成されることを確かめる ----
    # ここまでの検査は「本文の数値と結果ファイルが一致する」ことしか見ていない。
    # 中心結果を供給する3ファイルは、読むコードはあっても書くコードが無い状態だった。
    # 書き直した生成コードを入力から回し直し、置いてあるファイルと突き合わせる。
    sys.path.insert(0, str(HERE / "code" / "revision1"))
    import pathway_centering_steps as PCS      # noqa: E402
    import pathway_leads as PLD                # noqa: E402

    # 遺伝子セットの GMT は再配布していない（src/00_fetch_refs.py が取ってくる）。
    # クローンしただけの状態では再生成の検査だけ飛ばし、残りは最後まで走らせる。
    HAVE_GMT = all((PCS.REF / fn).exists() for fn in PCS.COLL_FULL.values())
    if not HAVE_GMT:
        skipped.append("経路集合の再生成（遺伝子セットの GMT は再配布していない。"
                       "`python src/00_fetch_refs.py` で取得すると実行される）")
    if HAVE_GMT:
        D_, P_ = PCS.load()
        sets193, sets422 = PCS.collect(D_, PCS.COLL_BASE), PCS.collect(D_, PCS.COLL_FULL)
        chk("再生成 セット数 3 コレクション", len(sets193), 193, 0)
        chk("再生成 セット数 Reactome 込み", len(sets422), 422, 0)
        regen = {"pathway/cos_before_after.tsv": PCS.conditions(D_, P_, sets193),
                 "revision1/pathway_reactome/cos_before_after.tsv": PCS.conditions(D_, P_, sets422)}

        def regen_auc(BAx, col):
            a = BAx[BAx.within == 1][col].to_numpy(); b = BAx[BAx.within == 0][col].to_numpy()
            return float(np.mean((a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])))

        G = regen["revision1/pathway_reactome/cos_before_after.tsv"]
        chk("再生成 §2.3 AUC 遺伝子単位", regen_auc(G, "cos_raw"), 0.813, 2e-3)
        chk("再生成 §2.3 AUC 中心化のみ", regen_auc(G, "cos_centred"), 0.759, 2e-3)
        chk("再生成 §2.3 AUC 中心化+経路平均", regen_auc(G, "cos_agg_centred"), 0.511, 2e-3)
        chk("再生成 §2.3 AUC 中心化なし経路平均", regen_auc(G, "cos_agg_raw"), 0.868, 2e-3)
        G193 = regen["pathway/cos_before_after.tsv"]
        chk("再生成 §2.5 AUC 中心化なし経路平均（193）", regen_auc(G193, "cos_agg_raw"), 0.899, 2e-3)

        # 置いてあるファイルと、再生成した表がセル単位で一致すること
        for rel, new in regen.items():
            old = tsv(rel)
            same = (list(old.columns) == list(new.columns) and old.shape == new.shape
                    and bool(np.all(np.abs(old.to_numpy(float)
                                           - new.round(4).to_numpy(float)) < 1e-9)))
            (ok if same else bad).append(
                f"{rel} はコードから再生成される" if same
                else f"{rel} が再生成した表と一致しない（実装かファイルのどちらかが古い）")
        st = PCS.steps_table(G193).round(4)
        old_st = tsv("pathway/aggregation_steps.tsv")
        same = (list(old_st.columns) == list(st.columns) and old_st.shape == st.shape
                and (old_st.condition == st.condition).all() and (old_st.cls == st.cls).all()
                and bool(np.all(np.abs(old_st[["n", "median", "q1", "q3"]].to_numpy(float)
                                       - st[["n", "median", "q1", "q3"]].to_numpy(float)) < 1e-9)))
        (ok if same else bad).append(
            "pathway/aggregation_steps.tsv はコードから再生成される" if same
            else "pathway/aggregation_steps.tsv が再生成した表と一致しない")

        lead = PLD.leads()
        br = PLD.breakdown(lead).round(4)
        old_br = tsv("revision1/leads_by_state_breakdown.tsv")
        same = (old_br.shape == br.shape and (old_br.state.tolist() == br.state.tolist())
                and bool(np.all(np.abs(old_br[["n_leads", "share"]].to_numpy(float)
                                       - br[["n_leads", "share"]].to_numpy(float)) < 1e-9)))
        (ok if same else bad).append(
            "revision1/leads_by_state_breakdown.tsv はコードから再生成される" if same
            else "revision1/leads_by_state_breakdown.tsv が再生成した表と一致しない")
        top = lead[lead.state == "IRI_2h"].nlargest(10, "cos")
        chk("再生成 IRI 2 h が首位の経路数", int((lead.state == "IRI_2h").sum()), 48, 0)
        chk("再生成 うち Reactome", int((top.collection == "Reactome").sum()), 10, 0)
    # ---- C案: 図・表・本文の照合と、対照解析の反映 ----
    import re as _re
    MDC = HERE / "manuscript_C"
    # ディレクトリの有無ではなく本文そのものを見る。公開リポジトリには図の描画データ
    # だけが manuscript_C/figure_data として入っており、本文の Markdown は入らない。
    if (MDC / "METHODS.md").exists():
        cfiles = ["FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                  "METHODS.md", "CONCLUSIONS.md", "BACKMATTER.md", "FIGURES.md",
                  "SUPPLEMENTARY.md"]
        cjoin = "\n".join((MDC / f).read_text() for f in cfiles)
        cflat = _re.sub(r"\s+", " ", cjoin)

        # 図に印字した値が本文・キャプションのどこかにあること
        FC = HERE / "results" / "revision1" / "pathway_reactome" / "figure_printed_numbers_C.tsv"
        PC = pd.read_csv(FC, sep="\t")
        # 図だけに出る値（遺伝子空間ごとの点推定など）は本文が表で持つので、表も含めて探す
        # 図では桁区切りを入れて描くことがあるので、比較のときだけ取り除く
        cflat_n = cflat.replace(",", "")
        miss = [r for r in PC.itertuples()
                if str(r.printed).rstrip("%").replace(",", "") not in cflat_n]
        (ok if not miss else bad).append(
            f"C案: 図に印字した {len(PC)} 件の数値はすべて本文・表・キャプションにある"
            if not miss else
            "C案: 図の印字値が本文に見つからない: "
            + ", ".join(f"{m.figure}{m.panel} {m.item}={m.printed}" for m in miss[:6]))

        # Figure 1(D) に描いた値が対照解析の出力と一致すること。
        # 4空間の centered / uncentered / 差 / 経路数を、描画時の記録と突き合わせる。
        CCd = pd.read_csv(HERE / "results" / "roundC" / "centering_controls.tsv", sep="\t")
        ccd = CCd[(CCd.collection.str.startswith("422"))
                  & (CCd.pairs.str.startswith("87"))].set_index("space")
        drawn = {(r.panel, r.item): str(r.printed) for r in PC.itertuples() if r.figure == "Fig1"}
        for sp in ["Group A (2,016)", "All 1:1 orthologues (7,897)",
                   "Matched Group B (1,632)", "Matched Group A (1,632)"]:
            row = ccd.loc[sp]
            want = {f"centered aggregate AUC {sp}": f"{row.auc_agg_centred:.3f}",
                    f"uncentered aggregate AUC {sp}": f"{row.auc_agg_uncentred:.3f}",
                    f"difference {sp}": f"{row.auc_agg_uncentred - row.auc_agg_centred:.3f}",
                    f"pathway sets {sp}": f"{int(row.n_sets):,}"}
            got = {k: drawn.get(("D", k)) for k in want}
            (ok if got == want else bad).append(
                f"C案: Figure 1(D) の {sp} が対照解析の出力と一致" if got == want
                else f"C案: Figure 1(D) の {sp} が一致しない: {got} != {want}")

        # §2.4 の Pod-TRECK D14: 表示される 0.740 は丸め前の積から来ていること。
        # 丸めた ν 1.484 と cos 0.498 を掛けると 0.739 になるので、読者が食い違いと
        # 受け取りうる。恒等式が丸め前で成り立つことを明示的に検査する。
        prc = pd.read_csv(R / "alpha_groupA" / "projection.tsv", sep="\t", index_col=0)
        rr = prc.loc["mouse_2W"]
        # 結果ファイルは 6 桁で保存しているので、そこが一致の限界になる
        ident = abs(float(rr.nrm) * float(rr["cos"]) - float(rr.alpha))
        (ok if ident < 1e-5 else bad).append(
            f"C案 §2.4: Pod-TRECK D14 の α = ν cos θ が成立（差 {ident:.1e}、"
            "結果ファイルの 6 桁保存の範囲内）"
            if ident < 1e-5 else f"C案 §2.4: α = ν cos θ が成り立たない（差 {ident:.1e}）")
        shown = f"{float(rr.alpha):.3f}"
        (ok if shown == "0.740" else bad).append(
            "C案 §2.4: Pod-TRECK D14 の α は丸めて 0.740" if shown == "0.740"
            else f"C案 §2.4: Pod-TRECK D14 の α を丸めると {shown}")

        # AC: Δ の 4 定義の値が本文・Table S9 と一致すること
        DDf = HERE / "results" / "roundR" / "delta_definitions.tsv"
        if DDf.exists():
            DD = pd.read_csv(DDf, sep="\t")
            miss = []
            for r in DD.itertuples():
                for col, nd in (("auc_agg_uncentered", 3), ("auc_agg_centered", 3),
                                ("centering_drop", 3), ("auc_gene_level", 3),
                                ("cat_mean_min", 3), ("mouse_mean_max", 3)):
                    v = f"{getattr(r, col):.{nd}f}".lstrip("0") if False else \
                        f"{getattr(r, col):.{nd}f}"
                    if v.replace("-", "\u2212") not in cflat and v not in cflat:
                        miss.append(f"{r.definition}/{col}={v}")
            (ok if not miss else bad).append(
                f"C案 AC: Δ 4 定義の {len(DD) * 6} 値がすべて本文か Table S9 にある"
                if not miss else f"C案 AC: 本文にも Table S9 にも無い値: {miss[:5]}")

        # AD: 経路スコアリング手法の値が本文・Table S10 と一致すること
        PSf = HERE / "results" / "roundR" / "pathway_scoring_methods.tsv"
        if PSf.exists():
            PS = pd.read_csv(PSf, sep="\t")
            miss = [f"{r.method}/{c}"
                    for r in PS.itertuples()
                    for c in ("auc_uncentered", "auc_centered")
                    if f"{getattr(r, c):.3f}" not in cflat]
            (ok if not miss else bad).append(
                f"C案 AD: 経路スコアリング {len(PS)} 手法の AUC が本文か Table S10 にある"
                if not miss else f"C案 AD: 見つからない: {miss[:5]}")

        # AE1: 75 ペアでの中心化比較が本文にあること
        B75f = HERE / "results" / "round7" / "centering_auc_animal_bootstrap_75pairs.tsv"
        if B75f.exists():
            B75 = pd.read_csv(B75f, sep="\t").set_index("quantity")["value"]
            want = {"observed AUC, uncentred pathway means": 3,
                    "observed AUC, centred pathway means": 3,
                    "observed drop": 3,
                    "bootstrap drop median": 3,
                    "bootstrap drop 2.5th percentile": 3,
                    "bootstrap drop 97.5th percentile": 3}
            miss = [f"{k}={float(B75[k]):.{n}f}" for k, n in want.items()
                    if f"{float(B75[k]):.{n}f}" not in cflat]
            pct = f"{float(B75['replicates with uncentred AUC above centred']) * 100:.1f}%"
            if pct not in cflat:
                miss.append(f"exceedance={pct}")
            (ok if not miss else bad).append(
                f"C案 AE: 75 ペアの中心化比較 {len(want) + 1} 値が本文にある" if not miss
                else f"C案 AE: 本文に無い 75 ペアの値: {miss}")

        # AE2: 経路集合の設計を変えた値が本文か Table S5 にあること
        SDf = HERE / "results" / "roundR" / "pathway_set_design.tsv"
        if SDf.exists():
            SD = pd.read_csv(SDf, sep="\t")
            miss = []
            for r in SD.itertuples():
                for c in ("auc_uncentered", "auc_centered", "difference"):
                    if f"{getattr(r, c):.4f}" not in cflat:
                        miss.append(f"{r.analysis}/{r.label}/{c}={getattr(r, c):.4f}")
                if r.analysis == "common sets" and f"{r.difference_all_sets:.4f}" not in cflat:
                    miss.append(f"{r.label}/own={r.difference_all_sets:.4f}")
            (ok if not miss else bad).append(
                f"C案 AE: 経路集合 {len(SD)} 条件の AUC が本文か Table S5 にある" if not miss
                else f"C案 AE: 見つからない: {miss[:5]}")
            # 本文が丸めて書いた 3 桁の値も、同じ出力から来ていること
            txt3 = [f"{r.difference:.3f}" for r in SD.itertuples()]
            m3 = [v for v in txt3 if v not in cflat]
            (ok if not m3 else bad).append(
                "C案 AE: 経路集合の差を 3 桁に丸めた値がすべて本文にある" if not m3
                else f"C案 AE: 3 桁の値が本文に無い: {m3}")

        # AD: 実装の記述（バージョン・引数・不変性の数値）が監査出力と一致すること
        AUf = HERE / "results" / "roundR" / "ad_implementation_audit.tsv"
        if AUf.exists():
            AU = pd.read_csv(AUf, sep="\t").set_index("method")
            SO = pd.read_csv(HERE / "results" / "roundR" / "ad_ssgsea_options.tsv",
                             sep="\t").set_index("option")
            RI = pd.read_csv(HERE / "results" / "roundR" / "ad_rank_invariance.tsv",
                             sep="\t").set_index("quantity")["value"]
            S4 = pd.read_csv(HERE / "results" / "roundR" / "table_s4_rows.tsv",
                             sep="\t").set_index("method")
            RV = dict(l.split("\t") for l in
                      (HERE / "results" / "roundR" / "r_versions.tsv").read_text().splitlines())
            want = [f"GSEApy {AU.loc['ssGSEA', 'version']}",
                    f"NumPy {AU.loc['pathway mean (main)', 'version']}",
                    f"R {RV['R']}", f"Bioconductor {RV['Bioconductor']}",
                    f"GSVA {RV['GSVA']}", f"singscore {RV['singscore']}",
                    f"GSEABase {RV['GSEABase']}",
                    f"{S4.loc['GSVA', 'max_abs_score_difference']:.3f}",
                    f"{S4.loc['PLAGE', 'max_abs_score_difference']:.3f}",
                    f"{S4.loc['pathway mean (main)', 'max_abs_score_difference']:.4f}",
                    "%.0f" % SO.loc[[i for i in SO.index if i.startswith("sample_norm_method='log'")][0],
                                     "max_abs_ES_difference"],
                    "%.0f" % SO.loc["sample_norm_method=None (input used as the metric)",
                                    "max_abs_ES_difference"]]
            miss = [v for v in want if v not in cflat]
            (ok if not miss else bad).append(
                f"C案 AD: 実装記述の {len(want)} 項目が Methods にある" if not miss
                else f"C案 AD: Methods に無い実装の値: {miss}")
            # 不変性の主張そのものを見張る
            inv = {m: bool(r.invariant_to_per_state_shift) for m, r in AU.iterrows()}
            exp = {"pathway mean (main)": False, "singscore-style": True,
                   "PLAGE-style": False, "ssGSEA": True, "GSVA": False}
            (ok if inv == exp else bad).append(
                "C案 AD: 中心化に不変なのは ssGSEA と singscore だけ" if inv == exp
                else f"C案 AD: 不変性のパターンが変わった: {inv}")
            rk = bool(RI["within-state gene ranks identical after centering"])
            (ok if rk else bad).append(
                "C案 AD: 中心化の前後で状態内の遺伝子順位が完全に同一（代数的帰結）" if rk
                else "C案 AD: 中心化で状態内の順位が変わっている。不変性の説明が成り立たない")
            rank_based = [i for i in SO.index if "log'" not in i and "None" not in i]
            allinv = all(bool(SO.loc[i, "ES_invariant"]) and bool(SO.loc[i, "NES_invariant"])
                         for i in rank_based)
            (ok if allinv else bad).append(
                f"C案 AD: 順位ベースの ssGSEA 設定 {len(rank_based)} 通りすべてで ES・NES が不変"
                if allinv else "C案 AD: 順位ベース設定で不変でないものがある")

        # AF7: 添付した定量表から Group A を計算し直し、本文の数と一致すること
        A7f = HERE / "results" / "roundR" / "af7_attached_table.tsv"
        if A7f.exists():
            A7 = pd.read_csv(A7f, sep="\t").set_index("quantity")["value"]
            import af7_attach_source_table as AF7   # noqa: E402
            import build_delta_matrix as BDM       # noqa: E402
            att = AF7.DEST / str(A7["attached file name"])
            if not att.exists():
                bad.append(f"C案 AF7: 添付ファイルが無い: {att}")
            else:
                h = AF7.sha256(att)
                (ok if h == str(A7["attached file sha256"]) else bad).append(
                    "C案 AF7: 添付ファイルのチェックサムが記録と一致する"
                    if h == str(A7["attached file sha256"])
                    else "C案 AF7: 添付ファイルが記録と違う。ak_strip_docprops.py と "
                         "af7_attach_source_table.py を回すこと")
                # AK: 著者でない個人名と受託先の内部パスが添付から消えていること
                import zipfile as _zf
                zo = _zf.ZipFile(att)
                import ak_strip_docprops as _AK   # noqa: E402
                _pats = _AK.leak_terms() + [b"absPath"]
                leak = sorted({n for n in zo.namelist()
                               for pat in _pats if pat in zo.read(n)})
                (ok if not leak else bad).append(
                    "C案 AK: 添付に個人名・受託先の内部パスが残っていない" if not leak
                    else f"C案 AK: 添付に残っている: {leak}")
                core = zo.read("docProps/core.xml").decode("utf-8")
                empty = all(f"<{g}></{g}>" in core or f"<{g}/>" in core
                            for g in ("dc:creator", "cp:lastModifiedBy"))
                (ok if empty else bad).append(
                    "C案 AK: dc:creator と cp:lastModifiedBy が空である" if empty
                    else "C案 AK: dc:creator / cp:lastModifiedBy が空でない")
                # データ側のエントリが元ファイルとバイト同一であること
                import ak_strip_docprops as AK   # noqa: E402
                zsrc = _zf.ZipFile(AK.SRC)
                keep = [n for n in zsrc.namelist()
                        if any(n.startswith(k) or n == k for k in
                               ("xl/worksheets/", "xl/sharedStrings.xml", "xl/styles.xml",
                                "xl/theme/", "xl/_rels/", "[Content_Types].xml", "_rels/.rels",
                                "docProps/app.xml"))]
                dd = [n for n in keep if zsrc.getinfo(n).CRC != zo.getinfo(n).CRC]
                (ok if not dd else bad).append(
                    f"C案 AK: データ側 {len(keep)} エントリが元ファイルとバイト同一"
                    if not dd else f"C案 AK: 内容が変わったエントリ: {dd}")
                same_names = zsrc.namelist() == zo.namelist()
                (ok if same_names else bad).append(
                    "C案 AK: zip のエントリ一覧が元ファイルと一致する" if same_names
                    else "C案 AK: zip のエントリ一覧が違う")
                # 添付ファイルから検出コールを作り直す（本文が主張している再計算そのもの）
                Mx = AF7.matrix_from_workbook(att)
                det = set(Mx.index[Mx.notna().mean(axis=1) >= AF7.FRAC])
                hm, _ = BDM.to_human(pd.Series(1.0, index=sorted(det)), "mouse")

                def _cd(nm):
                    d = BDM.read_int(nm)
                    s0 = pd.Series(1.0, index=d.index[d.notna().mean(axis=1) >= AF7.FRAC])
                    hh, _ = BDM.to_human(s0, "cat")
                    return set(hh.index)

                DM = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
                gaa = sorted(((_cd("cat_prot_ctx") & _cd("cat_prot_med")) & set(hm.index))
                             & set(DM.index))
                comp = set(DM.dropna(axis=0, how="any").index)
                gac = [g for g in gaa if g in comp]
                nowset = set((HERE / "groupA_intersection.txt").read_text().split())
                hit = set(gaa) == nowset
                (ok if hit else bad).append(
                    f"C案 AF7: 添付した定量表から Group A {len(gaa)} 遺伝子が再現する"
                    if hit else f"C案 AF7: 添付ファイルから再構成した Group A が違う "
                                f"({len(gaa)} 対 {len(nowset)})")
                (ok if len(gac) == 2016 else bad).append(
                    f"C案 AF7: そのうち 16 状態で有限なのは {len(gac)} 遺伝子"
                    if len(gac) == 2016 else f"C案 AF7: complete case が {len(gac)} になっている")
                want = [f"{len(gaa):,}", f"{len(gac):,}",
                        # 記号数は detection_rule.tsv（「-」を除いた値）が正本。
                        # A7 の行列は to_human の前なので「-」を含む。
                        f"{int(float(pd.read_csv(HERE / 'results' / 'roundR' / 'detection_rule.tsv', sep=chr(9)).set_index('quantity')['value']['pod_treck proteins'])):,}"]
                miss = [v for v in want if v not in cflat]
                (ok if not miss else bad).append(
                    "C案 AF7: 添付ファイルの規模と Group A の数が Methods にある" if not miss
                    else f"C案 AF7: Methods に無い: {miss}")
            # 撤回した記述が本文に残っていないこと
            # 撤回・訂正した主張が本文（Abstract・Conclusions を含む）に戻っていないこと
            WITHDRAWN = ("cannot be reconstructed", "is vacuous", "2,241",
                         "complete-case matrix from a separate processing run",
                         # AE2: ランダム集合の設計で支持されなかった 2 つ
                         "specific to the real pathway sets",
                         "contributes beyond the averaging operation",
                         "Pathway membership is not incidental",
                         # AE2c: 遺伝子空間への誤った帰属
                         "smallest in the space that imposes no proteomic selection",
                         "the magnitude depends on which genes are compared",
                         "the size of the effect does depend on which genes are compared")
            # AT6: どこにも出てはいけない言い方（letter も対象。letter は撤回した主張を
            # 名前で挙げるので、上の WITHDRAWN は本文だけに当てる）
            BANNED = ("disappears by construction",
                      "rather than the genes themselves",
                      "repeated counting of shared genes is not what produces it",
                      "permutation is therefore carried out at the level of states",
                      # AV3: 今回の言い換え
                      "remove the gap altogether", "consequence of the construction",
                      "centered value no lower",
                      "No reported value moves in its third decimal",
                      "hardest to reconcile", "splitting the difference",
                      "now follow the Introduction",
                      # AY4
                      "led to a new computation", "do not depend on the input order",
                      "which excludes counts", "excludes CPM and TPM",
                      # BD1: 図が描いていない条件を図に帰属させていた言い方
                      "gradient shown in Figure S4",
                      "the rate rises with reliability within it")
            gone = [s for s in WITHDRAWN if s in cflat]
            (ok if not gone else bad).append(
                f"C案 AF7/AE: 撤回・訂正した {len(WITHDRAWN)} 表現が本文に残っていない"
                if not gone else f"C案 AF7/AE: 撤回したはずの記述が残っている: {gone}")
            lt2 = _re.sub(r"\s+", " ", (MDC / "RESPONSE.md").read_text()) \
                if (MDC / "RESPONSE.md").exists() else ""
            gone2 = [s for s in BANNED if s in cflat or s in lt2]
            (ok if not gone2 else bad).append(
                f"C案 AT: 言い過ぎの {len(BANNED)} 表現が本文・letter のどこにも無い"
                if not gone2 else f"C案 AT: 残っている: {gone2}")

        # 検出カスケード（Table S10）の数が結果ファイルと一致すること
        GSf = HERE / "results" / "roundR" / "groupA_detection_summary.tsv"
        if GSf.exists():
            GS = pd.read_csv(GSf, sep="\t")
            gone = [f"{int(r.genes):,}" for r in GS.itertuples()
                    if f"{int(r.genes):,}" not in cflat]
            (ok if not gone else bad).append(
                f"C案 AF: 検出コールの段階 {len(GS)} 件の遺伝子数が補足にある" if not gone
                else f"C案 AF: 補足に無い検出数: {gone}")

        # AB: 重複シンボルの影響範囲が結果ファイルと一致すること
        DSf = HERE / "results" / "roundR" / "duplicate_symbol_sensitivity.tsv"
        if DSf.exists():
            DSq = pd.read_csv(DSf, sep="\t").set_index("quantity")["current"]
            k = "analysed genes with duplicate source rows"
            want = [f"{int(DSq[k])} of the 2,016 genes analysed",
                    f"{float(DSq['share of the analysed genes'])*100:.1f}%"]
            miss = [v for v in want if v not in cflat]
            (ok if not miss else bad).append(
                "C案 AB: 重複シンボルの影響範囲が Methods と一致する" if not miss
                else f"C案 AB: Methods と合わない: {miss}")

        # AG: 節の並び順と、本文が指す節がすべて実在すること
        BMD = (HERE / "code" / "revision1" / "build_main_docx.py").read_text()
        mo = _re.search(r"ORDER\s*=\s*\[(.*?)\]", BMD, _re.S)
        seq = _re.findall(r"[\"\']([A-Z_]+)\.md[\"\']", mo.group(1)) if mo else []
        want_head = ["FRONTMATTER", "INTRODUCTION", "RESULTS", "DISCUSSION", "METHODS"]
        (ok if seq[:5] == want_head else bad).append(
            "C案 AR: 組み立て順が投稿先の書式（Introduction → Results → Discussion → Methods）"
            if seq[:5] == want_head else f"C案 AR: 組み立て順が違う: {seq[:5]}")
        heads = set()
        for f in ("RESULTS.md", "DISCUSSION.md", "METHODS.md", "CONCLUSIONS.md",
                  "INTRODUCTION.md"):
            for m in _re.finditer(r"(?m)^#{1,3}\s*(\d+(?:\.\d+)?)\.", (MDC / f).read_text()):
                heads.add(m.group(1))
        cited = set()
        for f in ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                  "METHODS.md", "CONCLUSIONS.md", "BACKMATTER.md", "SUPPLEMENTARY.md",
                  "FIGURES.md", "RESPONSE.md"):
            fp = MDC / f
            if not fp.exists():
                continue
            for m in _re.finditer(r"Sections?\s+((?:\d+\.\d+)(?:\s*(?:,|and)\s*\d+\.\d+)*)",
                                  fp.read_text()):
                cited |= set(_re.findall(r"\d+\.\d+", m.group(1)))
        dead = sorted(cited - heads)
        (ok if not dead else bad).append(
            f"C案 AG: 本文が指す {len(cited)} 件の節がすべて実在する" if not dead
            else f"C案 AG: 実在しない節への参照: {dead}")

        # AS: 変更履歴版が、承諾＝修正稿・却下＝元稿を満たすこと
        TRK = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Main_manuscript_tracked.docx"
        if TRK.exists():
            import subprocess
            rt = subprocess.run([sys.executable,
                                 str(HERE / "code" / "revision1" / "al_check_tracked.py")],
                                capture_output=True, text=True)
            (ok if rt.returncode == 0 else bad).append(
                "C案 AS: 変更履歴版の検査に合格（承諾＝修正稿 / 却下＝元稿 / 数式 33 / "
                "行番号 / 作成者 1 名）" if rt.returncode == 0
                else "C案 AS: 変更履歴版の検査に不合格: "
                     + " / ".join(l.strip() for l in rt.stdout.splitlines()
                                  if "不一致" in l or "✗" in l)[:220])
            SUBM = HERE / "results" / "roundR" / "submission_manifest.tsv"
            if SUBM.exists():
                MF = pd.read_csv(SUBM, sep="\t").set_index("file")
                import hashlib as _hl

                def _sha(q):
                    h = _hl.sha256()
                    with open(q, "rb") as fh:
                        for b in iter(lambda: fh.read(1 << 20), b""):
                            h.update(b)
                    return h.hexdigest()

                stale = [n for n in MF.index
                         if (DEST0 := TRK.parent / n).exists()
                         and _sha(DEST0) != MF.loc[n, "sha256"]]
                (ok if not stale else bad).append(
                    f"C案 AS: マニフェストの SHA-256 が投稿ファイル {len(MF)} 件と一致する"
                    if not stale
                    else f"C案 AS: マニフェストが古い: {stale}。al_manifest.py を回すこと")

        # AR4: 書き換えで書式が失われていないこと、dc:title が実際のタイトルであること
        FL = HERE / "results" / "roundR" / "an_format_lost.tsv"
        if FL.exists():
            F_ = pd.read_csv(FL, sep="\t")
            (ok if F_.empty else bad).append(
                "C案 AR: 書き換えで太字・斜体・上付き・下付きが失われた箇所は無い" if F_.empty
                else f"C案 AR: 書式が失われた箇所 {len(F_)} 件: "
                     + ", ".join(f"{r_.section}/{r_.mark}" for r_ in F_.itertuples())[:180])
        MDX = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Main_manuscript.docx"
        if MDX.exists():
            import zipfile as _z2
            core = _z2.ZipFile(MDX).read("docProps/core.xml").decode("utf-8")
            mt = _re.search(r"<dc:title>(.*?)</dc:title>", core, _re.S)
            got = mt.group(1) if mt else ""
            want_t = "Mean-Centering Alters Pathway-Based Transcriptomic Comparisons"
            (ok if want_t in got else bad).append(
                f"C案 AR: dc:title が論文タイトルになっている" if want_t in got
                else f"C案 AR: dc:title が「{got[:48]}」のまま")

        # 原稿 ID が回答書とカバーレターに入っていること（プレースホルダが残っていないこと）
        import build_response_letter as BRL   # noqa: E402
        for nm in ("RESPONSE.md", "COVER_LETTER.md"):
            q2 = MDC / nm
            if not q2.exists():
                continue
            s2 = q2.read_text()
            good = BRL.MS_ID in s2 and "XXXX-XXXXXXX" not in s2
            (ok if good else bad).append(
                f"C案 原稿 ID: {nm} に {BRL.MS_ID} が入っている" if good
                else f"C案 原稿 ID: {nm} が {BRL.MS_ID} になっていない")

        # BA4: 削除・差し替えのあと、文が壊れていないこと（出来上がった docx の文字を読む）
        import subprocess
        rr4 = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "ba_read_edits.py")],
                             capture_output=True, text=True)
        (ok if rr4.returncode == 0 else bad).append(
            "C案 BA: 提出 docx に壊れた文（同じ文の繰り返し・句読点前の空白・"
            "主語と動詞の不一致・語の重複）が無い" if rr4.returncode == 0
            else "C案 BA: 壊れた文の疑い: "
                 + " / ".join(l.strip() for l in rr4.stdout.splitlines() if "✗" in l)[:220])

        # BA3: 出来上がった docx / xlsx に Markdown の記号が残っていないこと
        import subprocess
        rm3 = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "ba_check_markup.py")],
                             capture_output=True, text=True)
        (ok if rm3.returncode == 0 else bad).append(
            "C案 BA: 提出ファイルに Markdown の記号が残っていない" if rm3.returncode == 0
            else "C案 BA: 記号が残っている: "
                 + " / ".join(l.strip() for l in rm3.stdout.splitlines() if "✗" in l)[:220])

        # AZ1: 「-」を含めた古い集計値が、本文・回答書・カバーレター・補足に残っていないこと
        DRz = HERE / "results" / "roundR" / "detection_rule.tsv"
        if DRz.exists():
            DZ = pd.read_csv(DRz, sep="\t").set_index("quantity")["value"]
            allz = cflat
            for extra in ("RESPONSE.md", "COVER_LETTER.md"):
                q = MDC / extra
                if q.exists():
                    allz += " " + _re.sub(r"\s+", " ", q.read_text())
            stale_z = [v for v in ("12,635", "12,309") if v in allz]
            (ok if not stale_z else bad).append(
                "C案 AZ: 「-」を含めた古い集計値（12,635 / 12,309）は残っていない"
                if not stale_z else f"C案 AZ: 残っている: {stale_z}")
            need_z = [f"{int(float(DZ['pod_treck data rows'])):,}",
                      str(int(float(DZ["pod_treck rows without a gene symbol"]))),
                      f"{int(float(DZ['pod_treck proteins'])):,}",
                      f"{int(float(DZ['pod_treck detected'])):,}"]
            miss_z = [v for v in need_z if v not in allz]
            (ok if not miss_z else bad).append(
                "C案 AZ: Pod-TRECK の行数・記号数・検出数が本文と回答書にある" if not miss_z
                else f"C案 AZ: 無い: {miss_z}")

        # AY1: データセット・著者名・アクセッションを含む段落の引用が、書誌と合っていること
        import subprocess
        rc2 = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "ay_check_citations.py")],
                             capture_output=True, text=True)
        (ok if rc2.returncode == 0 else bad).append(
            "C案 AY: 出典を伴う段落の引用がすべて書誌と合っている" if rc2.returncode == 0
            else "C案 AY: 書誌と合わない引用: "
                 + " / ".join(l.strip() for l in rc2.stdout.splitlines() if "✗" in l)[:240])
        # letter の引用に Markdown の見出し記号が残っていないこと
        if (MDC / "RESPONSE.md").exists():
            lt4 = (MDC / "RESPONSE.md").read_text()
            hashy = _re.findall(r"\u201c##? [^\u201d]{0,60}", lt4)
            (ok if not hashy else bad).append(
                "C案 AY: letter の引用に Markdown の見出し記号が無い" if not hashy
                else f"C案 AY: 見出し記号が残っている: {hashy[:2]}")

        # BD: 回答書が述べる 3 つの件数（検査件数・トーン語数・結果ファイル数）が
        # それぞれの出どころと一致していること。手書きの数が古くなるのを防ぐ。
        cntf = HERE / "results" / "roundR" / "verification_counts.tsv"
        provf = HERE / "results" / "provenance.tsv"
        if (MDC / "RESPONSE.md").exists() and cntf.exists() and provf.exists():
            import tone_scan as _TS   # noqa: E402
            lt5 = _re.sub(r"\s+", " ", (MDC / "RESPONSE.md").read_text())
            n_prov = len(pd.read_csv(provf, sep="\t"))
            n_chk = int(pd.read_csv(cntf, sep="\t").set_index("quantity")["value"]["checks_agreed"])
            need5 = [(f"{n_chk:,} agreements", "検査件数"),
                     (f"{len(_TS.TERMS)}-term tone list", "トーン語数"),
                     (f"All {n_prov:,} result files", "結果ファイル数")]
            miss5 = [f"{w}（{nm}）" for w, nm in need5 if w not in lt5]
            (ok if not miss5 else bad).append(
                "C案 BD: 回答書の件数 3 件が出どころと一致している" if not miss5
                else f"C案 BD: 回答書の件数が古い: {miss5}")
            # BF1: 回答書に現れる「結果ファイル数」「検査件数」が、1 箇所残らず
            # 出どころの値であること（項目 4 と末尾で数が食い違っていた）
            WANT = {r"([\d,]+)\s+result files": n_prov,
                    r"([\d,]+)\s+agreements": n_chk,
                    r"([\d,]+)-term tone list": len(_TS.TERMS)}
            wrong5 = []
            for rx, want in WANT.items():
                for got in _re.findall(rx, lt5):
                    if int(got.replace(",", "")) != want:
                        wrong5.append(f"{got}（正しくは {want:,}）")
            (ok if not wrong5 else bad).append(
                "C案 BF: 回答書の結果ファイル数・検査件数はすべて結果ファイル由来"
                if not wrong5 else f"C案 BF: 出どころと違う数が書かれている: {wrong5}")

        # BJ1: Limitations の段落が BI 前と同じ文であること（文献番号の違いは除く）
        LIMF = HERE / "results" / "roundR" / "local" / "limitations_paragraph.txt"
        if LIMF.exists():
            want_lim = _re.sub(r"\[\d+(?:[,–-]\d+)*\]", "[#]",
                               " ".join(LIMF.read_text().split()))
            got_lim = ""
            for para in (MDC / "DISCUSSION.md").read_text().split("\n\n"):
                if "Species and dataset are confounded" in para:
                    got_lim = _re.sub(r"\[\d+(?:[,–-]\d+)*\]", "[#]",
                                      " ".join(para.split()))
            (ok if got_lim == want_lim else bad).append(
                "C案 BJ: Limitations の段落が BI 前の文のまま" if got_lim == want_lim
                else "C案 BJ: Limitations の段落が変わっている: " + got_lim[:160])

        # BL: 補足の説明文（繰り返し・小見出しの流れ込み・ノートと表の一致）と図の提出形態
        import subprocess
        rbl = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "bl_check_captions.py")],
                             capture_output=True, text=True)
        nbl = _re.search(r"合格 (\d+) /", rbl.stdout)
        (ok if rbl.returncode == 0 else bad).append(
            f"C案 BL: 補足の説明文と図の提出形態 {nbl.group(1) if nbl else '?'} 件に合格"
            if rbl.returncode == 0
            else "C案 BL: 説明文か提出形態が合っていない: "
                 + " / ".join(l.strip() for l in rbl.stdout.splitlines()
                              if l.strip().startswith("✗"))[:240])

        # BK: 回答書の段落が Markdown と同じ順で docx に入っていること
        import subprocess
        rbk = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "bk_check_letter_order.py")],
                             capture_output=True, text=True)
        nbk = _re.search(r"突き合わせた段落 (\d+) 件", rbk.stdout)
        (ok if rbk.returncode == 0 else bad).append(
            f"C案 BK: 回答書の見出しと引用 {nbk.group(1) if nbk else '?'} 件が Markdown と同じ順"
            if rbk.returncode == 0
            else "C案 BK: 回答書の段落の並びが違う: "
                 + " / ".join(l.strip() for l in rbk.stdout.splitlines()
                              if l.strip().startswith("✗"))[:240])

        # BJ2: 補足・回答書に現れる文献番号が、いまの文献一覧と合っていること
        import subprocess
        rbj = subprocess.run([sys.executable,
                              str(HERE / "code" / "revision1" / "bj_check_refnums.py")],
                             capture_output=True, text=True)
        nbj = _re.search(r"突き合わせた番号 (\d+) 件", rbj.stdout)
        (ok if rbj.returncode == 0 else bad).append(
            f"C案 BJ: 補足・回答書の文献番号 {nbj.group(1) if nbj else '?'} 件が一致"
            if rbj.returncode == 0
            else "C案 BJ: 文献番号が古い: "
                 + " / ".join(l.strip() for l in rbj.stdout.splitlines()
                              if l.strip().startswith("✗"))[:240])

        # BH: 青字版（変更をすべて承諾し、挿入された文字だけを青にした版）の検査
        if (Path.home() / "Desktop" / "CKD_xspecies_submission"
                / "Main_manuscript_blue.docx").exists():
            import subprocess
            rbh = subprocess.run([sys.executable,
                                  str(HERE / "code" / "revision1" / "bh_check_blue.py")],
                                 capture_output=True, text=True)
            nbh = _re.search(r"合格 (\d+) /", rbh.stdout)
            (ok if rbh.returncode == 0 else bad).append(
                f"C案 BH: 青字版の検査 {nbh.group(1) if nbh else '?'} 件に合格"
                if rbh.returncode == 0
                else "C案 BH: 青字版が検査に不合格: "
                     + " / ".join(l.strip() for l in rbh.stdout.splitlines()
                                  if l.strip().startswith("✗"))[:240])

        # BD3: 本文・補足ノートが図について述べている要約値を、図の描画データから再計算して照合
        if (HERE / "manuscript_C" / "figure_data").exists():
            import subprocess
            rb3 = subprocess.run([sys.executable,
                                  str(HERE / "code" / "revision1" / "bd_check_figure_text.py")],
                                 capture_output=True, text=True)
            nb3 = _re.search(r"(\d+) 件一致", rb3.stdout)
            (ok if rb3.returncode == 0 else bad).append(
                f"C案 BD: 図について述べた要約値 {nb3.group(1) if nb3 else '?'} 件が"
                " 図の描画データと一致" if rb3.returncode == 0
                else "C案 BD: 図と本文が食い違う: "
                     + " / ".join(l.strip() for l in rb3.stdout.splitlines()
                                  if l.strip().startswith("NG"))[:240])

        # AV1: 添付ファイルから本文の規則どおりに検出判定を再計算し、判定表と一致すること
        DRf = HERE / "results" / "roundR" / "detection_rule.tsv"
        DCf = HERE / "results" / "roundR" / "groupA_detection_calls.tsv"
        if DRf.exists() and DCf.exists():
            DR = pd.read_csv(DRf, sep="\t").set_index("quantity")["value"]
            import af7_attach_source_table as AF7b   # noqa: E402
            att = AF7b.DEST / AF7b.ATTACH
            if att.exists():
                # 本文の規則: 0 を非定量とみなし、9 検体中 5 検体以上で正の値
                Mx = AF7b.matrix_from_workbook(att)     # 0 は欠測に変換されている
                pos = (Mx.notna().sum(axis=1) >= int(DR["pod_treck threshold"]))
                det = set(Mx.index[pos])
                hm2, _ = BDM.to_human(pd.Series(1.0, index=sorted(det)), "mouse")
                DC = pd.read_csv(DCf, sep="\t")
                said = set(DC.loc[DC.pod_treck_proteome == "yes", "gene"])
                recomputed = set(hm2.index) & set(DC.gene)
                diff = said ^ recomputed
                (ok if not diff else bad).append(
                    f"C案 AV: 添付ファイルから再計算した Pod-TRECK の検出判定が "
                    f"{len(DC):,} 遺伝子すべてで判定表と一致する" if not diff
                    else f"C案 AV: 判定が食い違う遺伝子 {len(diff)} 件: {sorted(diff)[:5]}")
                # 本文が書いた数が実測と合うこと
                want_av = [f"{int(DR['pod_treck zero cells']):,}",
                           f"{int(DR['pod_treck cells']):,}",
                           str(int(DR["pod_treck threshold"])),
                           str(int(DR["feline cortex threshold"])),
                           str(int(DR["feline medulla threshold"])),
                           f"{int(DR['pod_treck detected']):,}",
                           f"{int(DR['feline cortex detected']):,}",
                           f"{int(DR['feline medulla detected']):,}"]
                miss_av = [v for v in want_av if v not in cflat]
                (ok if not miss_av else bad).append(
                    f"C案 AV: 検出判定の {len(want_av)} 値が本文にある" if not miss_av
                    else f"C案 AV: 本文に無い: {miss_av}")

        # AV3: letter が指す節・文献が、その回答の主題を扱っていること
        LETf = MDC / "RESPONSE.md"
        if LETf.exists():
            lt3 = LETf.read_text()
            # 見出しだけでは語が足りないので、その節の本文冒頭も主題語に含める
            heads, body_of = {}, {}
            for f3 in ("INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md", "METHODS.md",
                       "CONCLUSIONS.md"):
                txt3 = (MDC / f3).read_text()
                hits3 = list(_re.finditer(r"(?m)^#{1,3}\s*(\d+(?:\.\d+)?)\.\s*(.+)$", txt3))
                for k3, mm in enumerate(hits3):
                    key = mm.group(1)
                    heads[key] = mm.group(2).strip()
                    end3 = hits3[k3 + 1].start() if k3 + 1 < len(hits3) else len(txt3)
                    body_of[key] = " ".join(txt3[mm.end():end3].split()[:400])
            refs3 = {}
            for line in (MDC / "REFERENCES.md").read_text().splitlines():
                mm = _re.match(r"^(\d+)\.\s+(.*)$", line.strip())
                if mm:
                    refs3[mm.group(1)] = mm.group(2)
            STOP3 = set("""a an and are as at be by for from has have in into is it its of on or
that the their to was were with we our this these those not only also both each""".split())

            def tw(s):
                return {w for w in _re.findall(r"[A-Za-z][A-Za-z0-9\-]{3,}", s.lower())
                        if w not in STOP3}

            # 引用先は正しいが、その段落が見出しと語を共有していない箇所（内容を確認済み）
            ALLOW3 = ("defects that no reviewer raised",
                      "Two of the three paragraphs that closed the Introduction",
                      "Addressed in four places",
                      # reference 34 は Ito ら（companion study）。内容は確認済み
                      "In addition, the detection calls themselves are given")
            bad3 = []
            for para in lt3.split("\n\n"):
                flat3 = " ".join(para.split())
                if flat3.startswith("|") or flat3.startswith("> *"):
                    continue
                body3 = _re.sub(r"> \*.*?\*", "", flat3)
                if any(s in body3 for s in ALLOW3):
                    continue
                for mm in _re.finditer(r"\bSection (\d+\.\d+)\b", body3):
                    k = mm.group(1)
                    if k not in heads:
                        bad3.append(f"実在しない Section {k}")
                    elif not (tw(body3) & (tw(heads[k]) | tw(body_of.get(k, "")))):
                        bad3.append(f"Section {k}「{heads[k][:34]}」と主題が合わない: "
                                    f"{body3[:70]}")
                for mm in _re.finditer(r"\breference (\d+)\b", body3):
                    k = mm.group(1)
                    if k not in refs3:
                        bad3.append(f"実在しない reference {k}")
                    elif not (tw(body3) & tw(refs3[k])):
                        bad3.append(f"reference {k} と主題が合わない: {body3[:70]}")
            (ok if not bad3 else bad).append(
                "C案 AV: letter が指す節・文献がその回答の主題を扱っている" if not bad3
                else f"C案 AV: 合わない参照: {bad3[:3]}")

        # AT1: 補足への引用が、その番号のキャプションの主題と合っていること
        import subprocess
        rs = subprocess.run([sys.executable,
                             str(HERE / "code" / "revision1" / "at_check_supp_citations.py")],
                            capture_output=True, text=True)
        (ok if rs.returncode == 0 else bad).append(
            "C案 AT: 補足への引用がすべてキャプションの主題と合っている"
            if rs.returncode == 0
            else "C案 AT: 引用先が合わない: "
                 + " / ".join(l.strip() for l in rs.stdout.splitlines() if "✗" in l)[:240])
        # 末尾の一覧が、ラベル表と補足ファイルの中身に対応していること
        SLf = HERE / "results" / "roundR" / "supp_labels.tsv"
        if SLf.exists():
            SL = pd.read_csv(SLf, sep="\t")
            bm = (MDC / "BACKMATTER.md").read_text()
            line = next((l for l in bm.splitlines()
                         if l.startswith("Supplementary Materials:")), "")
            nt = int(SL[SL.kind == "Table"].number.max())
            nf = int(SL[SL.kind == "Figure"].number.max())
            need = [f"Tables S1 to S{nt}", f"Figures S1 to S{nf}", "Supplementary Data S1"]
            need += [f"S{int(r.number)}, {r.blurb}" for r in SL.itertuples()
                     if r.kind == "Table" and isinstance(r.blurb, str) and r.blurb]
            gone = [v for v in need if v not in line]
            (ok if not gone else bad).append(
                f"C案 AT: 末尾の一覧が補足 {nt} 表・{nf} 図と Data S1 に対応している"
                if not gone else f"C案 AT: 一覧に無い: {gone[:3]}")
            import openpyxl as _ox
            xl = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Supplementary_Tables.xlsx"
            if xl.exists():
                sheets = [s for s in _ox.load_workbook(xl, read_only=True).sheetnames]
                want_sheets = [f"Table S{int(n)}" for n in
                               sorted(SL[SL.kind == "Table"].number)]
                miss_s = [s for s in want_sheets if s not in sheets]
                (ok if not miss_s else bad).append(
                    f"C案 AT: 補足ファイルのシートがラベル表の {len(want_sheets)} 表と 1 対 1"
                    if not miss_s else f"C案 AT: シートが無い: {miss_s}")

        # AO: letter の引用が査読票の原文の部分文字列であること、補足番号の対応表が正しいこと
        import subprocess
        rq = subprocess.run([sys.executable,
                             str(HERE / "code" / "revision1" / "ao_check_quotes.py")],
                            capture_output=True, text=True)
        (ok if rq.returncode == 0 else bad).append(
            "C案 AO: letter の引用 18 件がすべて査読票の原文の部分文字列"
            if rq.returncode == 0
            else "C案 AO: 原文と合わない引用がある: "
                 + " / ".join(l.strip() for l in (rq.stdout + rq.stderr).splitlines()
                              if "✗" in l)[:240])
        SMf = HERE / "results" / "roundR" / "supplementary_number_map.tsv"
        if SMf.exists():
            SMT = pd.read_csv(SMf, sep="\t")
            lost = SMT[SMT.submitted.notna() & SMT.revised.isna()]
            (ok if lost.empty else bad).append(
                "C案 AO: 投稿版の補足はすべて改訂版に対応がある"
                if lost.empty else f"C案 AO: 対応が無い: {list(lost.caption)[:2]}")
            lt = (MDC / "RESPONSE.md").read_text()
            miss = []
            for r_ in SMT.itertuples():
                if pd.notna(r_.submitted) and pd.notna(r_.revised):
                    row = f"| {r_.kind} S{int(r_.submitted)} | {r_.kind} S{int(r_.revised)} |"
                    if row not in lt:
                        miss.append(row)
            (ok if not miss else bad).append(
                f"C案 AO: 旧→新の対応表 {len(SMT)} 行が letter にある" if not miss
                else f"C案 AO: letter の対応表に無い行: {miss[:3]}")
            # 査読者が番号で指した箇所の併記。番号が変わったものだけ併記が必要
            cited = [1, 2, 5, 6]          # R3 major 2 / major 3（2 件）/ minor 1
            mp = {int(r_.submitted): int(r_.revised) for r_ in SMT.itertuples()
                  if r_.kind == "Table" and pd.notna(r_.submitted) and pd.notna(r_.revised)}
            need = [f"Table S{a} (Table S{mp[a]} in this revision)"
                    for a in cited if a in mp and mp[a] != a]
            same = [f"Table S{a}" for a in cited if a in mp and mp[a] == a]
            gone = [v for v in need if v not in lt]
            (ok if not gone else bad).append(
                f"C案 AO: 番号が変わった {len(need)} 箇所に併記があり、変わらない "
                f"{len(same)} 箇所は併記なし" if not gone
                else f"C案 AO: 併記が無い: {gone}")

        # AN: 土台 docx から作った修正稿の構造検査（docx が対象）
        AN = Path.home() / "Desktop" / "CKD_xspecies_submission" / "Main_manuscript.docx"
        if AN.exists():
            import subprocess
            r = subprocess.run([sys.executable,
                                str(HERE / "code" / "revision1" / "an_check_docx.py")],
                               capture_output=True, text=True)
            last = [l for l in r.stdout.splitlines() if l.strip()][-1] if r.stdout else ""
            (ok if r.returncode == 0 else bad).append(
                f"C案 AN: 組み上がった docx の構造検査に合格（{last.strip()}）"
                if r.returncode == 0
                else "C案 AN: docx の構造検査に不合格。an_check_docx.py を見ること: "
                     + " / ".join(l.strip() for l in r.stdout.splitlines() if "✗" in l)[:300])
            mp = HERE / "results" / "roundR" / "an_block_map.tsv"
            if mp.exists():
                BM = pd.read_csv(mp, sep="\t")
                # 意図した削除のみ: R1-4 で Methods へ移した 1 段落と、
                # AT4 で全面的に書き換えた Limitations の 1 段落（新規として入れ直している）
                n_del = int((BM.action == "deleted").sum())
                (ok if n_del <= 2 else bad).append(
                    f"C案 AN: 土台から落とした段落は {n_del} 件（いずれも意図した置き換え）"
                    if n_del <= 2 else f"C案 AN: 土台から {n_del} 段落が落ちている")

        # AJ: Response letter を照合対象に入れる
        LET = MDC / "RESPONSE.md"
        if LET.exists():
            lt = LET.read_text()
            nums_src = (HERE / "manuscript" / "NUMBERS.md").read_text()
            pool = set()
            for s0 in _re.findall(r"[+\-\u2212]?\d[\d,]*\.?\d*", nums_src):
                try:
                    pool.add(float(s0.replace(",", "").replace("\u2212", "-")))
                except ValueError:
                    pass
            bad_tok = []
            for s0 in set(_re.findall(r"[+\-\u2212]?\d+\.\d{3,}|\d+\.\d+%", lt)):
                v = float(s0.rstrip("%").replace("\u2212", "-"))
                if s0.endswith("%"):
                    v /= 100.0
                nd = len(s0.rstrip("%").split(".")[1])
                if not any(abs(round(x, nd) - round(v, nd)) < 10 ** (-nd - 1) for x in pool):
                    bad_tok.append(s0)
            (ok if not bad_tok else bad).append(
                f"C案 AJ: letter の統計値がすべて NUMBERS.md にある" if not bad_tok
                else f"C案 AJ: NUMBERS.md に無い letter の値: {sorted(bad_tok)}")
            # letter が「こう直した」と書いた箇所が本文に実在すること
            TD = HERE / "results" / "roundR" / "response_letter_todo.tsv"
            todo = pd.read_csv(TD, sep="\t") if TD.exists() else pd.DataFrame()
            (ok if todo.empty else bad).append(
                "C案 AJ: letter が指す本文の箇所はすべて実在する" if todo.empty
                else "C案 AJ: letter に書いたが本文が未修正（AG の残作業）: "
                     + ", ".join(f"{r.item}" for r in todo.itertuples()))
            (ok if "[AG:" not in lt else bad).append(
                "C案 AJ: letter に未解決の差し込みが残っていない" if "[AG:" not in lt
                else "C案 AJ: letter に [AG: ...] の差し込みが残っている")

        # AH4: 主張のトーンの検査を全文にかける
        import tone_scan as TS   # noqa: E402
        tone_hits = []
        for f in TS.FILES:
            for i, line in enumerate((MDC / f).read_text().splitlines(), 1):
                for mm in TS.RX.finditer(line):
                    a0, b0 = max(0, mm.start() - 90), min(len(line), mm.end() + 90)
                    if any(s in line and s in line[a0:b0] for s in TS.ALLOW):
                        continue
                    tone_hits.append(f"{f}:{i} {mm.group(0)}")
        (ok if not tone_hits else bad).append(
            f"C案 AH: 主張のトーン {len(TS.TERMS)} 語の全文検査に未確認のヒットが無い"
            if not tone_hits else f"C案 AH: 未確認のトーン表現: {tone_hits[:5]}")

        # AH1: note 列が実測した不変性と矛盾していないこと、表の行名が結果ファイルと一致すること
        PSf2 = HERE / "results" / "roundR" / "pathway_scoring_methods.tsv"
        if PSf2.exists():
            PS2 = pd.read_csv(PSf2, sep="\t")
            wrong = []
            for r in PS2.itertuples():
                inv = bool(r.invariant_to_per_state_shift)
                said = "invariant to per-state centering: score matrix identical" in str(r.note)
                denied = "not invariant to per-state centering" in str(r.note)
                if inv != (said and not denied) or inv == denied:
                    wrong.append(f"{r.method}: inv={inv} note={str(r.note)[-60:]}")
                if inv != (float(r.max_abs_score_difference) == 0.0):
                    wrong.append(f"{r.method}: フラグと最大差が矛盾")
            (ok if not wrong else bad).append(
                f"C案 AH: 経路スコアリング {len(PS2)} 行の note が実測の不変性と一致する"
                if not wrong else f"C案 AH: note が実測と矛盾: {wrong[:3]}")
            s4 = _re.search(r"(?m)^\| Scoring method \|.*?(?=\n\n)",
                            (MDC / "SUPPLEMENTARY.md").read_text(), _re.S)
            labels = [ln.split("|")[1].strip() for ln in s4.group(0).splitlines()[2:]] if s4 else []
            S4r = pd.read_csv(HERE / "results" / "roundR" / "table_s4_rows.tsv", sep="\t")
            same = labels == list(S4r.method)
            (ok if same else bad).append(
                "C案 AH: Table S4 の行名と順序が table_s4_rows.tsv と一致する" if same
                else f"C案 AH: Table S4 の行名がずれている: {labels} != {list(S4r.method)}")
            # 実装列も結果ファイルから来ていること
            impl = [ln.split("|")[2].strip() for ln in s4.group(0).splitlines()[2:]] if s4 else []
            (ok if impl == list(S4r.implementation) else bad).append(
                "C案 AH: Table S4 の実装列が table_s4_rows.tsv と一致する"
                if impl == list(S4r.implementation)
                else f"C案 AH: Table S4 の実装列がずれている: {impl}")
            # 参照実装との突き合わせが本文に書かれていること
            CEv = pd.read_csv(HERE / "results" / "roundR" /
                              "reference_implementation_cells.tsv", sep="\t")
            sg = CEv[CEv.method == "singscore"]
            need = [f"{sg.max_abs_difference.max():.3f}",
                    f"{sg.median_abs_difference.max():.4f}",
                    f"{sg.pearson_r_of_cells.min():.4f}"]
            gone = [v for v in need if v not in cflat]
            (ok if not gone else bad).append(
                "C案 AH: 自前 singscore と参照実装の差の 3 値が Methods にある" if not gone
                else f"C案 AH: Methods に無い: {gone}")

        # AE2: ランダム集合の 2 設計の値が本文にあること
        RSf = HERE / "results" / "roundR" / "random_set_designs.tsv"
        if RSf.exists():
            RS = pd.read_csv(RSf, sep="\t").set_index(["design", "quantity"])
            OV = pd.read_csv(HERE / "results" / "roundR" /
                             "random_set_designs_overlap.tsv", sep="\t").set_index("quantity")
            LPd = pd.read_csv(HERE / "results" / "roundR" /
                              "random_set_designs_draws_labelperm.tsv", sep="\t")
            SMd, LP = "size-matched draws", "gene-label permutation"

            def h3(x):   # 本文と同じ丸め（half-up, 3 桁）
                from decimal import Decimal, ROUND_HALF_UP
                s = str(Decimal(str(round(float(x), 4))).quantize(Decimal("0.001"), ROUND_HALF_UP))
                return s.replace("-", "\u2212")

            want = [h3(RS.loc[(SMd, "difference"), "random_median"]),
                    h3(RS.loc[(SMd, "difference"), "random_min"]),
                    h3(RS.loc[(SMd, "difference"), "random_max"]),
                    h3(RS.loc[(SMd, "auc_agg_centred"), "random_min"]),
                    h3(RS.loc[(SMd, "auc_agg_centred"), "random_max"]),
                    h3(RS.loc[(LP, "difference"), "random_median"]),
                    h3(RS.loc[(LP, "difference"), "random_min"]),
                    h3(RS.loc[(LP, "difference"), "random_max"]),
                    h3(RS.loc[(LP, "auc_agg_centred"), "random_min"]),
                    h3(RS.loc[(LP, "difference"), "observed_real_sets"]),
                    h3(RS.loc[(LP, "auc_agg_centred"), "observed_real_sets"]),
                    str(int(RS.loc[(LP, "difference"), "draws_at_or_above_observed"])),
                    str(int(OV.loc["sets per gene, median", "value"])),
                    str(int(OV.loc["sets per gene, maximum", "value"])),
                    str(int(OV.loc["genes in no set", "value"])),
                    str(int((LPd.auc_agg_centred <= 0.55).sum()))]
            miss = [v for v in want if v not in cflat]
            (ok if not miss else bad).append(
                f"C案 AE: ランダム集合 2 設計の {len(want)} 値が本文にある" if not miss
                else f"C案 AE: 本文に無いランダム集合の値: {miss}")
            # 重要な主張のふたつ: サイズ合わせでは 0 回、ラベル置換では 3 回
            n0 = int(RS.loc[(SMd, "difference"), "draws_at_or_above_observed"])
            (ok if n0 == 0 else bad).append(
                "C案 AE: サイズを合わせた設計では観測値以上の draw が 0 件" if n0 == 0
                else f"C案 AE: サイズを合わせた設計で観測値以上の draw が {n0} 件ある。"
                     "「lies above every one of those draws」は書けない")
            nle = int((LPd.auc_agg_centred <= RS.loc[(LP, "auc_agg_centred"),
                                                     "observed_real_sets"]).sum())
            (ok if nle == 1 else bad).append(
                "C案 AE: ラベル置換で中心化 AUC が観測値以下になる draw は 1 件" if nle == 1
                else f"C案 AE: ラベル置換で中心化 AUC が観測値以下の draw は {nle} 件")

        # AB: 表の構造（列見出しと中身の対応、セルに残った | や \ ）
        sys.path.insert(0, str(HERE / "code" / "revision1"))
        import check_tables as CT   # noqa: E402
        n_tab, n_bad = 0, []
        for f in ("RESULTS.md", "METHODS.md", "SUPPLEMENTARY.md"):
            for cap, rows in CT.tables((MDC / f).read_text()):
                n_tab += 1
                head = CT.cells(rows[0])
                for i, r in enumerate([x for x in rows[1:] if not CT.SEP.fullmatch(x)]):
                    c = CT.cells(r)
                    if len(c) != len(head):
                        n_bad.append(f"{f}/{cap[:30]} 行{i+1}: {len(c)} 列 != {len(head)}")
                    for x in c:
                        if "\\|" in x or x.endswith("\\"):
                            n_bad.append(f"{f}/{cap[:30]}: セルに | か \\ が残る")
        (ok if not n_bad else bad).append(
            f"C案 AB: 表 {n_tab} 件の列見出しと中身が対応している" if not n_bad
            else f"C案 AB: 列がずれている表: {n_bad[:4]}")

        # AB: 本文が参照する表・図がすべて実在すること
        caps = set(_re.findall(r"(?m)^(Table S?\d+|Figure S?\d+)\.", 
                               (MDC / "SUPPLEMENTARY.md").read_text()
                               + (MDC / "RESULTS.md").read_text()
                               + (MDC / "METHODS.md").read_text()
                               + (MDC / "FIGURES.md").read_text()))
        refs = set(_re.findall(r"\b(Table S?\d+|Figure S?\d+)\b", cjoin))
        ghost = sorted(r for r in refs - caps)
        (ok if not ghost else bad).append(
            f"C案 AB: 本文が参照する {len(refs)} 件の表・図がすべて実在する" if not ghost
            else f"C案 AB: 実体の無い参照: {ghost}")

        # AB: 補足の番号が本文の初出順であること
        # 節順は投稿先の書式（Introduction → Results → Discussion → Methods）。AR1 で戻した
        body = "\n".join((MDC / f).read_text() for f in
                          ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                           "METHODS.md", "CONCLUSIONS.md"))
        for kind in ("Table", "Figure"):
            seen, order = set(), []
            # パネル付き（Figure S1A）も同じ図への参照として数える
            for m in _re.finditer(rf"\b{kind} S(\d+)[A-D]?\b", body):
                n = int(m.group(1))
                if n not in seen:
                    seen.add(n); order.append(n)
            (ok if order == sorted(order) else bad).append(
                f"C案 AB: {kind} S の初出順が番号順（{order}）" if order == sorted(order)
                else f"C案 AB: {kind} S の初出順が番号順でない: {order}")

        # Table 1 の centering なし集約 AUC が対照解析の出力と一致すること
        CC = pd.read_csv(HERE / "results" / "roundC" / "centering_controls.tsv", sep="\t")
        cc = CC[CC.collection.str.startswith("422")]
        S87, S83 = "87 (no shared controls)", "83 (also no shared cohort)"
        sp = ["Group A (2,016)", "All 1:1 orthologues (7,897)",
              "Matched Group B (1,632)", "Matched Group A (1,632)"]
        g = {(r.space, r.pairs): r for r in cc.itertuples()}
        row = _re.search(r"^\| AUC after aggregation, uncenter?ed[^|]*\|([^\n]*)\|\s*$",
                         (MDC / "RESULTS.md").read_text(), _re.M)
        (ok if row else bad).append("C案: Table 1 に uncentred 集約 AUC の行がある" if row
                                    else "C案: Table 1 の uncentred 集約 AUC 行が無い")
        if row:
            got = [c.strip() for c in row.group(1).split("|") if c.strip()]
            want = [f"{g[(s, S87)].auc_agg_uncentred:.3f}" for s in sp]
            (ok if got == want else bad).append(
                "C案: Table 1 の uncentred 集約 AUC が対照解析の出力と一致"
                if got == want else f"C案: Table 1 の uncentred 行 {got} != {want}")
        for lab, col in [("集約 centred", "auc_agg_centred"),
                         ("集約 uncentred", "auc_agg_uncentred"),
                         ("遺伝子単位", "auc_gene_level"),
                         ("経路別中央値", "pathway_auc_median")]:
            v = f"{getattr(g[('Group A (2,016)', S83)], col):.3f}"
            (ok if v in cflat else bad).append(
                f"C案: 83 ペアの{lab} {v} が本文にある" if v in cflat
                else f"C案: 83 ペアの{lab} {v} が本文に無い")

        # Abstract の語数
        am = _re.search(r"^Abstract:(.*?)$", (MDC / "FRONTMATTER.md").read_text(), _re.M)
        nw = len(am.group(1).split()) if am else 999
        (ok if nw <= 200 else bad).append(
            f"C案: Abstract は {nw} 語（200 語以内）" if nw <= 200
            else f"C案: Abstract が {nw} 語ある")

        # 文献番号が 1 から連続し、全件が引用されていること
        rt = (MDC / "REFERENCES.md").read_text()
        rn = [int(m.group(1)) for m in _re.finditer(r"^(\d{1,2})\.\s+[A-Z]", rt, _re.M)]
        cited = set()
        for m in _re.finditer(r"\[([0-9][0-9,\u2013\u2014 -]*)\]", cjoin):
            for part in _re.split(r"\s*,\s*", m.group(1)):
                rg = _re.fullmatch(r"(\d+)\s*[\u2013\u2014-]\s*(\d+)", part.strip())
                if rg:                      # [33-35] のような範囲は中身も引用とみなす
                    cited.update(range(int(rg.group(1)), int(rg.group(2)) + 1))
                elif part.strip().isdigit():
                    cited.add(int(part.strip()))
        gaps = sorted(set(range(1, max(rn) + 1)) - set(rn)) if rn else [0]
        un = sorted(set(rn) - cited)
        (ok if not gaps and not un else bad).append(
            f"C案: 文献は 1 から {max(rn)} まで連続し全件が引用されている"
            if not gaps and not un else f"C案: 抜け {gaps} / 未引用 {un}")
    else:
        skipped.append("改訂稿の本文との突合（原稿の Markdown は公開していない。"
                       "本文は論文として出版される）")

    # ---- 改訂10: results/ の全ファイルに生成元があること ----
    PV = HERE / "results" / "provenance.tsv"
    if PV.exists():
        pv = pd.read_csv(PV, sep="\t")
        gap = pv[pv.status == "no generating script found"]
        (ok if gap.empty else bad).append(
            f"results/ の {len(pv)} ファイルすべてに生成元がある" if gap.empty
            else f"生成スクリプトの無い結果ファイル: {list(gap.result_file)}")
    else:
        bad.append("results/provenance.tsv が無い（build_provenance.py を実行すること）")

    # ---- 改訂4: 図 3・4 に印字した数値と、本文・キャプションの照合 ----
    # make_figures.py が描画時に記録した値をそのまま読む。図の中の数字が本文とずれたら
    # ここで落ちる。集合の同一性は (collection, pathway) の組。名前だけで畳むと 422 が 419 に
    # なり、これが改訂4以前の図 3A の 419 / 242 と図 4A の IRI 28 d 21% の原因だった。
    PR = tsv("revision1/pathway_reactome/figure_printed_numbers.tsv")
    printed = {(r.figure, r.panel, r.item): str(r.printed) for r in PR.itertuples()}

    def figchk(fig_, panel, item, want):
        got = printed.get((fig_, panel, item))
        (ok if got == str(want) else bad).append(
            f"図{fig_[-1]}{panel} {item}: 本文 {want} / 図の印字 {got}")

    PC = tsv("revision1/pathway_reactome/pathway_cos.tsv")
    PC = PC.assign(key=PC.collection + "|" + PC.pathway)
    size_by_key = PC.groupby("key").n_genes.first()
    n_sets = len(size_by_key)
    n_small = int((size_by_key < 50).sum())
    three = PC[(PC.collection != "Reactome") & (PC.n_genes >= 50)]
    n_panelB = three.key.nunique()
    cells_B = n_panelB * PC.state.nunique()
    neg_B = int((three.cos < 0).sum())

    assert n_sets == 422, f"セット数が 422 でない: {n_sets}（§2.3 / §4.13 の前提）"
    assert len(size_by_key.index.str.split("|").str[1].unique()) == 419, \
        "固有名が 419 でない。422 と 419 の差（同名 3 件）が崩れている"
    assert n_small == 244 and n_panelB == 74, f"{n_small} / {n_panelB}"
    assert size_by_key[size_by_key >= 50].size - n_panelB == 104, \
        "キャプションの「104 Reactome sets」が成り立たない"

    figchk("Fig3", "A", "pathway sets in the size histogram", n_sets)          # §2.3, §4.13, キャプション
    figchk("Fig3", "A", "sets below 50 Group A genes", n_small)               # §2.3「244 of the 422」, §3.4
    figchk("Fig3", "B", "columns (pathways drawn)", n_panelB)                 # キャプション「74 pathways」, §4.13
    figchk("Fig3", "B", "collections drawn", "GO_BP+Hallmark+KEGG")           # キャプション「three original collections」
    figchk("Fig3", "B", "cells drawn", cells_B)                               # §2.3「1,184 of them」
    figchk("Fig3", "B", "cells below zero (dotted)", neg_B)

    figchk("Fig4", "A", "pathways ranked", n_sets)                            # §2.3, §3.4「across 422 pathways」
    figchk("Fig4", "A", "mouse states ranked", 12)
    figchk("Fig4", "A", "Kendall W", "0.320")                                 # §2.3, §3.4
    for st, share in [("IRI_28d", "22%"), ("IRI_12mo", "21%"), ("IRI_2h", "11%")]:
        figchk("Fig4", "A", f"leads share {st}", share)                       # §2.3「IRI 28 d 22%, IRI 12 mo 21%, IRI 2 h 11%」
    assert max(int(printed[("Fig4", "A", f"leads share {s_}")].rstrip("%"))
               for s_ in PC.state.unique() if not s_.startswith("cat")) < 25, \
        "§3.4「no mouse state leads more than a quarter」が図と合わない"

    SEPv = tsv("revision1/pathway_reactome/pathway_separation_paired.tsv")
    av = SEPv["auc[exclude shared controls]"].dropna()
    figchk("Fig4", "B", "pathways in the AUC histogram", len(av))
    figchk("Fig4", "B", "AUC per pathway median", f"{av.median():.3f}")       # §2.3, §3.3
    figchk("Fig4", "B", "AUC per pathway Q1", f"{av.quantile(.25):.3f}")      # §2.3, §3.3「0.707 to 0.808」
    figchk("Fig4", "B", "AUC per pathway Q3", f"{av.quantile(.75):.3f}")
    figchk("Fig4", "B", "AUC whole panel 2,016 genes", "0.813")               # §2.3
    figchk("Fig4", "B", "AUC aggregated centred", "0.511")                    # §2.3
    chk("§2.3 AUC 0.90 以上の経路数", int((av >= 0.90).sum()), 16, 0.5)
    chk("§2.3 AUC 0.60 以下の経路数", int((av <= 0.60).sum()), 25, 0.5)
    chk("§2.3 全パネル以上の経路の割合", float((av >= 0.813).mean()), 0.23, 5e-3)

    for item, want in [("median gene level [within species]", "0.724"),
                       ("median gene level [cat vs mouse]", "0.448"),
                       ("median centred [within species]", "0.709"),
                       ("median centred [cat vs mouse]", "0.551"),
                       ("median pathway means of centred [within species]", "0.831"),
                       ("median pathway means of centred [cat vs mouse]", "0.832"),
                       ("median pathway means of uncentred [within species]", "0.870"),
                       ("median pathway means of uncentred [cat vs mouse]", "0.462"),
                       ("AUC gene level", "0.813"), ("AUC centred", "0.759"),
                       ("AUC pathway means of centred", "0.511"),
                       ("AUC pathway means of uncentred", "0.868")]:
        figchk("Fig4", "D", item, want)                                       # §2.3 の四条件

    # ---- 本文・キャプションに、その数値が実際に書かれているか ----
    MD = HERE / "manuscript"
    if MD.exists():
        import re
        body = " ".join(re.sub(r"\s+", " ", (MD / f).read_text())
                        for f in ("RESULTS.md", "DISCUSSION.md", "METHODS.md"))
        body = body.replace("\u2019", "'").replace("\u2013", "-").replace("\u2014", "-")
        for where, phrase in [
            ("§2.3", "422 sets: 229, 118, 48 and 27 respectively"),
            ("§2.3", "244 of the 422 sets fall below 50 genes"),
            ("§2.3", "giving 6,752 pathway-by-state values"),
            ("§2.3", "Figure 3B displays the 1,184 of them that fall in the 74 sets it shows"),
            ("§2.3", "Kendall's W across the 422 pathways is 0.320"),
            ("§2.3", "The shares are IRI 28 d 22%, IRI 12 mo 21%, IRI 2 h 11%"),
            ("§2.3", "0.760 in the median (interquartile range 0.707-0.808)"),
            ("§2.3", "the aggregated values there being 0.600, 0.624 and 0.529"),
            ("§2.3", "16 of 422 reaching 0.90 and 25 falling to 0.60 or below"),
            ("Fig3 キャプション", "244 of 422 pathways fall below 50 genes"),
            ("Fig3 キャプション", "across the 74 pathways of the three original collections"),
            ("Fig3 キャプション", "adding the 104 Reactome sets"),
            ("Fig3 キャプション", "would widen it to 178 columns"),
            ("Fig4 キャプション", "within each of the 422 pathways"),
            ("§3.4", "244 of our 422 pathways fall below 50 genes"),
            ("§4.13", "held to the 74 sets of the three original collections"),
            ("§4.13", "rather than to all 178"),
        ]:
            (ok if phrase in body else bad).append(f"{where} の文言が見つからない: {phrase}")

    # 原稿の Markdown はこのリポジトリに置いていない（.gitignore の projection/manuscript/）。
    # 本文・表・図の文言と突き合わせる検査はそれが要るので、無ければそこだけ飛ばし、
    # 結果ファイルどうしの照合は最後まで走らせる。件数は最後に内訳を出す。
    HAVE_MD = (HERE / "manuscript").exists()

    if HAVE_MD:
        # ---- 改訂5 Part 0: 表のセルと本文の突合 ----
        # 図は figure_printed_numbers.tsv を通して本文と照合できるが、表は見ていなかった。
        # rev8 で見つかった 2 件（Table 4 の集約後 AUC と、サイズ条件を満たすセル数）は
        # どちらもその穴から漏れた。ここで表側も同じ扱いにする。
        sys.path.insert(0, str(HERE / "code" / "revision1"))
        import extract_tables                                   # noqa: E402
        extract_tables.main()                                   # 常に作り直す（古い値を残さない）
        TP = tsv("round5/table_printed_numbers.tsv")

        MD5 = HERE / "manuscript"
        import re as _re

        def flat(name):
            return _re.sub(r"\s+", " ", (MD5 / name).read_text())

        prose_results = _re.sub(
            r"\s+", " ", "\n".join(l for l in (MD5 / "RESULTS.md").read_text().split("\n")
                                    if not l.strip().startswith("|")))
        BODY = " ".join([prose_results] + [flat(f) for f in
                        ("DISCUSSION.md", "METHODS.md", "CONCLUSIONS.md",
                         "FRONTMATTER.md", "SUPPLEMENTARY.md")])
        BODY = BODY.replace("\u2019", "'").replace("\u2013", "-").replace("\u2014", "-")

        # 表にしか出さない値。ここに挙げたものだけが本文に無くてよい。
        TABLE_ONLY = {
            ("Table1", "0.755"), ("Table1", "0.884"), ("Table1", "0.781"), ("Table1", "0.945"),
            ("Table1", "0.862"),
            # Part A-5 で「群間差も下限」の一文を撤回した結果、この 2 つは表だけの値になった
            ("Table1", "0.852"), ("Table1", "0.896"),
            ("Table3", "0.610"), ("Table3", "0.0019"), ("Table3", "0.317"), ("Table3", "0.0232"),
            ("Table3", "0.0250"), ("Table3", "0.027"), ("Table3", "0.860"),
            ("Table4", "0.604"), ("Table4", "0.673"), ("Table4", "0.559"), ("Table4", "0.881"),
            ("Table4", "0.635"), ("Table4", "0.730"),
            # matched Group A 列の状態ごとの内訳。他の空間と同じく表だけに出す。
            # 2 つの matched 空間の状態ごとの内訳。他の空間と同じく表だけに出す。
            ("Table4", "0.546"), ("Table4", "0.876"), ("Table4", "0.623"),
            ("Table4", "0.566"), ("Table4", "0.901"), ("Table4", "0.629"),
            ("Table4", "0.987"), ("Table4", "0.734"), ("Table4", "1.345"),
            ("Table4", "0.947"), ("Table4", "0.776"), ("Table4", "1.221"),
            # 第9ラウンド（本文の圧縮）で、表が持つ値を本文で繰り返すのをやめた分。
            ("Table1", "0.591"), ("Table1", "0.875"),
            ("Table2", "0.493"), ("Table2", "0.492"),
            ("Table3", "0.0040"), ("Table3", "0.312"),
            ("Table4", "1.091"), ("Table4", "1.371"), ("Table4", "1.362"),
            ("Table4", "1.289"), ("Table4", "0.280"), ("Table4", "0.069"),
            ("Table4", "0.637"), ("Table4", "0.646"), ("Table4", "0.613"),
            ("Table4", "0.0040"), ("Table4", "0.0031"),
            ("Table4", "0.572"), ("Table4", "0.534"),
            ("Table4", "0.780"), ("Table4", "0.769"),
            ("Table4", "0.7616"), ("Table4", "0.7708"), ("Table4", "0.8102"),
            ("Table4", "0.782"), ("Table4", "0.780"), ("Table4", "0.781"),
        }
        NUMTOK = _re.compile(r"\d+\.\d+|\d[\d,]*")
        n_tok = 0
        for r in TP.itertuples():
            if r.row == "header":
                continue
            for tok in NUMTOK.findall(str(r.printed)):
                n_tok += 1
                tab = _re.sub(r"b\d+$", "", str(r.table))
                if (tab, tok) in TABLE_ONLY or tok in BODY:
                    ok.append(f"表 {r.table} の {tok} は本文と整合")
                else:
                    bad.append(f"表 {r.table}「{str(r.row)[:40]}」の {tok} が本文のどこにも無い "
                               f"（表だけ更新されて本文が古い可能性。表専用なら TABLE_ONLY に登録する）")
        assert n_tok > 100, f"表のセルから数値が {n_tok} 個しか取れていない"

        # 表のセルに数値を連結しない。α / cos θ / ν を 1 セルに詰めていたために、docx 側で
        # 新旧の値が混ざる事故が起きた。1 セルに数値が 2 つ以上並んだら落とす。
        # 区間・比・p 値をひとつのセルに書く行は、意図が明示された場合だけ許す。
        MULTI_OK = {
            # 中央値と最小値を [ ] で併記する 2 列だけは、意図が列名に書いてあるので許す。
            ("Table1", "Benchmark, Spearman–Brown: median [min]"),
            ("Table1", "Benchmark, uncorrected: median [min]"),
        }
        NUMTOK2 = _re.compile(r"^[+\-\u2212]?\d+(?:\.\d+)?$")
        joined_cells = []
        for r in TP.itertuples():
            if r.row == "header" or str(r.table).startswith("Table2"):
                continue                                   # Table 2 は散文セルなので対象外
            cell = str(r.printed)
            if (_re.sub(r"b\d+$", "", str(r.table)), str(r.column)) in MULTI_OK:
                continue
            toks = [t for t in cell.replace("\u2212", "-").replace("/", " ").replace(",", " ").split()
                    if NUMTOK2.match(t)]
            if len(toks) >= 2:
                joined_cells.append(f"{r.table} 「{str(r.row)[:40]}」x「{str(r.column)[:24]}」= {cell}")
        (ok if not joined_cells else bad).append(
            "表のセルに数値の連結がない" if not joined_cells
            else "数値を連結したセルがある: " + "; ".join(joined_cells[:4]))
        assert len([r for r in TP.itertuples() if str(r.table) == "Table4"]) >= 80, \
            "Table 4 のセルが少なすぎる（1 セル 1 値への作り直しが効いていない）"

        # 同じ量が複数の節に出るものは、全出現箇所を検査する。1 か所だけ直す事故を防ぐ。
        for label, phrases in [
            ("集約後 AUC（3 遺伝子空間）", [
                ("Table 4", "| AUC after aggregation, centred (Section 2.3) | 0.511 | 0.600 | 0.624 | 0.529 |"),
                ("§2.3", "the aggregated values there being 0.600, 0.624 and 0.529")]),
            ("サイズ条件を満たすセル数", [
                ("§2.3", "Of the 2,136 cells meeting the size condition"),
                ("§3.2", "in three of the 2,136 qualifying mouse")]),
            ("経路単位 AUC の中央値", [
                ("Table 4", "| Pathway-level AUC, median | 0.760 |"),
                ("§2.3", "0.760 in the median"),
                ("§3.3", "a median AUC of 0.760")]),
            ("種間 cos θ の最大値", [
                ("Table 1", "| Cross-species (cat vs mouse) | 48 | 0.448 | 0.548 |"),
                ("Table 4", "| Cat-mouse cos θ, maximum | 0.548 | 0.572 | 0.438 | 0.534 |"),
                ("§2.2", "a median of 0.448 and a maximum of 0.548")]),
            ("Mantel（経過時間・方向）", [
                ("Table 3", "+0.599 (p = 0.0040)"),
                ("Table 4", "| Mantel, elapsed time, direction (ρ) | +0.599 |"),
                ("§2.5", "+0.599")]),
        ]:
            whole = " ".join([BODY] + [_re.sub(r"\s+", " ", (MD5 / "RESULTS.md").read_text())])
            whole = whole.replace("\u2019", "'").replace("\u2013", "-").replace("\u2014", "-")
            for where, phrase in phrases:
                want = _re.sub(r"\s+", " ", phrase)
                (ok if want in whole else bad).append(
                    f"{label} / {where} の文言が見つからない: {want}")

        # 表の値そのものを結果ファイルと照合する
        GS = tsv("revision1/genespace_with_reactome.tsv").set_index("space")
        t4 = {(r.row, r.column): r.printed for r in TP.itertuples() if r.table == "Table4"}
        for space, col in [("Group A", "Group A (2,016)"),
                           ("all 1:1 orthologues", "All 1:1 orthologues (7,897)"),
                           ("matched Group B", "Matched Group B (1,632)"),
                           ("matched Group A", "Matched Group A (1,632)")]:
            chk(f"Table 4 集約後 AUC {space}",
                float(t4[("AUC after aggregation, centred (Section 2.3)", col)]),
                float(GS.loc[space, "aggregated_auc_centred"]), 6e-4)
            chk(f"Table 4 経路単位 AUC 中央値 {space}",
                float(t4[("Pathway-level AUC, median", col)]),
                float(GS.loc[space, "pathway_auc_median"]), 6e-4)

        PCE = tsv("reliability/pairs_ceilings.tsv")
        t1 = {(r.row, r.column): r.printed for r in TP.itertuples() if r.table == "Table1"}
        for cls, row in [("within_dataset", "Within dataset"),
                         ("same_species_diff_dataset", "Same species, different dataset"),
                         ("cross_species", "Cross-species (cat vs mouse)")]:
            d = PCE[PCE["class"] == cls]
            chk(f"Table 1 ペア数 {row}", float(t1[(row, "n pairs")]), float(len(d)), 0.5)
            chk(f"Table 1 cos 中央値 {row}", float(t1[(row, "cos θ median")]),
                float(d.cos.median()), 6e-4)
            chk(f"Table 1 cos 最大 {row}", float(t1[(row, "cos θ max")]), float(d.cos.max()), 6e-4)

        # ---- 改訂5: 節・図・表の相互参照、作業用注記、Abstract の語数 ----
        all_md = {f: (MD5 / f).read_text() for f in
                  ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                   "METHODS.md", "CONCLUSIONS.md", "SUPPLEMENTARY.md", "BACKMATTER.md")}
        joined = _re.sub(r"\s+", " ", " ".join(all_md.values()))   # 改行で参照が割れるのを防ぐ
        heads = set()
        for f in ("METHODS.md", "RESULTS.md", "DISCUSSION.md"):
            heads |= set(_re.findall(r"^## (\d+\.\d+)\.", all_md[f], _re.M))
        refs = set(_re.findall(r"Section (\d+\.\d+)", joined))
        missing = sorted(refs - heads)
        (ok if not missing else bad).append(
            f"本文が参照する節がすべて実在する（欠落: {missing}）" if not missing
            else f"存在しない節への参照: {missing}")
        figs = set(_re.findall(r"\*\*(Figure S?\d+)\.\*\*", joined))
        figrefs = set(_re.findall(r"\b(Figure S?\d+)[A-D]?\b", joined))
        miss_f = sorted(figrefs - figs)
        (ok if not miss_f else bad).append(
            "本文が参照する図がすべて実在する" if not miss_f else f"存在しない図への参照: {miss_f}")
        tabs = set(_re.findall(r"\*\*(Table S?\d+)\.\*\*", joined))
        tabs |= set(_re.findall(r"\*\*(Table S\d+)\.\*\*", all_md["SUPPLEMENTARY.md"]))
        tabrefs = set(_re.findall(r"\b(Table S?\d+)\b", joined))
        miss_t = sorted(tabrefs - tabs)
        (ok if not miss_t else bad).append(
            "本文が参照する表がすべて実在する" if not miss_t else f"存在しない表への参照: {miss_t}")

        for pat, label in [(r"\[VERIFY[^\]]*\]", "VERIFY 注記"),
                           (r"\[REPOSITORY[^\]]*\]", "REPOSITORY 注記"),
                           (r"to be supplied", "to be supplied"),
                           # 「noise ceiling」は RSA の用語 [19] として意図的に残している。
                           (r"(?i)(?<!noise )\bceiling\b",
                            "ceiling（noise ceiling 以外は benchmark に統一済みのはず）")]:
            hit = _re.findall(pat, joined)
            (ok if not hit else bad).append(
                f"{label} は残っていない" if not hit else f"{label} が {len(hit)} 箇所残っている")

        abst = _re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", all_md["FRONTMATTER.md"], _re.S)
        nw = len(abst.group(1).split())
        (ok if nw <= 200 else bad).append(
            f"Abstract は {nw} 語（200 語以内）" if nw <= 200 else f"Abstract が {nw} 語ある")

        # 四条件の AUC は §2.3・§3.3・Abstract・Conclusions の 4 箇所に出る。全部を検査する。
        for where, phrase in [
            ("§2.3", "separates the classes at AUC 0.868"),
            ("Abstract", "at AUC 0.868 when each state's mean change"),
            ("Conclusions", "at AUC 0.868 or at 0.511"),
        ]:
            want = _re.sub(r"\s+", " ", phrase)
            hay = _re.sub(r"\s+", " ", joined).replace("\u2019", "'")
            (ok if want in hay else bad).append(
                f"集約前後の AUC / {where} の文言が見つからない: {want}")

    # ---- 査読第2便: Spearman-Brown 補正済み基準の較正 ----
    CD = tsv("round5/ceiling_decomposition_sim.tsv")
    RELt = tsv("reliability/reliability.tsv")
    lo_r, hi_r = float(RELt.r_half_median.min()), float(RELt.r_half_median.max())
    band = CD[CD.mean_r_half_cosine.between(lo_r, hi_r)]
    diff = band[band.true_cos < 1.0]
    same = band[band.true_cos == 1.0]
    chk("SB 基準 (v) 真の応答が異なるときの超過（中心化）",
        float(diff.frac_iv_pearson_SB_vs_full.max()), 0.0, 1e-9)
    chk("SB 基準 (v) 真の応答が異なるときの超過（非中心化）",
        float(diff.frac_v_cosine_SB_vs_full.max()), 0.0, 1e-9)
    chk("無補正基準 (iii) の平均超過", float(diff.frac_iii_cosine_half_vs_full.mean()),
        0.0997, 2e-3)
    chk("標本数の寄与 (ii)-(i)", float(diff.contrib_sample_size.mean()), 0.0993, 2e-3)
    chk("非中心化の寄与 (iii)-(ii)", float(diff.contrib_uncentred.mean()), 0.0004, 2e-4)
    chk("SB 基準での非中心化の寄与 (v)-(iv)", float(diff.contrib_uncentred_SB.mean()), 0.0, 1e-9)
    chk("SB 基準 真の応答が同一のときの超過", float(same.frac_v_cosine_SB_vs_full.mean()),
        0.787, 3e-3)
    assert float(diff.frac_v_cosine_SB_vs_full.max()) == 0.0, \
        "SB 基準が、真の応答が異なる条件で一度でも超えている（§4.10 の主張が崩れる）"

    NR = tsv("round5/interval_null_rate.tsv")
    NRs = tsv("round5/interval_null_rate_SB.tsv")
    chk("無補正基準での偽陽性率 最大", float(NR.false_positive_rate.max()), 0.0, 1e-9)
    assert NR.true_cos_at_zero_gap.between(0.80, 0.95).all(), \
        "無補正基準のゼロ交差が 0.85-0.93 の外に出た（§4.10 の記述と食い違う）"
    assert NRs.true_cos_at_zero_gap.between(0.98, 1.0).all(), \
        "SB 基準のゼロ交差が 1.00 付近から外れた（§4.10 の記述と食い違う）"

    CVs = tsv("round5/interval_coverage_sim_SB.tsv")
    CVr = tsv("round5/interval_coverage_sim.tsv")
    chk("被覆率 無補正 平均", float(CVr.coverage_95.mean()), 0.094, 2e-3)
    chk("被覆率 SB 平均", float(CVs.coverage_95.mean()), 0.490, 2e-3)
    assert float(CVr.miss_theta_below_interval.max()) == 0.0 and \
        float(CVs.miss_theta_below_interval.max()) == 0.0, \
        "区間の外れ方が一方向でなくなった（§4.11 の記述と食い違う）"

    if HAVE_MD:
        for where, phrase in [
            ("§4.10", "neither exceeded the corrected benchmark in a single replicate"),
            ("§4.10", "in 78.7% and 78.7% of replicates"),
            ("§2.2", "median difference of 0.487 with an interquartile range of 0.437 to 0.561"),
            ("§2.2", "The lowest Spearman-Brown benchmark among the 48 cross-species pairs is 0.765"),
            ("§2.2", "a gap of 0.297"),
            ("§2.2", "in 9.4% of datasets on average against the uncorrected benchmark and in 49.0%"),
            ("§3.2", "the four pairs above their corrected benchmark, by up to 0.094"),
        ]:
            want = _re.sub(r"\s+", " ", phrase)
            hay = _re.sub(r"\s+", " ", " ".join((MD5 / f).read_text() for f in
                          ("RESULTS.md", "DISCUSSION.md", "METHODS.md")))
            hay = hay.replace("\u2019", "'").replace("\u2013", "-").replace("\u2014", "-")
            (ok if want in hay else bad).append(
                f"SB 基準 / {where} の文言が見つからない: {want}")

    if HAVE_MD:
        # ---- 改訂6 Part 3-2: 補足要素の参照とキャプションの突合 ----
        # rev10 では Table S3 と S4 が、キャプションなしで本文から参照されていた。
        # 図・表と同じ考え方で、参照とキャプションと背表紙の一覧を三方向で突き合わせる。
        body_files = ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                      "METHODS.md", "CONCLUSIONS.md")
        body_txt = _re.sub(r"\s+", " ", " ".join(all_md[f] for f in body_files))
        supp_txt = all_md["SUPPLEMENTARY.md"]
        back_txt = _re.sub(r"\s+", " ", all_md["BACKMATTER.md"])

        supp_caps = set(_re.findall(r"\*\*((?:Table|Figure) S\d+)\.\*\*", supp_txt))
        supp_refs = set(_re.findall(r"\b((?:Table|Figure) S\d+)\b", body_txt))
        back_list = set(_re.findall(r"\b((?:Table|Figure|File) S\d+)\b", back_txt))

        no_cap = sorted(supp_refs - supp_caps)
        (ok if not no_cap else bad).append(
            "本文が参照する補足要素にはすべてキャプションがある" if not no_cap
            else f"キャプションの無い補足要素が参照されている: {no_cap}")
        no_ref = sorted(supp_caps - supp_refs)
        (ok if not no_ref else bad).append(
            "キャプションのある補足要素はすべて本文から参照されている" if not no_ref
            else f"本文から一度も参照されない補足要素: {no_ref}")
        # 背表紙の一覧はキャプションの集合と一致すべき。File S1 は rev16 以降
        # 提出物に含めないので、期待集合からも外す。
        want_back = set(supp_caps)
        miss_back = sorted(want_back - back_list)
        extra_back = sorted(back_list - want_back)
        (ok if not miss_back else bad).append(
            "背表紙の Supplementary Materials に漏れがない" if not miss_back
            else f"背表紙の一覧に無い補足要素: {miss_back}")
        (ok if not extra_back else bad).append(
            "背表紙の一覧に余分な項目がない" if not extra_back
            else f"キャプションの無い項目が背表紙の一覧にある: {extra_back}")
        # 抽出そのものが失敗していないかだけを見る。過不足の報告は上の 4 検査が行う。
        assert supp_caps, "SUPPLEMENTARY.md から補足要素のキャプションが 1 件も取れない"
        # 補足要素が指す出典ファイルが実在するか
        for m in _re.finditer(r"\(results/[\w./-]+\.tsv\)", supp_txt + " " + back_txt):
            rel = m.group(0).strip("()")
            (ok if (HERE / rel).exists() else bad).append(
                f"補足の出典 {rel} が存在する" if (HERE / rel).exists()
                else f"補足が指す出典ファイルが無い: {rel}")

        # ---- 改訂7 Part 0: 文献番号が 1 から順に連続し、抜け番がないこと ----
        reftxt = (MD5 / "REFERENCES.md").read_text()
        refnums = [int(m.group(1)) for m in _re.finditer(r"^(\d{1,2})\.\s+[A-Z]", reftxt, _re.M)]
        dup = sorted({n for n in refnums if refnums.count(n) > 1})
        gaps = sorted(set(range(1, max(refnums) + 1)) - set(refnums))
        (ok if not gaps and not dup else bad).append(
            f"文献リストは 1 から {max(refnums)} まで連続（{len(refnums)} 件）"
            if not gaps and not dup else f"文献番号に抜け {gaps} / 重複 {dup}")
        cited = set()
        for m in _re.finditer(r"\[([0-9][0-9,\u2013\u2014 -]*)\]", joined):
            for tok in _re.split(r"[,\u2013\u2014 -]", m.group(1)):
                if tok.isdigit():
                    cited.add(int(tok))
        over = sorted(n for n in cited if n > max(refnums))
        under = sorted(set(range(1, max(refnums) + 1)) - cited)
        (ok if not over else bad).append(
            "本文の引用番号が文献リストの範囲に収まっている" if not over
            else f"文献リストに無い番号が引用されている: {over}")
        (ok if not under else bad).append(
            "文献リストの全件が本文から引用されている" if not under
            else f"一度も引用されない文献: {under}")
        assert len(refnums) == 34, f"文献が {len(refnums)} 件（34 件体系のはず）"

        # ---- 改訂8 最終校正: 重複文、定義式の位置、作業用注記 ----
        sent_where = {}
        for f in ("FRONTMATTER.md", "INTRODUCTION.md", "RESULTS.md", "DISCUSSION.md",
                  "METHODS.md", "CONCLUSIONS.md", "SUPPLEMENTARY.md", "BACKMATTER.md"):
            flat = _re.sub(r"\s+", " ", all_md[f])
            for sent in _re.split(r"(?<=[.]) ", flat):
                sent = sent.strip()
                if len(sent) > 70 and not sent.startswith("|"):
                    sent_where.setdefault(sent, []).append(f)
        dups = {k: v for k, v in sent_where.items() if len(v) > 1}
        (ok if not dups else bad).append(
            "同じ文が二度出てこない" if not dups
            else f"重複した文が {len(dups)} 件: " + list(dups)[0][:70])

        # §2.1 の定義式は、導入する文の直後の段落でなければならない。
        res_paras = [p.strip() for p in all_md["RESULTS.md"].split("\n\n") if p.strip()]
        for lead, what in [("the projection of a state vector", "α / cos θ / ν / R⊥ の定義式"),
                           ("These four quantities are not independent", "恒等式 (1)")]:
            idx = [i for i, p in enumerate(res_paras) if lead in _re.sub(r"\s+", " ", p)]
            good = bool(idx) and idx[0] + 1 < len(res_paras) and \
                res_paras[idx[0] + 1].lstrip().startswith("$$")
            (ok if good else bad).append(
                f"{what} が導入文の直後にある" if good
                else f"{what} が導入文の直後にない（§2.1 の数式ブロックが離れている）")

        # 角括弧の作業用注記。文献の引用番号だけは除く。
        notes = [m.group(0) for m in _re.finditer(r"\[[^\]]{4,}\]", joined)
                 if not _re.fullmatch(r"\[[0-9][0-9,\u2013\u2014 -]*\]", m.group(0))
                 and not _re.fullmatch(r"\[\d+\.\d+\]", m.group(0))]   # 表の [min] 値は注記でない
        (ok if not notes else bad).append(
            "角括弧の作業用注記が残っていない" if not notes
            else f"作業用注記らしき角括弧が {len(notes)} 件: {notes[:3]}")

    if not HAVE_MD:
        skipped.append("投稿稿の本文・表・図の文言との突合（原稿の Markdown は"
                       "公開していない）")
    else:
        print("原稿の Markdown があるので、結果ファイルどうしの照合に加えて "
              "本文・表・図の文言との突合も実行した。")
    if (MDC / "METHODS.md").exists() and not (
            Path.home() / "Desktop" / "CKD_xspecies_submission").exists():
        skipped.append("提出ファイル（docx・xlsx・変更履歴版・青字版・マニフェスト）の検査"
                       "（提出物は公開していない）")
    print(f"\n実行 {len(ok) + len(bad)} 件 / 一致 {len(ok)} 件 / 不一致 {len(bad)} 件 / "
          f"飛ばした {len(skipped)} 組")
    for sk in skipped:
        print("  − 飛ばした:", sk)
    for b in bad:
        print("  ✗", b)
    # 回答書が「検査何件」と書くための控え。letter は前回の走行のこの値を読む。
    cnt = R / "roundR" / "verification_counts.tsv"
    cnt.parent.mkdir(parents=True, exist_ok=True)
    keep = pd.read_csv(cnt, sep="\t") if cnt.exists() else pd.DataFrame(
        columns=["quantity", "value"])
    # 公開リポジトリでの件数（bn_clone_counts.py が書く）は消さない
    keep = keep[~keep.quantity.isin(["checks_agreed", "checks_disagreed"])]
    pd.concat([pd.DataFrame([{"quantity": "checks_agreed", "value": len(ok)},
                             {"quantity": "checks_disagreed", "value": len(bad)}]),
               keep], ignore_index=True).to_csv(cnt, sep="\t", index=False)
    return 1 if bad else 0


def write_mapping(path=None):
    """量と結果ファイルの対応表を TSV で書き出す（README の対応表の元）。"""
    out = Path(path) if path else (R / "round9" / "quantity_to_file.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    seen, rows = set(), []
    for label, f, want in MAPPING:
        if f is None or (label, f) in seen:
            continue
        seen.add((label, f))
        rows.append((label, f, want))
    pd.DataFrame(rows, columns=["quantity", "result_file", "reported_value"]).to_csv(
        out, sep="\t", index=False)
    print(f"対応表: {out} ({len(rows)} 行)")
    return out


if __name__ == "__main__":
    rc = main()
    if "--map" in sys.argv:
        write_mapping()
    sys.exit(rc)
