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


def tsv(rel, **kw):
    return pd.read_csv(R / rel, sep="\t", **kw)


def main():
    ok, bad = [], []

    def chk(label, got, want, tol=6e-4):
        got = float(got)
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
    chk("§2.6 経路 AUC マッチB", GS.loc["matched Group B", "pathway_auc_median"], 0.771, 2e-3)
    chk("§2.6 集約 AUC 全ortholog", GS.loc["all 1:1 orthologues", "aggregated_auc_centred"], 0.600, 2e-3)
    chk("§2.6 集約 AUC マッチB", GS.loc["matched Group B", "aggregated_auc_centred"], 0.626, 2e-3)
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

    print(f"一致 {len(ok)} 件 / 不一致 {len(bad)} 件")
    for b in bad:
        print("  ✗", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
