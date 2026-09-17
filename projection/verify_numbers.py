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
            ("§2.3", "16 of 422 reaching 0.90 and 25 falling to 0.60 or below"),
            ("Fig3 キャプション", "244 of 422 pathways fall below 50 genes"),
            ("Fig3 キャプション", "across the 74 pathways of the three original collections"),
            ("Fig3 キャプション", "adding the 104 Reactome sets"),
            ("Fig3 キャプション", "would widen it to 178 columns"),
            ("Fig4 キャプション", "within each of the 422 pathways"),
            ("§3.3", "the interquartile range across the 422 pathways is 0.707 to 0.808"),
            ("§3.4", "Kendall's W across 422 pathways is 0.320"),
            ("§3.4", "244 of our 422 pathways fall below 50 genes"),
            ("§4.13", "held to the 74 sets of the three original collections"),
            ("§4.13", "rather than to all 178"),
        ]:
            (ok if phrase in body else bad).append(f"{where} の文言が見つからない: {phrase}")

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
            ("§2.3", "the aggregated values there being 0.600, 0.624 and 0.529"),
            ("§2.6", "in the other three spaces it falls to 0.600, 0.624 and 0.529")]),
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
            ("§2.2", "a median of 0.448 and a maximum of 0.548"),
            ("§2.6", "0.548 vs 0.765; 0.572 vs 0.780; 0.438 vs 0.765; 0.534 vs 0.769")]),
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
                       (r"(?i)\bceiling\b", "ceiling（benchmark に統一済みのはず）")]:
        hit = _re.findall(pat, joined)
        (ok if not hit else bad).append(
            f"{label} は残っていない" if not hit else f"{label} が {len(hit)} 箇所残っている")

    abst = _re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", all_md["FRONTMATTER.md"], _re.S)
    nw = len(abst.group(1).split())
    (ok if nw <= 200 else bad).append(
        f"Abstract は {nw} 語（200 語以内）" if nw <= 200 else f"Abstract が {nw} 語ある")

    # 四条件の AUC は §2.3・§3.3・Abstract・Conclusions の 4 箇所に出る。全部を検査する。
    for where, phrase in [
        ("§2.3", "separating at AUC 0.868"),
        ("§3.3", "pathway means separate at 0.868"),
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
    # 背表紙の一覧は、キャプションの集合に File S1 を足したものと一致すべき
    want_back = supp_caps | {"File S1"}
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

    print(f"一致 {len(ok)} 件 / 不一致 {len(bad)} 件")
    for b in bad:
        print("  ✗", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
