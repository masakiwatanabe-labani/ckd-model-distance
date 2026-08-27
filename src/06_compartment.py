"""皮質 vs 髄質の対応ありデザイン解析。

探索段階では 18 頭中 17 頭が皮質と共通の個体だった。個体差を除いた比較ができるため、
「進行軸が直交する」という主張は対応ありで検定し直すこと。

出力:
  - 対応のある個体だけを使った皮質／髄質の進行コントラスト
  - 進行軸どうしの Spearman（+ ブートストラップ CI）
  - 区画交互作用（同一個体内で、病期変化が区画で異なる遺伝子）
  - 尿濃縮軸と MR 下流軸を区画別に追跡
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from lib_stats import (load_config, get_logger, moderated_ttest, bh,  # noqa: E402
                       leave_one_out_contrast, append_summary)

log = get_logger("06_compartment")
cfg = load_config()
ROOT = Path(cfg["_root"])
INT = ROOT / cfg["paths"]["interim"]
RES = ROOT / cfg["paths"]["results"]
SEED = cfg["seed"]


def load_modules():
    with open(ROOT / "config" / "modules.yaml", encoding="utf-8") as fh:
        return {k: str(v).split() for k, v in yaml.safe_load(fh).items()}


def main():
    ctx = pd.read_parquet(INT / "cat_rna_ctx.parquet")
    med = pd.read_parquet(INT / "cat_rna_med.parquet")
    gctx = pd.read_parquet(INT / "cat_rna_ctx_grp.parquet").squeeze()
    gmed = pd.read_parquet(INT / "cat_rna_med_grp.parquet").squeeze()

    shared = sorted(set(ctx.columns) & set(med.columns))
    log.info("対応のある個体: %d 頭", len(shared))
    genes = ctx.index.intersection(med.index)
    ctx_p, med_p = ctx.loc[genes, shared], med.loc[genes, shared]
    grp = gctx[shared]
    payload = {"対応のある個体数": len(shared), "共通遺伝子数": len(genes),
               "群構成": grp.value_counts().to_dict()}

    # ---------- 対応あり個体だけでの進行コントラスト ----------
    a = grp[grp == "CKD3/4"].index
    b = grp[grp == "CKD1/2"].index
    prog_ctx = moderated_ttest(ctx_p[a], ctx_p[b])
    prog_med = moderated_ttest(med_p[a], med_p[b])
    prog_ctx.to_csv(RES / "paired_progression_cortex.csv")
    prog_med.to_csv(RES / "paired_progression_medulla.csv")

    idx = prog_ctx.index.intersection(prog_med.index)
    x, y = prog_ctx.loc[idx, "lfc"], prog_med.loc[idx, "lfc"]
    ok = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    rho, p = stats.spearmanr(x[ok], y[ok])
    # ブートストラップ CI（遺伝子リサンプリング）
    rng = np.random.default_rng(SEED)
    xv, yv = x[ok].to_numpy(), y[ok].to_numpy()
    boot = np.array([stats.spearmanr(*(lambda i: (xv[i], yv[i]))(rng.integers(0, len(xv), len(xv))))[0]
                     for _ in range(1000)])
    payload["進行軸 皮質×髄質 rho"] = f"{rho:.3f} (95%CI {np.percentile(boot,2.5):.3f}–{np.percentile(boot,97.5):.3f}, n={ok.sum()})"
    log.info("進行軸 rho=%.3f", rho)

    # 参考：各病期での皮質×髄質
    for g, label in [("CKD1/2", "早期"), ("CKD3/4", "晩期")]:
        cols = grp[grp == g].index
        ctrl = grp[grp == "Control"].index
        dc = ctx_p[cols].mean(axis=1) - ctx_p[ctrl].mean(axis=1)
        dm = med_p[cols].mean(axis=1) - med_p[ctrl].mean(axis=1)
        o = dc.notna() & dm.notna()
        payload[f"{label} 皮質×髄質 rho"] = round(float(stats.spearmanr(dc[o], dm[o])[0]), 3)

    # ---------- 区画交互作用（対応あり差分の差分） ----------
    ctrl = grp[grp == "Control"].index
    d_ctx = ctx_p.sub(ctx_p[ctrl].mean(axis=1), axis=0)
    d_med = med_p.sub(med_p[ctrl].mean(axis=1), axis=0)
    diff = d_ctx - d_med          # 同一個体内の区画差
    dis = grp[grp != "Control"].index
    inter = moderated_ttest(diff[dis], diff[ctrl])
    inter["q"] = bh(inter["p"].to_numpy())
    inter.sort_values("p").to_csv(RES / "compartment_interaction.csv")
    payload["区画交互作用 q<0.05"] = int((inter.q < 0.05).sum())

    # ---------- 主要軸の区画別追跡 ----------
    mods = load_modules()
    focus = ["Osmotic_TonEBP", "Urea_concentration", "TAL", "CD_principal",
             "CD_intercalated", "MR_downstream", "MR_regulators",
             "Lymphocyte", "B_plasma", "T_cell", "ECM_fibrosis", "Osteoclast_like_Mac"]
    rows = []
    for name in focus:
        gs = [g.upper() for g in mods.get(name, [])]
        for tis, mat in [("cortex", ctx_p), ("medulla", med_p)]:
            m = mat.copy()
            m.index = [str(i).upper() for i in m.index]
            m = m[~m.index.duplicated()]
            g = [x for x in gs if x in m.index]
            sub = m.loc[g]
            sub = sub[sub.std(axis=1) > 0]
            if len(sub) < 3:
                continue
            ev = ((sub.T - sub.mean(axis=1)) / sub.std(axis=1)).mean(axis=1)
            rec = {"module": name, "tissue": tis, "n_genes": len(sub)}
            for gg in ["Control", "CKD1/2", "CKD3/4"]:
                rec[gg] = round(float(ev[grp[grp == gg].index].mean()), 3)
            aa = ev[a].to_numpy(); bb = ev[b].to_numpy()
            obs = aa.mean() - bb.mean()
            pool = np.concatenate([aa, bb])
            rr = np.random.default_rng(SEED)
            null = np.array([(lambda q: q[:len(aa)].mean() - q[len(aa):].mean())(rr.permutation(pool))
                             for _ in range(cfg["permutation"]["n_perm_module"])])
            rec["prog_delta"] = round(float(obs), 3)
            rec["prog_p_perm"] = round(float((np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)), 4)
            loo = [np.delete(aa, k).mean() - bb.mean() for k in range(len(aa))] + \
                  [aa.mean() - np.delete(bb, k).mean() for k in range(len(bb))]
            rec["loo_min"], rec["loo_max"] = round(min(loo), 3), round(max(loo), 3)
            rec["sign_stable"] = bool(np.sign(min(loo)) == np.sign(max(loo)))
            rows.append(rec)
    fdf = pd.DataFrame(rows)
    fdf.to_csv(RES / "focus_modules_by_compartment.csv", index=False)
    payload["主要軸（進行 p<0.05 かつ LOO 安定）"] = (
        fdf.query("prog_p_perm<0.05 and sign_stable")[["module", "tissue", "prog_delta", "prog_p_perm"]]
        .to_dict("records"))

    append_summary("06_compartment", payload, cfg)


if __name__ == "__main__":
    main()
