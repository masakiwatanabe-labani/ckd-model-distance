"""Part 6. Group A の検出判定を、遺伝子ごとに追えるかたちで書き出す。

Group A は全解析の基盤なので、遺伝子リストだけでなく「どのプロテオームで、何サンプル中
何サンプルで非欠測だったか」を出す。判定基準は build_delta_matrix.group_a_sets と同じ
（そのデータセットのサンプルの 50% 以上で非欠測）で、ここではその中間量を保存する。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import build_delta_matrix as B  # noqa: E402

OUT = HERE / "results" / "round7"
OUT.mkdir(parents=True, exist_ok=True)
FRAC = 0.5


def main() -> int:
    rows = {}
    meta = []
    detected_sets = {}
    for name, species, label in [("cat_prot_ctx", "cat", "cat_cortex"),
                                 ("cat_prot_med", "cat", "cat_medulla"),
                                 ("mouse_prot", "mouse", "podtreck")]:
        d = B.read_int(name)
        n_samp = d.shape[1]
        nonmiss = d.notna().sum(axis=1)
        frac = d.notna().mean(axis=1)
        native = pd.DataFrame({"n_nonmissing": nonmiss, "frac_nonmissing": frac.round(4)})
        native["detected"] = native.frac_nonmissing >= FRAC
        # 判定は build_delta_matrix.group_a_sets と同一の手順で行う。検出された遺伝子だけを
        # ヒト空間へ写し、その index を「検出された」とする。
        det = pd.Series(1.0, index=native.index[native.detected])
        hd, _ = B.to_human(det, species)
        # 併せて、判定の根拠として非欠測割合もヒト空間へ写す（複数が写る場合は平均）。
        h, _ = B.to_human(native.frac_nonmissing, species)
        hn, _ = B.to_human(native.n_nonmissing.astype(float), species)
        rows[f"{label}_frac_nonmissing"] = h
        rows[f"{label}_n_nonmissing"] = hn
        rows[f"{label}_n_samples"] = pd.Series(n_samp, index=h.index)
        detected_sets[label] = set(hd.index)
        meta.append((label, name, d.shape[0], n_samp, int(native.detected.sum()), len(h)))
        print(f"{label:12s} native genes {d.shape[0]:6d}  samples {n_samp:3d}  "
              f"detected {int(native.detected.sum()):6d}  in human space {len(h):6d}")

    T = pd.DataFrame(rows)
    for label in ("cat_cortex", "cat_medulla", "podtreck"):
        T[f"{label}_detected"] = T.index.isin(detected_sets[label])
    T["groupA_intersection"] = (T.cat_cortex_detected & T.cat_medulla_detected
                                & T.podtreck_detected)
    T["groupA_union"] = ((T.cat_cortex_detected | T.cat_medulla_detected)
                         & T.podtreck_detected)
    T.index.name = "gene"

    D = pd.read_csv(HERE / "delta_matrix.tsv", sep="\t", index_col=0)
    # 公表リストは「検出の交差」をさらに Δ 行列に存在する遺伝子へ絞ったもの
    # （build_delta_matrix.py の inter_in）。同じ絞りをかけて照合する。
    inter = sorted(g for g in T.index[T.groupA_intersection] if g in D.index)
    published = {g.strip() for g in (HERE / "groupA_intersection.txt").read_text().split()
                 if g.strip()}
    print(f"\n再構成した Group A intersection {len(inter)} / 公表リスト {len(published)} "
          f"/ 一致 {len(set(inter) & published)}")
    assert set(inter) == published, (
        f"再構成した Group A が公表リストと一致しない: 余分 "
        f"{sorted(set(inter) - published)[:5]} 不足 {sorted(published - set(inter))[:5]}")

    T["in_delta_matrix"] = T.index.isin(D.index)
    T["complete_in_16_states"] = T.index.isin(D.dropna(axis=0, how="any").index)
    T.round(4).to_csv(OUT / "groupA_detection_calls.tsv", sep="\t")

    M = pd.DataFrame(meta, columns=["proteome", "source_table", "genes_in_native_space",
                                    "n_samples", "genes_detected_native",
                                    "genes_mapped_to_human"])
    M.to_csv(OUT / "groupA_detection_sources.tsv", sep="\t", index=False)
    n_used = int((T.groupA_intersection & T.complete_in_16_states).sum())
    print(f"Group A intersection のうち 16 状態すべてで有限値: {n_used}")
    print(M.to_string(index=False))
    print(f"\n書き出し: {OUT / 'groupA_detection_calls.tsv'}, "
          f"{OUT / 'groupA_detection_sources.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
