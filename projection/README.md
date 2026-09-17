# 猫CKD参照軸に対するマウスモデルの射影解析

猫自然発症CKD（GSE303653 / Li et al., *Commun Biol* 2025）を固定参照軸として、
マウスCKDモデル（Pod-TRECK、IRI/GSE98622）がその軸をどれだけ再現するかを、
**方向（cos θ）と応答振幅（ν）に分けて**定量する。

このディレクトリは自己完結している。上位の `ckd_xspecies` からは
`src/lib_stats.py`（moderated t）、`config/config.yaml`、`data/interim/*.parquet`、
`data/ref/*`（オルソログ表・遺伝子セット）を読む。

---

## 構成

```
manuscript/        本文（IJMS 形式・英語）
  INTRODUCTION.md    §1
  RESULTS.md         §2.1–2.6 ＋ 図表キャプション
  DISCUSSION.md      §3.1–3.9
  METHODS.md         §4（数式つき）
  NUMBERS.md         ★数値の正本。本文の数値はここからのみ引く
  figures/           Fig 1–6（PNG 300dpi ＋ PDF、フォント埋め込み TrueType）
                       Fig 1 §2.1 幾何模式 + α/cos/ν の分解
                       Fig 2 §2.3 時間差 vs 方向/振幅
                       Fig 3 §2.4 共有対照バイアス + 観測 cos と天井
                       Fig 4 §2.5 経路サイズのノイズ床 + 状態x経路ヒートマップ
                       Fig 5 §2.5 経路間の順位 + 経路単位 vs 集約 + サイズ非依存
                       Fig 6 §2.6 3つの交絡候補

notes/             判断の記録（本文には出さない）
  FINDINGS.md              全解析の判定
  LITERATURE.md            先行研究調査（P）と差分
  MANTEL_S1.md             時間軸（L）
  PATHWAY.md               経路粒度（O）
  PROTEIN.md               蛋白層（Q、Supplementary）
  HUMAN_EXPLORATORY.md     ヒト外挿（K）— 不採用の経緯
  KPMP_TRANSCRIPTOMICS_SURVEY.md  KPMP のデータ調査（M）

results/           全出力（TSV ＋ 生成レポート）
*.py               解析スクリプト
*.tsv, *.txt       中間データ（Δ行列、遺伝子リスト、共変量）
```

## 主要な結果

| # | 主張 | 根拠 |
|---|---|---|
| §2.1 | **α は類似度として解釈できない** | 猫髄質CKD1/2 の α = 1.020 が基準自身の 1.000 を超える。同一データセット・同一コホート・同一プラットフォーム内 |
| §2.2 | α = ν·cos θ、R⊥ = sin θ。**α は必ず cos と ν に分けて報告する** | 恒等式は実データ16状態で assert 検証 |
| §2.3 | **軌道は回る。基準軸への角度と時間の関連は検出されない（Group A 空間）** | ペア間 Mantel は時間と +0.596（方向 +0.599 / 振幅 +0.312）。一方、固定軸への角度は時間と +0.074。全 orthologue 空間では +0.627 で再現しない（§2.7） |
| §2.4 | **異種の cos は 0.548 を超えない。測定限界ではない** | 最も保守的な天井でも最小 0.620、最接近でも余裕 0.171、48組中 0 組が超過 |
| §2.5 | **単一の最良モデルは存在しない／集約は種差を捨てる** | Kendall W = 0.362。§2.4 と同一設計（全120ペアの cos、対照共有を除外）で、全パネル 0.813 に対し経路単位 中央 0.788、集約後は 0.511（偶然水準） |
| §2.6 | 写像品質・Δ精度・平均発現量のいずれも種差を説明しない | 3交絡を個別に検査 |

## 再現手順

```bash
cd ckd_xspecies/projection
V=../.venv/bin/python

$V build_delta_matrix.py        # Δ行列（16状態）＋ 進行軸補助行列
$V build_precision.py           # Δ精度の共変量
$V alpha_decomp.py --demo --out /tmp/demo      # 自己検証（恒等式と方向/振幅の判別）
$V alpha_decomp.py --delta delta_matrix.tsv --ref cat_CKD34 \
     --genes groupA_intersection.txt --time time_map.tsv --out results/alpha_groupA
$V cos_matrix.py                # 対照解析（基準軸を固定しない cos 行列）
$V reliability.py               # 減衰の天井（約3分）
$V ceiling_validation.py        # 天井推定の妥当性検査
$V auc_permutation.py           # 種差の状態レベル並べ替え
$V mantel_s1.py                 # 時間軸（Mantel の方向/振幅分解）
$V pathway_cos.py && $V pathway_analysis.py    # 経路粒度
$V protein_decomp.py && $V protein_ceiling_check.py   # 蛋白層（Supplementary）
$V check1b_recipe.py            # Δ定義の感度分析（8変種）
$V pathway_separation.py        # 経路単位の種分離（§2.4 と同一設計）
$V make_figures.py              # Fig 1–6
$V ceiling_validation.py        # 天井の妥当性検査（pairs_ceilings / ceiling_comparison）
./run_aa.sh && ./run_aa4.sh     # AA. 遺伝子空間の感度解析（3空間）
$V aa_summary.py                # 同 集計表 → results/aa_genespace/
$V code/revision1/pair_uncertainty_both_ceilings.py   # 改訂1: 天井2種の個体ブートストラップ
$V code/revision1/gap_by_group.py                    # 改訂1: gap の3群集計と状態ラベル置換
$V code/revision1/cosine_attenuation_sim.py          # 改訂1: 減衰式の当てはまり（Figure S1）
cd code/revision1 && $V pathway_cos.py && $V pathway_analysis.py && $V pathway_separation.py
$V code/revision1/genespace_with_reactome.py         # 改訂1: 遺伝子空間 x Reactome
$V code/revision1/feline_filter_check.py             # 改訂1: 猫側フィルタの感度
$V verify_numbers.py            # ★NUMBERS.md と結果ファイルの照合（152件）
```

Group A/B のマッチング（共変量を差し替えて3回）:

```bash
for e in control_log2cpm delta_se_moderated delta_se_splithalf; do
  $V check1_geneset.py --delta delta_matrix.tsv --x cat_CKD34 --y mouse_2W \
    --core groupA_intersection.txt --core-label GroupA --added-label GroupB \
    --expr $e.tsv --out results/groupAB_$e
done
```

## 守っている規則

- **Δ の定義は S1 に固定**: Ensembl one2one オルソログ / 群平均 FPKM ≥ 1 / 素の log2FC。
  8変種の感度分析は `results/check1/supplementary_sensitivity.tsv`。
- **本文の数値は `manuscript/NUMBERS.md` からのみ引く。** 執筆後に `verify_numbers.py` を通す。
  この照合は主要な主張の前提（例: 猫髄質CKD1/2 の α > 基準、異種ペアの天井非超過、
  caliper 2 での符号反転）を assert で固定しているので、壊れたら落ちる。
- **α 単独の表を作らない。** 必ず cos と ν を併記する。
- **観測 cos と減衰補正後 cos を混同しない。** 本文は原則すべて観測 cos。
- 個別経路の検定はしない。n ≥ 50 かつランダム帯 97.5% 超のセルのみ言及する。
- 状態数が一桁〜十数個の相関に「有意」と書かない。CI を示す。

## 不採用にしたもの（経緯は notes/ に残す）

- **ヒト外挿（§2.5 候補）**: cat × human 0.465 > mouse × human 0.390（AUC 0.898, p=0.0016）が
  出たが、写像品質を遺伝子ごとに揃える（|Δperc_id| ≤ 2）と **符号が反転**（−0.058, AUC 0.331）。
  片側のみが写像を経る非対称な比較のため成立しない。Limitations に記載。
- **状態ごとの時間相関**: n=12 で CI が [−0.7, +0.8]。帰無の支持と検出力不足を区別できないため
  結論に用いない。
