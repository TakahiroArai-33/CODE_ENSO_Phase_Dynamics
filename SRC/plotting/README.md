# 論文用の図

Arai・Kawamura の原稿と実際の EPS 図（2026-10-30時点）で使われている44パネルだけを出力します。保存先は `SRC/figure/`、形式は PNG・EPS（計88ファイル）です。ファイル名とフォルダ名は原稿の図番号に揃えています。

推定から作図まで一括実行：

```bash
bash SRC/enso_ERSSTv5_2026/run_all.sh
```

推定済みデータから作図だけを実行：

```bash
bash SRC/enso_ERSSTv5_2026/run_all.sh --figures-only
```

環境構築・NOAA データ取得は [ルートREADME](../../README.md) を参照してください。MATLAB の計算結果5ファイルを使い、Python が推定・作図します。

| 出力フォルダ | 原稿で使うパネル | パネル数 |
| --- | --- | ---: |
| `fig1` | Annual・ENSO の複素平面と各8枚のSST位相マップ | 18 |
| `fig2` | Qe、ODE軌道、El Niño・La Niña のSDE通過時間 | 4 |
| `fig3` | Q 平均化の位相差方程式（3:1） | 1 |
| `fig4` | 直接推定の位相差方程式とヒストグラム（3:1） | 2 |
| `fig5` | Niño3.4 と再構成時系列、2期間を緑色表示 | 1 |
| `fig6` | El Niño・La Niña の4月・9月開始軌道、カラーバー付き | 4 |
| `fig_s1` | 振幅フィルター10〜50%の Qe | 5 |
| `fig_s2` | 3:1 の直接推定の振幅別重ね描き | 1 |
| `fig_s3` | 2:1・4:1・5:1・5:2 の位相差方程式とヒストグラム | 8 |

比率は Annual:ENSO の順です。位相差は原稿と同じ ψ/(2π) 表示を用います。NPZ の単位は [単位の説明](../enso_ERSSTv5_2026/SRC_py/README.md#npz-units) に従います。
原稿のEPSは複数パネルを配置した完成図です。ここではその構成パネルを出力します。

全パネルが揃ってから保存先を置き換えるので、途中で失敗した場合は前回の完成済み出力を保持します。PNG と EPS は直接保存し、Ghostscript は使いません。

通過時間の SDE は各開始月5000試行です。動作確認には次を使えます。

```bash
SDE_TRIAL_COUNT=5 bash SRC/enso_ERSSTv5_2026/run_all.sh --figures-only --output-dir /tmp/ersstv5-figure-check
```

少数試行の結果は確認用です。実際の試行数と入力ハッシュは `run_metadata.json`、パネルと生成元の対応は `figure_manifest.tsv` に保存します。
