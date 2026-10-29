# 論文図に必要な位相方程式の推定

`../myFolder/` の MATLAB 出力から、二次元 Q・Q Averaging・Direct Averaging を計算します。
annual は g2、ENSO は g9。`common_fixed.py` で MATLAB datenum を変換し、annual を1月、ENSO を Niño3.4 ピークに合わせます。

論文図まで一括で作成する入口は、ひとつ上の `run_all.sh` です。環境と NOAA SST の取得後、プロジェクトのルートから実行します。

```bash
bash SRC/enso_ERSSTv5_2026/run_all.sh
```

下表の中間結果を計算し、原稿で使用する図だけを `SRC/figure/` に PNG・EPS で保存します。
既存の NPZ で作図だけをやり直す場合は `--figures-only` を付けます。
`--output-dir /tmp/ersstv5-figures` で保存先を変更できます。相対パスは実行時の作業ディレクトリを基準にします。
入力や依存パッケージが不足している場合は推定の前に停止します。NOAA の元データは自動取得しません。

| 保存先 | 推定条件 |
| --- | --- |
| `Interim_data_Q/` | M=3、振幅フィルター10・20・30・40・50%（5ファイル） |
| `Interim_Averaging/` | M=3、30%、ratio=1:3（1ファイル） |
| `Interim_data_DirectAveraging/` | M=3、30%で全5比率、1:3のみ10・20・40・50%も追加（9ファイル） |

`ratio` は ENSO:Annual の順です。本文図の Annual:ENSO=3:1 は `ratio=1:3` に対応します。
既定の Conda 環境は `virtual_env`。`CONDA_ENV="" bash .../run_all.sh` では現在の `python` を使用します。
MATLAB 計算データは変更しません。環境・NOAA 入力の準備は [ルートREADME](../../../README.md) を参照してください。

`SRC_py/run_all_fixed.sh` は中間結果だけを作る内部用の入口です。必要な場合に限り、ルートから次を実行します。

```bash
bash SRC/enso_ERSSTv5_2026/SRC_py/run_all_fixed.sh
```

`--check-inputs` を付けると、推定せずに MATLAB 入力と推定用の依存関係だけを確認します。

単位は [解析結果NPZの単位](#npz-units) を確認してください。Q は rad/year、平均化結果の `y1`〜`y3` は cycles/year です。

## NPZ units

Time is measured in years. Do not infer units from key names alone.

| Stored quantity | Units |
| --- | --- |
| Phases and phase differences: `phi*`, `x` | rad |
| Q estimates: `Q`, `omega`; Direct Averaging: `omega1`, `omega2` | rad/year |
| Q Averaging and Direct Averaging: `y1`–`y3`, `y*_err` | cycles/year |
| Q noise covariance: `noise_cov` (E = 2D) | rad²/year |
| Direct Averaging: `D1`–`D3` (√(2D)) | rad/√year |
| Q Averaging: `Eaa`, `Eee`, `Eae`, `E_delta` | cycles²/year |

Divide radian-based drift and angular frequency by 2π, and noise covariance by (2π)², to obtain cycle-based units. Averaged `y*` values are already in cycles/year and must not be converted again.
