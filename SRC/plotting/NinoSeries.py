#!/usr/bin/env python3
# coding: utf-8

# %% 



"""
nino3.4時系列と固有関数のプロット。
ただし、左固有関数を用いている。
(froyland et al. 2021, Nat. Commun.からFig. 10, Mode Reconstructionを参照。)
"""

from plot_config import FIGURE_OUTPUT_DIR

import datetime
from pathlib import Path
import os
import numpy as np


import matplotlib.pyplot as plt
from Common import (
    get_phase_series,
    get_data,
    get_SST_phase_composite,
)

plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams['font.family'] = 'Times New Roman' # font familyの設定
plt.rcParams['font.size'] = 10 # 全体のフォントサイズが変更されます。
plt.rcParams['xtick.labelsize'] = 10 # 軸だけ変更されます.
plt.rcParams['ytick.labelsize'] = 10 # 軸だけ変更されます.
plt.rcParams['figure.facecolor']= 'white'


DATAKEY = "ERSSTv5_2026"
S_val = int(8)
AMP_val = float(0.3)
EMB_WINDOW = int(48) # froyland論文におけるQのこと、Takens delay embedding window

# 網掛けする期間。終了日は含まない（1986年9月〜1988年1月）。
DOUBLE_PEAK_PERIODS = (
    (datetime.datetime(1986, 9, 1), datetime.datetime(1988, 2, 1)),
)
# ラニーニャ側の2つの谷（1998年7月〜2000年7月）。
ADDITIONAL_DOUBLE_PEAK_PERIODS = (
    (datetime.datetime(1998, 7, 1), datetime.datetime(2000, 8, 1)),
)


def mode_reconstruction(target, right_modes, left_modes, mu, emb_window):
    """
    ensoの固有関数を用いてnino3.4を再構成する。
    (froyland et al. 2021, Nat. Commun.からFig. 10, Mode Reconstructionを参照。)
    """
    y = np.asarray(target)
    squeeze_output = y.ndim == 1
    if squeeze_output:
        y = y[:, None]

    g = np.asarray(right_modes)
    g_left = np.asarray(left_modes)
    mu = np.asarray(mu)

    weighted_left = np.conjugate(g_left) * mu[:, None]
    gram_inv = np.linalg.inv(weighted_left.T @ g)

    n = y.shape[0]
    lags = np.arange(emb_window) - emb_window // 2
    reconstruction = np.zeros((n, y.shape[1]), dtype=np.result_type(g, y, complex))
    counts = np.zeros(n, dtype=float)

    for lag in lags:
        if lag >= 0:
            dual_slice = slice(lag, n)
            target_slice = slice(0, n - lag)
            mode_slice = slice(lag, n)
            output_slice = slice(0, n - lag)
        else:
            dual_slice = slice(0, n + lag)
            target_slice = slice(-lag, n)
            mode_slice = slice(0, n + lag)
            output_slice = slice(-lag, n)

        coeff = gram_inv @ (weighted_left[dual_slice].T @ y[target_slice])
        reconstruction[output_slice] += g[mode_slice] @ coeff
        counts[output_slice] += 1.0

    reconstruction /= counts[:, None]
    reconstruction = np.real_if_close(reconstruction)
    return reconstruction[:, 0] if squeeze_output else reconstruction



#%%
def main():

    output_dir = FIGURE_OUTPUT_DIR

    #* datas key: "annual_enso", "all"
    #* datas["annual_enso"] key: "date", "nino3.4", "Z_annual", "Z_enso",
    #*     "Z_annual_conj", "Z_enso_conj", "raw_Z_annual", "raw_Z_enso",
    #*     "raw_Z_annual_conj", "raw_Z_enso_conj", "raw_Z_annual_left",
    #*     "raw_Z_enso_left", "raw_Z_annual_conj_left", "raw_Z_enso_conj_left",
    #*     "lambda_annual", "lambda_enso", "lambda_annual_conj",
    #*     "lambda_enso_conj", "gamma_annual", "gamma_enso",
    #*     "gamma_annual_conj", "gamma_enso_conj", "mode_info", "phase_origin"
    #* datas["all"] key: "key", "dataset_dir", "myfolder", "date", "Z",
    #*     "Z_left", "nino34", "Lambda", "mu", "mode_info", "phase_origin"
    datas = get_data(key=DATAKEY)

    #* phase_datas key: "date", "phase_annual", "phase_enso",
    #*     "phase_annual_conj", "phase_enso_conj", "amp_annual", "amp_enso",
    #*     "amp_annual_conj", "amp_enso_conj"
    phase_datas = get_phase_series(datas)

    #* sst_phase_composite key: "lon", "lat", "composite", "counts",
    #*     "phase_index", "key1", "key2", "S", "amp", "amp_percentile"


    #* -----
    annual_enso = datas["annual_enso"]
    all_data = datas["all"]

    date = annual_enso["date"]
    nino34 = np.real(np.asarray(annual_enso["nino3.4"]))
    raw_Z_annual = np.asarray(annual_enso["raw_Z_annual"])
    raw_Z_enso = np.asarray(annual_enso["raw_Z_enso"])
    raw_Z_enso_conj = np.asarray(annual_enso["raw_Z_enso_conj"])
    raw_Z_enso_left = np.asarray(annual_enso["raw_Z_enso_left"])
    raw_Z_enso_conj_left = np.asarray(annual_enso["raw_Z_enso_conj_left"])


    # ENSO mode とその共役を並べ、複素共役 pair 全体で Nino3.4 を再構成する。
    # (right eigenfunctiondであり、これを回転したZ_ensoから位相を計算している。)
    enso_modes = np.column_stack((
        raw_Z_enso,
        raw_Z_enso_conj,
    ))

    # 左固有関数も同じ mode pair の順序で並べ、双対基底として射影に使う。
    enso_left_modes = np.column_stack((
        raw_Z_enso_left,
        raw_Z_enso_conj_left,
    ))

    # 右/左固有関数、測度 mu、遅延窓を使って Fig.10 型の mode reconstruction を計算する。
    nino34_enso_reconstruction = mode_reconstruction(
        nino34,
        enso_modes,
        enso_left_modes,
        all_data["mu"],
        EMB_WINDOW,
    )

    # 数値誤差で残る虚部を落とし、以降は実数時系列として扱う。
    nino34_enso_reconstruction = np.real(np.asarray(nino34_enso_reconstruction))
    # Nino3.4 と再構成値の両方が有限な時刻だけをプロット対象にする。
    valid_plot = np.isfinite(nino34) & np.isfinite(nino34_enso_reconstruction)
    if not np.any(valid_plot):
        raise ValueError("No finite samples found for Nino3.4 reconstruction plot.")



    #* --- プロット2 ----
    x = date[valid_plot]
    y1 = nino34[valid_plot]
    y2 = nino34_enso_reconstruction[valid_plot]
    y1_color = "k"
    y2_color = "tomato"
    lw = 1.0
    xmin = date[0]
    xmax = date[-1]


    fig, ax = plt.subplots(figsize=(4, 2.5))
    ax.plot(x, y1, color=y1_color, lw=lw, label="Nino3.4", zorder=5)
    ax.plot(x, y2, color=y2_color, lw=lw, label=r"$Y(\omega_\ell)$", zorder=4)

    ax.set_xlabel("year")
    ax.set_ylabel("[°C]")
    ax.grid(linestyle="-.", alpha=0.5)
    ax.set_xlim(xmin, xmax)
    xtick_years = [1980, 1990, 2000, 2010, 2020]
    ax.set_xticks([datetime.datetime(year, 1, 1) for year in xtick_years])
    ax.set_xticklabels([str(year) for year in xtick_years])
    ax.set_yticks([-3, -2, -1, 0, 1, 2, 3]) 

    ax.legend(
        fontsize=8,
        loc="lower left",
        ncol=2
    )

    # Fig5は両期間を重ねた版だけを使用する。
    for start, end in DOUBLE_PEAK_PERIODS + ADDITIONAL_DOUBLE_PEAK_PERIODS:
        ax.axvspan(
            start, end,
            # limegreenを白で薄め、PNG・EPSで同じ淡い色にする。
            facecolor="#c2f0c2", edgecolor="none",
            linewidth=0.0, zorder=0.5,
        )

    for extension in ("eps", "png"):
        output = output_dir / f"nino_reconstruction_double_peaks.{extension}"
        fig.savefig(output, dpi=300, format=extension, bbox_inches="tight")
        print(f"saved figure: {output}")

    plt.close(fig)


if __name__ == "__main__":
    main()
