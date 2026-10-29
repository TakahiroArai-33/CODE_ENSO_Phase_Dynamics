#!/usr/bin/env python3
# coding: utf-8
"""Fig4・FigS3の位相差方程式と分布をcycles表示で出力する。"""

from plot_config import FIGURE_OUTPUT_DIR

from typing import get_args

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import Plot_DirectAveraging as original
from PhasePlot_PsiCycles import (
    LABEL_FONTSIZE, MATERIALS_DIR, PANEL_FIGSIZE, TWOPI, phase_title, phase_xlabel, plot_phase_equations,
    psi_symbol, setup_phase_axis, variant_stem, load_averaging_ylim_reference,
)


def main(DATAKEY=original.DATAKEY, AMP_val=original.AMP_val,
         DA_M_val=original.DA_M_val, RATIO_val=original.RATIO_val,
         savedir=None, png_output=True, eps_output=True, add_shade=True):
    amp = int(round(AMP_val * 100.0 if AMP_val <= 1.0 else AMP_val))
    path = original.direct_result_path(original.cached_data(DATAKEY), int(DA_M_val), amp, RATIO_val)
    result = original.load_result(path)
    ratio = original.result_ratio_pair(result, RATIO_val)
    if savedir is None:
        folder = f"ratio{original.ratio_filename(RATIO_val)}"
        savedir = MATERIALS_DIR / "DirectAveraging" / folder / f"DirectAveraging_amp{amp}_M{DA_M_val}_{folder}_{DATAKEY}"
    ylim_reference = None
    if ratio == (1, 3):
        ylim_reference = load_averaging_ylim_reference(DATAKEY, DA_M_val, 30, "1:3")
    saved = plot_phase_equations(original, result, ratio, savedir, add_shade, png_output, eps_output,
                                ylim_reference=ylim_reference)

    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    bins = np.asarray(result["distribution_bins"], dtype=float) / TWOPI
    # 横軸をradからcyclesへ変換し、密度も新しい横軸に対して正規化する。
    ax.hist(np.asarray(result["delta_phi_data"]) / TWOPI, bins=bins, density=True,
            color="0.75", edgecolor="k", linewidth=0.5)
    ax.hist(np.asarray(result["delta_phi_langevin"]) / TWOPI, bins=bins, density=True,
            histtype="step", color="blue", linewidth=1.0)
    setup_phase_axis(ax, phase_xlabel(ratio))
    ax.set_title(phase_title(ratio), fontsize=LABEL_FONTSIZE)
    ax.set_ylim(0.0, 2.5)
    ax.set_yticks([0.0, 1.25, 2.5])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.set_ylabel(rf"$P\left(\frac{{{psi_symbol(ratio)}}}{{2\pi}}\right)$", fontsize=LABEL_FONTSIZE)
    fig.tight_layout()
    saved.extend(original.save_figure(fig, variant_stem(savedir), "Hist", png_output, eps_output))
    return saved


if __name__ == "__main__":
    for ratio in get_args(original.RatioKey):
        main(RATIO_val=ratio)
