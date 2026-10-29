#!/usr/bin/env python3
# coding: utf-8
"""FigS2の位相差方程式。各ampの推定値をcycles表示で重ね描きする。"""

from plot_config import FIGURE_OUTPUT_DIR

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import Plot_DirectAveraging as base
import Plot_DirectAveraging_AmpOverlay as original
from PhasePlot_PsiCycles import (
    LABEL_FONTSIZE, MATERIALS_DIR, PANEL_FIGSIZE, TWOPI, phase_title, phase_xlabel, phase_rate_ylabel,
    psi_symbol, result_in_psi, setup_phase_axis, variant_stem,
    load_averaging_ylim_reference, set_phase_ylim,
)


def main(DATAKEY=original.DATAKEY, AMP_VALUES=original.AMP_VALUES,
         DA_M_val=original.DA_M_val, RATIO_val=original.RATIO_val,
         savedir=None, png_output=True, eps_output=True,
         ylim_reference_amp=original.YLIM_REFERENCE_AMP, panels=("Psi",)):
    results = [(amp, result_in_psi(result)) for amp, result in original.load_amp_results(
        DATAKEY=DATAKEY, AMP_VALUES=AMP_VALUES, DA_M_val=DA_M_val, RATIO_val=RATIO_val,
    )]
    _, reference = original.select_ylim_reference_result(results, ylim_reference_amp)
    ratio = base.result_ratio_pair(reference, RATIO_val)
    ratio_key = base.ratio_key_from_pair(ratio)
    ylim_base, ylim_result = base, reference
    if ratio == (1, 3):
        ylim_base, ylim_result = load_averaging_ylim_reference(DATAKEY, DA_M_val, 30, "1:3")
    if savedir is None:
        folder = f"ratio{base.ratio_filename(RATIO_val)}"
        savedir = MATERIALS_DIR / "DirectAveraging" / folder / f"DirectAveraging_overlay_M{DA_M_val}_{folder}_{DATAKEY}"
    saved = []
    for panel, key, phase in (
        ("Annual", "y1", r"\theta_a"),
        ("ENSO", "y2", r"\theta_e"),
        ("Psi", "y3", psi_symbol(ratio)),
    ):
        if panel not in panels:
            continue
        fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
        for index, (amp, result) in enumerate(results):
            ax.plot(np.asarray(result["x"]) / TWOPI, result[key],
                    color=original.overlay_color(index, amp),
                    linewidth=original.overlay_linewidth(amp), label=f"{amp}%")
        setup_phase_axis(ax, phase_xlabel(ratio))
        ax.set_title(phase_title(ratio), fontsize=LABEL_FONTSIZE)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.set_ylabel(phase_rate_ylabel(phase), fontsize=LABEL_FONTSIZE)
        set_phase_ylim(ax, ylim_base, ylim_result, key, ratio_key, original.YLIM_REFERENCE_ADD_SHADE)
        fig.tight_layout()
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.30),
                  ncols=original.overlay_legend_ncols(len(results)), fontsize=6,
                  frameon=False, borderaxespad=0.0, handlelength=1.3,
                  handletextpad=0.35, columnspacing=0.8)
        saved.extend(base.save_figure(fig, variant_stem(savedir), panel, png_output, eps_output))
    return saved


if __name__ == "__main__":
    main()
