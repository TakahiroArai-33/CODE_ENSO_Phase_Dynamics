#!/usr/bin/env python3
# coding: utf-8
"""Fig3の位相差方程式。横軸はpsi/(2pi)、縦軸はcycle/year。"""

from plot_config import FIGURE_OUTPUT_DIR

import Plot_Averaging as original
from PhasePlot_PsiCycles import MATERIALS_DIR, plot_phase_equations


def main(DATAKEY=original.DATAKEY, AMP_val=original.AMP_val,
         AVG_M_val=original.AVG_M_val, RATIO_val=original.RATIO_val,
         savedir=None, png_output=True, eps_output=True, add_shade=False):
    amp = int(round(AMP_val * 100.0 if AMP_val <= 1.0 else AMP_val))
    path = original.averaging_result_path(original.cached_data(DATAKEY), int(AVG_M_val), amp, RATIO_val)
    result = original.load_result(path)
    ratio = original.result_ratio_pair(result, RATIO_val)
    if savedir is None:
        folder = f"ratio{original.ratio_filename(RATIO_val)}"
        savedir = MATERIALS_DIR / "Averaging" / folder / f"Averaging_amp{amp}_M{AVG_M_val}_{folder}_{DATAKEY}"
    return plot_phase_equations(original, result, ratio, savedir, add_shade, png_output, eps_output)


if __name__ == "__main__":
    main()
