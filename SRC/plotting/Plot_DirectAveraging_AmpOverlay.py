
from plot_config import FIGURE_OUTPUT_DIR
#!/usr/bin/env python3
# coding: utf-8

# %%

from typing import get_args

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from Common import DatasetKey
from Plot_DirectAveraging import (
    BX_MIN_YSPAN,
    CX_MIN_YSPAN,
    PANEL_FIGSIZE,
    PSI_FIXED_YHALFSPAN,
    RatioKey,
    cached_data,
    delta_phi_expr_from_pair,
    direct_result_path,
    estimate_band_values,
    load_result,
    phase_rate_ylabel,
    ratio_filename,
    ratio_key_from_pair,
    result_ratio_pair,
    reverse_delta_phi_expr_from_pair,
    save_figure,
    set_centered_ylim,
    set_dynamic_ylim,
    set_three_yticks,
    setup_phase_axis,
)


DATAKEY = "ERSSTv5_2026"
#* 重ね描きしたい amp をここで任意に指定する。
AMP_VALUES: tuple[int, ...] = (10, 20, 30, 40, 50)
DA_M_val = int(3)
RATIO_val: RatioKey = "1:3"
YLIM_REFERENCE_AMP = int(30)
YLIM_REFERENCE_ADD_SHADE = True
OVERLAY_DEFAULT_LINEWIDTH = 0.8
OVERLAY_REFERENCE_LINEWIDTH = 1.0
OVERLAY_AMP_COLORS = {
    10: "#E69F00",
    20: "#009E73",
    30: "blue",
    40: "#CC79A7",
    50: "#D55E00",
}
OVERLAY_FALLBACK_COLORS = ("#F0E442", "#7F3C8D", "#11A579", "#A95AA1", "#3969AC")


def amp_rank(amp_val):
    return int(round(amp_val * 100.0 if amp_val <= 1.0 else amp_val))


def overlay_color(index, amp_rank_val):
    if amp_rank_val in OVERLAY_AMP_COLORS:
        return OVERLAY_AMP_COLORS[amp_rank_val]
    return OVERLAY_FALLBACK_COLORS[index % len(OVERLAY_FALLBACK_COLORS)]


def overlay_linewidth(amp_rank_val):
    if amp_rank_val == amp_rank(YLIM_REFERENCE_AMP):
        return OVERLAY_REFERENCE_LINEWIDTH
    return OVERLAY_DEFAULT_LINEWIDTH


def overlay_legend_ncols(total):
    return min(5, max(3, int(total)))


def select_ylim_reference_result(results, reference_amp=YLIM_REFERENCE_AMP):
    reference_amp_rank = amp_rank(reference_amp)
    for amp_rank_val, result in results:
        if amp_rank_val == reference_amp_rank:
            return amp_rank_val, result
    return results[len(results) // 2]


def load_amp_results(DATAKEY=DATAKEY, AMP_VALUES=AMP_VALUES, DA_M_val=DA_M_val, RATIO_val=RATIO_val):
    datas = cached_data(DATAKEY)
    results = []
    for amp_val in AMP_VALUES:
        amp_rank_val = amp_rank(amp_val)
        path = direct_result_path(datas, int(DA_M_val), amp_rank_val, RATIO_val)
        if not path.exists():
            print(f"skip missing: {path}")
            continue
        results.append((amp_rank_val, load_result(path)))
    if not results:
        raise FileNotFoundError(f"No DirectAveraging results found: {DATAKEY=}, {DA_M_val=}, {RATIO_val=}")
    return results


def plot_overlay_panel(results, y_key, xlabel, ylabel, savedir, panel_name, title=None, min_span=None,
                       fixed_ylim=None, fixed_yhalfspan=None, ylim_reference_result=None, err_key=None,
                       png_output=True, eps_output=False):
    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    y_values = []

    for i, (amp_rank, result) in enumerate(results):
        x = np.asarray(result["x"], dtype=float)
        y = np.asarray(result[y_key], dtype=float)
        y_values.append(y)
        ax.plot(
            x,
            y,
            color=overlay_color(i, amp_rank),
            linewidth=overlay_linewidth(amp_rank),
            label=f"{amp_rank}%",
        )

    if title:
        ax.set_title(title)
    setup_phase_axis(ax, xlabel)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    if fixed_ylim is not None:
        ax.set_ylim(*fixed_ylim)
    elif fixed_yhalfspan is not None:
        if ylim_reference_result is None:
            ylim_values = np.concatenate(y_values)
        else:
            ylim_values = estimate_band_values(
                ylim_reference_result,
                y_key,
                err_key,
                YLIM_REFERENCE_ADD_SHADE,
            )
        set_centered_ylim(ax, ylim_values, fixed_yhalfspan)
    elif min_span is not None:
        if ylim_reference_result is None:
            ylim_values = np.concatenate(y_values)
        else:
            ylim_values = estimate_band_values(
                ylim_reference_result,
                y_key,
                err_key,
                YLIM_REFERENCE_ADD_SHADE,
            )
        set_dynamic_ylim(ax, ylim_values, min_span)
    set_three_yticks(ax)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.20),
        ncols=overlay_legend_ncols(len(results)),
        fontsize=6,
        frameon=False,
        borderaxespad=0.0,
        handlelength=1.3,
        handletextpad=0.35,
        columnspacing=0.8,
    )
    return save_figure(fig, savedir, panel_name, png_output, eps_output)


def main(
        DATAKEY=DATAKEY,
        AMP_VALUES: tuple[int, ...] | list[int] = AMP_VALUES,
        DA_M_val=DA_M_val,
        RATIO_val: RatioKey = RATIO_val,
        savedir="./",
        png_output=True,
        eps_output=False,
        ylim_reference_amp=YLIM_REFERENCE_AMP,
):
    results = load_amp_results(
        DATAKEY=DATAKEY,
        AMP_VALUES=AMP_VALUES,
        DA_M_val=DA_M_val,
        RATIO_val=RATIO_val,
    )
    _, ylim_reference_result = select_ylim_reference_result(results, ylim_reference_amp)
    #* DirectAveraging の y1/y2/y3 は保存時に 2π で割ってあり、すでに yr^{-1} 単位。
    loaded_ratio = result_ratio_pair(ylim_reference_result, RATIO_val)
    ratio_key = ratio_key_from_pair(loaded_ratio)
    phase_expr = delta_phi_expr_from_pair(loaded_ratio)
    xlabel = "$" + phase_expr + "$"
    reverse_xlabel = "$" + reverse_delta_phi_expr_from_pair(loaded_ratio) + "$"

    saved = []
    saved.extend(plot_overlay_panel(
        results,
        y_key="y1",
        xlabel=xlabel,
        ylabel=r"$\frac{d}{dt}\theta_a$ [yr$^{-1}$]",
        savedir=savedir,
        panel_name="Annual",
        fixed_ylim=(0.9, 1.1),
        png_output=png_output,
        eps_output=eps_output,
    ))
    saved.extend(plot_overlay_panel(
        results,
        y_key="y2",
        xlabel=reverse_xlabel,
        ylabel=r"$\frac{d}{dt}\theta_e$ [yr$^{-1}$]",
        savedir=savedir,
        panel_name="ENSO",
        min_span=BX_MIN_YSPAN[ratio_key],
        ylim_reference_result=ylim_reference_result,
        err_key="y2_err",
        png_output=png_output,
        eps_output=eps_output,
    ))
    saved.extend(plot_overlay_panel(
        results,
        y_key="y3",
        xlabel=xlabel,
        ylabel=phase_rate_ylabel(phase_expr),
        savedir=savedir,
        panel_name="Psi",
        fixed_yhalfspan=PSI_FIXED_YHALFSPAN.get(ratio_key),
        min_span=CX_MIN_YSPAN[ratio_key],
        ylim_reference_result=ylim_reference_result,
        err_key="y3_err",
        png_output=png_output,
        eps_output=eps_output,
    ))
    return saved





if __name__ == "__main__":
    main(RATIO_val="1:3", savedir=f"DirectAveraging_overlay_M3_ratio1to3_{DATAKEY}",
         png_output=True, eps_output=True)
