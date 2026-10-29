#!/usr/bin/env python3
# coding: utf-8
"""変更後の位相方程式図用。psiの向きとcycles表示を揃える共通処理。"""

from plot_config import FIGURE_OUTPUT_DIR

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


TWOPI = 2.0 * np.pi
PANEL_FIGSIZE = (2.4, 1.7)
LABEL_FONTSIZE = 9
MATERIALS_DIR = FIGURE_OUTPUT_DIR / "素材"


def psi_symbol(ratio):
    # NPZのratioはENSO:Annual。論文の添字はAnnual:ENSO。
    n, m = ratio
    return rf"\psi_{{{m:d}\,{n:d}}}"


def phase_title(ratio):
    n, m = ratio
    coefficient = str(m) if n == 1 else rf"\frac{{{m:d}}}{{{n:d}}}"
    return rf"${psi_symbol(ratio)}:=\theta_a-{coefficient}\theta_e$"


def phase_xlabel(ratio):
    return rf"$\frac{{{psi_symbol(ratio)}}}{{2\pi}}$"


def phase_rate_ylabel(phase):
    return rf"$\frac{{d}}{{dt}}\left(\frac{{{phase}}}{{2\pi}}\right)$ [cycle / year]"


def setup_phase_axis(ax, xlabel):
    ax.set_xlim(0.0, 1.0)
    ax.set_xticks(np.linspace(0.0, 1.0, 5), ["0", "0.25", "0.5", "0.75", "1"])
    ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE)
    ax.grid(linestyle="--", alpha=0.6)


def result_in_psi(result):
    """y2と誤差帯を-psiからpsiへ反転する。速度の単位と符号は維持。"""
    x = np.asarray(result["x"], dtype=float)
    reflected_x = np.mod(-x, TWOPI)
    converted = dict(result)
    for key in ("y2", "y2_err"):
        if key in result:
            converted[key] = np.interp(reflected_x, x, np.asarray(result[key], dtype=float))
    return converted


def variant_stem(savedir):
    path = Path(savedir)
    if path.suffix in {".png", ".eps"}:
        path = path.with_suffix("")
    if not path.name.endswith("_psi_cycles"):
        path = path.with_name(path.name + "_psi_cycles")
    return path


def load_averaging_ylim_reference(datakey, m_val, amp, ratio):
    """直接推定版の比較に使うQ平均化の結果（速度は既にcycle/year）。"""
    import Plot_Averaging as averaging

    path = averaging.averaging_result_path(averaging.cached_data(datakey), int(m_val), int(amp), ratio)
    return averaging, result_in_psi(averaging.load_result(path))


def set_phase_ylim(ax, base, result, key, ratio_key, add_shade):
    values = base.estimate_band_values(result, key, key + "_err", add_shade)
    if key == "y1":
        ax.set_ylim(0.9, 1.1)
    elif key == "y2":
        base.set_dynamic_ylim(ax, values, base.BX_MIN_YSPAN[ratio_key])
    elif ratio_key in base.PSI_FIXED_YHALFSPAN:
        base.set_centered_ylim(ax, values, base.PSI_FIXED_YHALFSPAN[ratio_key])
    else:
        base.set_dynamic_ylim(ax, values, base.CX_MIN_YSPAN[ratio_key])
    base.set_three_yticks(ax)


def plot_phase_equations(base, result, ratio, savedir, add_shade, png_output, eps_output,
                         ylim_reference=None, panels=("Psi",)):
    # README.mdの単位説明: y1/y2/y3は既にcycles/year。
    result = result_in_psi(result)
    ylim_base, ylim_result = ylim_reference if ylim_reference is not None else (base, result)
    x = np.asarray(result["x"], dtype=float) / TWOPI
    ratio_key = base.ratio_key_from_pair(ratio)
    saved = []
    for panel, key, phase in (
        ("Annual", "y1", r"\theta_a"),
        ("ENSO", "y2", r"\theta_e"),
        ("Psi", "y3", psi_symbol(ratio)),
    ):
        if panel not in panels:
            continue
        fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
        base.plot_estimate(ax, x, result, key, key + "_err", add_shade)
        setup_phase_axis(ax, phase_xlabel(ratio))
        ax.set_title(phase_title(ratio), fontsize=LABEL_FONTSIZE)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.set_ylabel(phase_rate_ylabel(phase), fontsize=LABEL_FONTSIZE)
        set_phase_ylim(ax, ylim_base, ylim_result, key, ratio_key, add_shade)
        fig.tight_layout()
        saved.extend(base.save_figure(fig, variant_stem(savedir), panel, png_output, eps_output))
    return saved
