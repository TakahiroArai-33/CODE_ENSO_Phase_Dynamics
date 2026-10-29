
from plot_config import FIGURE_OUTPUT_DIR
#!/usr/bin/env python3
# coding: utf-8

# %% 

from pathlib import Path
from typing import Literal, get_args

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from Common import DatasetKey, get_data

plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.titlesize"] = 9
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 9
plt.rcParams["figure.facecolor"] = "white"


RatioKey = Literal["1:2", "1:3", "1:4", "1:5", "2:5"]


DATAKEY = "ERSSTv5_2026"
AMP_val = int(30)
DA_M_val = int(3)
RATIO_val: RatioKey = "1:3"



TWOPI = 2.0 * np.pi
_DATA_CACHE = {}


PANEL_FIGSIZE = (2.4, 1.7)
ESTIMATE_LINE_COLOR = "blue"
ESTIMATE_SHADE_COLOR = "#bfe6ff"
ESTIMATE_SHADE_ALPHA = 1.0

BX_MIN_YSPAN: dict[RatioKey, float] = {
    "1:2": 0.10,
    "1:3": 0.16,
    "1:4": 0.08,
    "1:5": 0.08,
    "2:5": 0.08,
}
CX_MIN_YSPAN: dict[RatioKey, float] = {
    "1:2": 0.20,
    "1:3": 0.45,
    "1:4": 0.16,
    "1:5": 0.28,
    "2:5": 0.18,
}
PSI_FIXED_YHALFSPAN: dict[RatioKey, float] = {
    "1:2": 0.10,
    "1:4": 0.10,
    "1:5": 0.10,
    "2:5": 0.10,
}


def ratio_filename(ratio):
    return str(ratio).replace(":", "to")


def ratio_pair(ratio):
    return tuple(int(value) for value in str(ratio).split(":"))


def result_ratio_pair(result, fallback_ratio):
    if "ratio" in result:
        ratio = np.asarray(result["ratio"], dtype=int).reshape(-1)
        if ratio.size >= 2:
            return int(ratio[0]), int(ratio[1])
    return ratio_pair(fallback_ratio)


def ratio_key_from_pair(ratio):
    m_val, n_val = ratio
    return f"{m_val:d}:{n_val:d}"


def delta_phi_expr_from_pair(ratio):
    m_val, n_val = ratio
    if m_val == 1:
        return rf"\theta_a - {n_val:d}\theta_e"
    return rf"\theta_a - \frac{{{n_val:d}}}{{{m_val:d}}}\theta_e"


def reverse_delta_phi_expr_from_pair(ratio):
    m_val, n_val = ratio
    if m_val == 1:
        return rf"{n_val:d}\theta_e - \theta_a"
    return rf"\frac{{{n_val:d}}}{{{m_val:d}}}\theta_e - \theta_a"


def delta_phi_expr(ratio):
    return delta_phi_expr_from_pair(ratio_pair(ratio))


def reverse_delta_phi_expr(ratio):
    return reverse_delta_phi_expr_from_pair(ratio_pair(ratio))


def phase_rate_ylabel(phase_expr):
    return rf"$\frac{{d}}{{dt}}\left({phase_expr}\right)$ [yr$^{{-1}}$]"


def phase_density_ylabel(phase_expr):
    return rf"$P\left({phase_expr}\right)$"


def setup_phase_axis(ax, xlabel):
    ax.set_xlim(0.0, TWOPI)
    ax.set_xticks(np.linspace(0.0, TWOPI, 5))
    ax.set_xticklabels([r"$0$", r"$0.5\pi$", r"$\pi$", r"$1.5\pi$", r"$2\pi$"])
    ax.set_xlabel(xlabel)
    ax.grid(linestyle="--", alpha=0.6)


def set_dynamic_ylim(ax, values, min_span, margin_ratio=0.15):
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return

    ymin = float(np.nanmin(finite_values))
    ymax = float(np.nanmax(finite_values))
    center = 0.5 * (ymin + ymax)
    span = ymax - ymin
    span = max(span * (1.0 + 2.0 * margin_ratio), float(min_span))
    ax.set_ylim(center - 0.5 * span, center + 0.5 * span)


def set_centered_ylim(ax, values, half_span):
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return

    center = 0.5 * (float(np.min(finite_values)) + float(np.max(finite_values)))
    ax.set_ylim(center - float(half_span), center + float(half_span))


def set_three_yticks(ax):
    ymin, ymax = ax.get_ylim()
    ax.set_yticks(np.linspace(ymin, ymax, 3))


def estimate_band_values(result, y_key, err_key, add_shade):
    y = np.asarray(result[y_key], dtype=float)
    if not add_shade or err_key not in result:
        return y
    y_err = np.asarray(result[err_key], dtype=float)
    return np.concatenate((y - y_err, y + y_err))


def plot_estimate(ax, x, result, y_key, err_key, add_shade):
    y = np.asarray(result[y_key], dtype=float)
    if add_shade and err_key in result:
        y_err = np.asarray(result[err_key], dtype=float)
        ax.fill_between(
            x,
            y - y_err,
            y + y_err,
            color=ESTIMATE_SHADE_COLOR,
            alpha=ESTIMATE_SHADE_ALPHA,
            linewidth=0.0,
            zorder=1,
        )
    ax.plot(x, y, color=ESTIMATE_LINE_COLOR, linewidth=1.1, zorder=2)


def load_result(path):
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def cached_data(key):
    if key not in _DATA_CACHE:
        _DATA_CACHE[key] = get_data(key=key)
    return _DATA_CACHE[key]


def direct_result_path(datas, M_val, amp_rank, ratio):
    return (
        datas["all"]["dataset_dir"]
        / "SRC_py"
        / "Interim_data_DirectAveraging"
        / f"DirectAvaraging_result_M{M_val:d}_amp{amp_rank:d}_ratio{ratio}.npz"
    )


def output_base_path(savedir, panel_name=None):
    output_base = Path(savedir)
    if output_base.suffix in {".png", ".eps"}:
        output_base = output_base.with_suffix("")
    if not output_base.is_absolute():
        output_base = FIGURE_OUTPUT_DIR / output_base
    if panel_name is not None:
        output_base = output_base.with_name(f"{output_base.name}_{panel_name}")
    output_base.parent.mkdir(parents=True, exist_ok=True)
    return output_base


def save_figure(fig, savedir, panel_name, png_output, eps_output):
    output_base = output_base_path(savedir, panel_name)
    png_path = output_base.with_suffix(".png")
    eps_path = output_base.with_suffix(".eps")
    saved = []
    if png_output:
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        print(f"saved figure: {png_path}")
        saved.append(png_path)
    if eps_output:
        fig.savefig(eps_path, format="eps", bbox_inches="tight")
        print(f"saved figure: {eps_path}")
        saved.append(eps_path)
    plt.close(fig)
    return saved


def main(
        DATAKEY=DATAKEY,
        AMP_val=AMP_val,
        DA_M_val=DA_M_val,
        RATIO_val: RatioKey = RATIO_val,
        savedir="./",
        png_output=True,
        eps_output=False,
        add_shade=True,
):
    datas = cached_data(DATAKEY)
    amp_rank = int(round(AMP_val * 100.0 if AMP_val <= 1.0 else AMP_val))
    result_path = direct_result_path(datas, int(DA_M_val), amp_rank, RATIO_val)
    if not result_path.exists():
        print(f"skip missing: {result_path}")
        return None

    result = load_result(result_path)
    #* DirectAveraging の y1/y2/y3 は保存時に 2π で割ってあり、すでに yr^{-1} 単位。
    x = np.asarray(result["x"], dtype=float)
    loaded_ratio = result_ratio_pair(result, RATIO_val)
    ratio_key = ratio_key_from_pair(loaded_ratio)
    phase_expr = delta_phi_expr_from_pair(loaded_ratio)
    xlabel = "$" + phase_expr + "$"
    reverse_xlabel = "$" + reverse_delta_phi_expr_from_pair(loaded_ratio) + "$"
    saved = []

    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    plot_estimate(ax, x, result, "y1", "y1_err", add_shade)
    setup_phase_axis(ax, xlabel)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.set_ylim(0.9, 1.1)
    set_three_yticks(ax)
    ax.set_ylabel(r"$\frac{d}{dt}\theta_a$ [yr$^{-1}$]")
    fig.tight_layout()
    saved.extend(save_figure(fig, savedir, "Annual", png_output, eps_output))

    fig, bx = plt.subplots(figsize=PANEL_FIGSIZE)
    plot_estimate(bx, x, result, "y2", "y2_err", add_shade)
    setup_phase_axis(bx, reverse_xlabel)
    bx.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    set_dynamic_ylim(bx, estimate_band_values(result, "y2", "y2_err", add_shade), BX_MIN_YSPAN[ratio_key])
    set_three_yticks(bx)
    bx.set_ylabel(r"$\frac{d}{dt}\theta_e$ [yr$^{-1}$]")
    fig.tight_layout()
    saved.extend(save_figure(fig, savedir, "ENSO", png_output, eps_output))

    fig, cx = plt.subplots(figsize=PANEL_FIGSIZE)
    plot_estimate(cx, x, result, "y3", "y3_err", add_shade)
    setup_phase_axis(cx, xlabel)
    cx.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cx.set_ylabel(phase_rate_ylabel(phase_expr))
    psi_values = estimate_band_values(result, "y3", "y3_err", add_shade)
    if ratio_key in PSI_FIXED_YHALFSPAN:
        set_centered_ylim(cx, psi_values, PSI_FIXED_YHALFSPAN[ratio_key])
    else:
        set_dynamic_ylim(cx, psi_values, CX_MIN_YSPAN[ratio_key])
    set_three_yticks(cx)
    fig.tight_layout()
    saved.extend(save_figure(fig, savedir, "Psi", png_output, eps_output))

    fig, dx = plt.subplots(figsize=PANEL_FIGSIZE)
    bins = result["distribution_bins"]
    dx.hist(
        result["delta_phi_data"],
        bins=bins,
        density=True,
        color="gray",
        edgecolor="k",
        linewidth=0.5,
        alpha=0.5,
    )
    dx.hist(
        result["delta_phi_langevin"],
        bins=bins,
        density=True,
        histtype="step",
        color="blue",
        linewidth=1.0,
    )
    dx.set_xlim(0.0, TWOPI)
    dx.set_ylim(0.0, 0.45)
    set_three_yticks(dx)
    dx.set_xticks([0.0, np.pi, TWOPI])
    dx.set_xticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    dx.set_xlabel(xlabel)
    dx.set_ylabel(phase_density_ylabel(phase_expr))
    dx.grid(linestyle="--", alpha=0.6)
    fig.tight_layout()
    saved.extend(save_figure(fig, savedir, "Hist", png_output, eps_output))
    return saved





if __name__ == "__main__":
    for ratio in get_args(RatioKey):
        main(RATIO_val=ratio, savedir=f"DirectAveraging_amp30_M3_ratio{ratio_filename(ratio)}_{DATAKEY}",
             png_output=True, eps_output=True)
