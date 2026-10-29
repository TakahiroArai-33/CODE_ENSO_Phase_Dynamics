#!/usr/bin/env python3
# coding: utf-8

# %%%
from __future__ import annotations

from plot_config import FIGURE_OUTPUT_DIR

import os
from pathlib import Path


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Bbox

from Common import (
    get_data,
    get_phase_series,
    get_SST_phase_composite,
)


plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 10
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["figure.facecolor"] = "white"


DATAKEY = "ERSSTv5_2026"
S_VAL = 8
AMP_VAL = 0.3
OUTPUT_DPI = 300
INDIVIDUAL_MARGIN_MM = 1.0
TWOPI = 2.0 * np.pi
PHASE_COLORS = ["red", "lime", "navy", "magenta", "blue", "darkgreen", "gold", "cornflowerblue"]
MAP_LON = (100.0, 300.0)
MAP_LAT = (-40.0, 40.0)
AXIS_TICK_SIZE = 12
COLORBAR_TICK_SIZE = AXIS_TICK_SIZE
OUTPUT_DIR = FIGURE_OUTPUT_DIR
INDIVIDUAL_OUTPUT_DIR = OUTPUT_DIR / "素材" / "PhaseComposit個別"


def amp_label() -> int:
    return int(round(AMP_VAL * 100.0 if AMP_VAL <= 1.0 else AMP_VAL))


def phase_mask(phase: np.ndarray, i: int, valid: np.ndarray) -> np.ndarray:
    width = TWOPI / S_VAL
    center = width * i
    return valid & (np.abs(np.mod(phase - center + np.pi, TWOPI) - np.pi) <= 0.5 * width)


def composite_fields(composite_data: dict[str, object]) -> dict[int, np.ma.MaskedArray]:
    return {
        phase_number - 1: np.ma.asarray(field, dtype=float)
        for phase_number, field in composite_data["composite"]
    }


def draw_map(ax, lon, lat, field, title: str, vlim: float):
    lon2d, lat2d = np.meshgrid(lon, lat)
    land_mask = np.ma.getmaskarray(field)
    image = ax.pcolormesh(lon, lat, field, vmin=-vlim, vmax=vlim, shading="auto", cmap="seismic")
    if np.any(land_mask) and not np.all(land_mask):
        ax.contour(lon2d, lat2d, land_mask.astype(float), levels=[0.5], colors="k", linewidths=0.45)
    ax.set_ylim(*MAP_LAT)
    ax.set_yticks([-40, -20, 0, 20, 40])
    ax.set_yticklabels(
        [r"$40^\circ$S", r"$20^\circ$S", r"$0^\circ$", r"$20^\circ$N", r"$40^\circ$N"],
        fontsize=AXIS_TICK_SIZE,
    )
    ax.set_xlim(*MAP_LON)
    ax.set_xticks([120, 180, 240, 300])
    ax.set_xticklabels(
        [r"$120^\circ$E", r"$180^\circ$", r"$120^\circ$W", r"$60^\circ$W"],
        fontsize=AXIS_TICK_SIZE,
    )
    ax.set_title(title, fontsize=12)
    return image


def add_phase_colorbar(fig, image, ax, vlim: float):
    cbar_ticks = [-vlim, 0.0, vlim]
    cbar = fig.colorbar(image, ax=ax, ticks=cbar_ticks, fraction=0.046, pad=0.02)
    cbar.set_ticklabels([f"{tick:.1f}" for tick in cbar_ticks])
    cbar.ax.tick_params(labelsize=COLORBAR_TICK_SIZE, length=2, pad=1)
    return cbar


def draw_phase_composite_grid(composite_data: dict[str, object], vlim: float):
    lon = np.asarray(composite_data["lon"], dtype=float)
    lat = np.asarray(composite_data["lat"], dtype=float)
    fields = composite_fields(composite_data)

    fig, axes = plt.subplots(figsize=(6.5, 6.8), nrows=4, ncols=2)
    colorbar_axes = []

    for i in range(S_VAL):
        row = i // 2
        col = i % 2
        image = draw_map(axes[row, col], lon, lat, fields[i], f"j={i + 1:d}", vlim)
        cbar = add_phase_colorbar(fig, image, axes[row, col], vlim)
        colorbar_axes.append(cbar.ax)

    fig.tight_layout(h_pad=0.8, w_pad=0.7)
    fig.canvas.draw()
    return fig, axes, colorbar_axes, lon, lat, fields


def panel_layouts_from_grid(fig, axes, colorbar_axes):
    renderer = fig.canvas.get_renderer()
    layouts = []
    for ax, cax in zip(axes.ravel(), colorbar_axes):
        panel_bbox = Bbox.union([ax.get_tightbbox(renderer), cax.get_tightbbox(renderer)])
        ax_bbox = ax.get_window_extent(renderer)
        cax_bbox = cax.get_window_extent(renderer)

        panel = panel_bbox.transformed(fig.dpi_scale_trans.inverted())
        ax_box = ax_bbox.transformed(fig.dpi_scale_trans.inverted())
        cax_box = cax_bbox.transformed(fig.dpi_scale_trans.inverted())

        def relative(box):
            return [
                (box.x0 - panel.x0) / panel.width,
                (box.y0 - panel.y0) / panel.height,
                box.width / panel.width,
                box.height / panel.height,
            ]

        layouts.append((panel.width, panel.height, relative(ax_box), relative(cax_box)))
    return layouts


def save_individual_phase_composites(
        output_stem: str,
        lon: np.ndarray,
        lat: np.ndarray,
        fields: dict[int, np.ma.MaskedArray],
        vlim: float,
        layouts,
) -> None:
    INDIVIDUAL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    margin = INDIVIDUAL_MARGIN_MM / 25.4
    for i, (fig_width, fig_height, ax_rect, cax_rect) in enumerate(layouts):
        output_width = fig_width + margin
        output_height = fig_height + margin

        def with_top_right_margin(rect):
            return [
                rect[0] * fig_width / output_width,
                rect[1] * fig_height / output_height,
                rect[2] * fig_width / output_width,
                rect[3] * fig_height / output_height,
            ]

        fig = plt.figure(figsize=(output_width, output_height))
        ax = fig.add_axes(with_top_right_margin(ax_rect))
        cax = fig.add_axes(with_top_right_margin(cax_rect))
        image = draw_map(ax, lon, lat, fields[i], f"j={i + 1:d}", vlim)
        cbar_ticks = [-vlim, 0.0, vlim]
        cbar = fig.colorbar(image, cax=cax, ticks=cbar_ticks)
        cbar.set_ticklabels([f"{tick:.1f}" for tick in cbar_ticks])
        cbar.ax.tick_params(labelsize=COLORBAR_TICK_SIZE, length=2, pad=1)

        output = INDIVIDUAL_OUTPUT_DIR / f"{output_stem}_j{i + 1:d}.eps"
        fig.savefig(output, dpi=OUTPUT_DPI, format="eps")
        fig.savefig(output.with_suffix(".png"), dpi=OUTPUT_DPI, format="png")
        plt.close(fig)
        print(f"saved figure: {output}")


def plot_phase_composite_grid(
        composite_data: dict[str, object], output_name: str, vlim: float,
        save_grid: bool = False,
) -> None:
    # グリッドは個別パネルの寸法計算に使う。論文で使う個別図だけを保存する。
    fig, axes, colorbar_axes, lon, lat, fields = draw_phase_composite_grid(composite_data, vlim)
    layouts = panel_layouts_from_grid(fig, axes, colorbar_axes)

    if save_grid:
        output = OUTPUT_DIR / output_name
        fig.savefig(output, dpi=OUTPUT_DPI, format="eps")
        fig.savefig(output.with_suffix(".png"), dpi=OUTPUT_DPI, format="png")
        print(f"saved figure: {output}")
    plt.close(fig)

    output_stem = Path(output_name).with_suffix("").name
    save_individual_phase_composites(output_stem, lon, lat, fields, vlim, layouts)


def plot_center(
        ax,
        z: np.ndarray,
        phase: np.ndarray,
        valid: np.ndarray,
        label: str,
        show_phase_guides: bool = False,
        show_origin_axes: bool = True,
        trajectory_color: str = "0.78",
        gray_point_size: float = 8.0,
) -> None:
    lim = np.nanpercentile(np.abs(z), 99.0) * 1.08
    if not np.isfinite(lim) or lim <= 0.0:
        lim = 1.0

    ax.plot(np.real(z), np.imag(z), color=trajectory_color, linewidth=0.25, zorder=0)
    ax.scatter(np.real(z[~valid]), np.imag(z[~valid]), s=gray_point_size, color="0.65", linewidths=0, zorder=2)
    if show_origin_axes:
        ax.axhline(0.0, color="0.25", linewidth=0.7, zorder=1)
        ax.axvline(0.0, color="0.25", linewidth=0.7, zorder=1)

    if show_phase_guides:
        for i in range(S_VAL // 2):
            th = TWOPI / S_VAL * (i + 0.5)
            x = np.array([-lim, lim])
            ax.plot(x, np.tan(th) * x, color="k", linestyle="--", linewidth=0.5, zorder=1)

    for i in range(S_VAL):
        mask = phase_mask(phase, i, valid)
        ax.scatter(np.real(z[mask]), np.imag(z[mask]), s=8, color=PHASE_COLORS[i], linewidths=0, zorder=3)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(rf"Re ${label}$", fontsize=13, labelpad=1)
    ax.set_ylabel(rf"Im ${label}$", fontsize=13, labelpad=1)


def save_center_plot(
        name: str,
        z: np.ndarray,
        phase: np.ndarray,
        valid: np.ndarray,
        label: str,
        show_phase_guides: bool = False,
        show_origin_axes: bool = True,
        trajectory_color: str = "0.78",
        gray_point_size: float = 8.0,
) -> None:
    fig, ax = plt.subplots(figsize=(3.0, 3.0))
    plot_center(
        ax, z, np.mod(phase, TWOPI), valid, label,
        show_phase_guides=show_phase_guides,
        show_origin_axes=show_origin_axes,
        trajectory_color=trajectory_color,
        gray_point_size=gray_point_size,
    )
    fig.tight_layout(pad=0.35)
    output = OUTPUT_DIR / f"ComplexPlane_{name}_amp{amp_label():d}_fixed.eps"
    fig.savefig(output, dpi=OUTPUT_DPI, format="eps", bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=OUTPUT_DPI, format="png", bbox_inches="tight")
    plt.close(fig)
    print(f"saved figure: {output}")


def main() -> None:
    datas = get_data(key=DATAKEY)
    phase_datas = get_phase_series(datas)
    annual_enso = datas["annual_enso"]
    mode_numbers = annual_enso["mode_info"]["numbers"]

    enso_composite = get_SST_phase_composite(key1=DATAKEY, key2="ENSO", S=S_VAL, amp=AMP_VAL)
    annual_composite = get_SST_phase_composite(key1=DATAKEY, key2="Annual", S=S_VAL, amp=AMP_VAL)

    amp = amp_label()
    plot_phase_composite_grid(
        enso_composite,
        f"PhaseComposite_ENSO_amp{amp:d}_fixed.eps",
        vlim=2.5,
    )
    plot_phase_composite_grid(
        annual_composite,
        f"PhaseComposite_Annual_amp{amp:d}_fixed.eps",
        vlim=3.0,
    )

    save_center_plot(
        f"g{mode_numbers['annual']:d}",
        np.asarray(annual_enso["Z_annual"]),
        phase_datas["phase_annual"],
        np.asarray(annual_composite["phase_index"]) > 0,
        label="g_a",
        show_phase_guides=True,
        show_origin_axes=False,
        trajectory_color="0.90",
    )
    save_center_plot(
        f"g{mode_numbers['enso']:d}",
        np.asarray(annual_enso["Z_enso"]),
        phase_datas["phase_enso"],
        np.asarray(enso_composite["phase_index"]) > 0,
        label="g_e",
        gray_point_size=2.5,
    )


if __name__ == "__main__":
    main()
