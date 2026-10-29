#!/usr/bin/env python3
# coding: utf-8

# %% 
from __future__ import annotations

from plot_config import FIGURE_OUTPUT_DIR

from pathlib import Path
from typing import get_args

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Common import DatasetKey, get_data


plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 10
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["figure.facecolor"] = "white"


TWOPI = 2.0 * np.pi
MONTH_TICKS = np.linspace(0.0, TWOPI, 13)
MONTH_TICKLABELS = ["J.", "F.", "M.", "A.", "M.", "J.", "J.", "A.", "S.", "O.", "N.", "D.", "J."]

#* これで指定した結果を使います。
DATAKEY = "ERSSTv5_2026"
AMP_val = int(30)
Q_M_val = int(3)

#* ODE properties
ODE_DT = 1.0 / 12.0
ODE_YEARS = 10.0
ODE_STEPS = int(round(ODE_YEARS / ODE_DT)) + 1
ODE_INITIAL_PHASES = (
    (0.0, 0.0),
    (0.0, 0.5 * np.pi),
    (0.0, np.pi),
    (0.0, 1.5 * np.pi),
)

#* mpl properties
ODE_LINE_COLORS = ("#00C853", "#FFB000", "#00AEEF", "#FF4FA3")



def amp_rank(amp_val):
    return int(round(amp_val * 100.0 if amp_val <= 1.0 else amp_val))


def q_result_path(datas, Q_M_val, AMP_val):
    return (
        datas["all"]["dataset_dir"]
        / "SRC_py"
        / "Interim_data_Q"
        / f"Q_result_fixed_M{int(Q_M_val)}_amp{amp_rank(AMP_val)}.npz"
    )


def load_result(path):
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def bilinear_periodic(phi_a, phi_e, z):
    phi_a = float(np.mod(phi_a, TWOPI))
    phi_e = float(np.mod(phi_e, TWOPI))
    n_e, n_a = z.shape
    da = TWOPI / (n_a - 1)
    de = TWOPI / (n_e - 1)

    ia = min(int(phi_a / da), n_a - 2)
    ie = min(int(phi_e / de), n_e - 2)
    ta = (phi_a - ia * da) / da
    te = (phi_e - ie * de) / de

    z00 = z[ie, ia]
    z10 = z[ie, ia + 1]
    z01 = z[ie + 1, ia]
    z11 = z[ie + 1, ia + 1]
    return float((1 - ta) * (1 - te) * z00 + ta * (1 - te) * z10 + (1 - ta) * te * z01 + ta * te * z11)


def simulate_q_ode(q_result, initial_phase, step_count=ODE_STEPS, dt=ODE_DT):
    q = np.asarray(q_result["Q"], dtype=float)

    t = np.arange(max(2, int(step_count)), dtype=float) * float(dt)
    phi = np.empty((t.size, 2), dtype=float)
    phi[0] = np.asarray(initial_phase, dtype=float).reshape(2)
    dt = np.diff(t)

    for i, dt_i in enumerate(dt):
        phi_a, phi_e = phi[i]
        drift = np.array([
            bilinear_periodic(phi_a, phi_e, q[0]),
            bilinear_periodic(phi_a, phi_e, q[1]),
        ])
        phi[i + 1] = phi[i] + drift * dt_i

    return phi[:, 0], phi[:, 1]


def plot_wrapped_trajectory(ax, phi_a, phi_e, color, label):
    x = np.mod(np.asarray(phi_a, dtype=float), TWOPI)
    y = np.mod(np.asarray(phi_e, dtype=float), TWOPI)
    segments = []

    for i in range(x.size - 1):
        sx = float(x[i])
        sy = float(y[i])
        dx = float(x[i + 1] - sx)
        dy = float(y[i + 1] - sy)
        if dx > np.pi:
            dx -= TWOPI
        elif dx < -np.pi:
            dx += TWOPI
        if dy > np.pi:
            dy -= TWOPI
        elif dy < -np.pi:
            dy += TWOPI

        rem_x = dx
        rem_y = dy
        for _ in range(4):
            end_x = sx + rem_x
            end_y = sy + rem_y
            tx = 1.0
            ty = 1.0
            cross_x = False
            cross_y = False

            if end_x < 0.0 and rem_x < 0.0:
                tx = (0.0 - sx) / rem_x
                cross_x = True
            elif end_x > TWOPI and rem_x > 0.0:
                tx = (TWOPI - sx) / rem_x
                cross_x = True
            if end_y < 0.0 and rem_y < 0.0:
                ty = (0.0 - sy) / rem_y
                cross_y = True
            elif end_y > TWOPI and rem_y > 0.0:
                ty = (TWOPI - sy) / rem_y
                cross_y = True

            t = min(tx, ty, 1.0)
            if t <= 1e-12:
                if cross_x and tx <= 1e-12:
                    sx = 0.0 if rem_x > 0.0 else TWOPI
                if cross_y and ty <= 1e-12:
                    sy = 0.0 if rem_y > 0.0 else TWOPI
                continue

            bx = sx + rem_x * t
            by = sy + rem_y * t
            bx_plot = float(np.clip(bx, 0.0, TWOPI))
            by_plot = float(np.clip(by, 0.0, TWOPI))
            if abs(bx_plot - sx) > 1e-12 or abs(by_plot - sy) > 1e-12:
                segments.append([(sx, sy), (bx_plot, by_plot)])

            if t >= 1.0:
                break

            rem_x *= 1.0 - t
            rem_y *= 1.0 - t
            sx = bx_plot
            sy = by_plot
            if cross_x and abs(t - tx) <= 1e-12:
                sx = 0.0 if rem_x > 0.0 else TWOPI
            if cross_y and abs(t - ty) <= 1e-12:
                sy = 0.0 if rem_y > 0.0 else TWOPI

    if segments:
        lines = LineCollection(
            segments,
            colors=color,
            linewidths=0.8,
            alpha=1.0,
            label=label,
            zorder=4,
        )
        ax.add_collection(lines)

    ax.scatter(x[0], y[0], s=12, marker="o", facecolors=color, edgecolors="white", linewidth=0.45, zorder=5)
    ax.scatter(x[-1], y[-1], s=16, marker="x", color="white", linewidth=1.4, zorder=5)
    ax.scatter(x[-1], y[-1], s=13, marker="x", color=color, linewidth=0.9, zorder=6)


def save_figure(fig, savedir, png_output=True, eps_output=False):
    output_base = Path(savedir)
    if output_base.suffix in {".png", ".eps"}:
        output_base = output_base.with_suffix("")
    if not output_base.is_absolute():
        output_base = FIGURE_OUTPUT_DIR / output_base
    output_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_base.with_suffix(".png")
    eps_path = output_base.with_suffix(".eps")

    if png_output:
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        print(f"saved figure: {png_path}")
    if eps_output:
        fig.savefig(eps_path, format="eps", bbox_inches="tight")
        print(f"saved figure: {eps_path}")


def main(
        DATAKEY=DATAKEY,
        AMP_val=AMP_val,
        Q_M_val=Q_M_val,
        ODE_STEPS=ODE_STEPS,
        ODE_INITIAL_PHASES=ODE_INITIAL_PHASES,
        ODE_LINE_COLORS=ODE_LINE_COLORS,
        savedir="./",
        png_output=True,
        eps_output=False,
):
    datas = get_data(key=DATAKEY)
    q_result = load_result(q_result_path(datas, Q_M_val, AMP_val))

    xx = np.asarray(q_result["phi1"], dtype=float)
    yy = np.asarray(q_result["phi2"], dtype=float)
    # Q は rad/year で保存されているので、等高線表示だけ 2π で割って yr^{-1} に変換する。
    Q = np.asarray(q_result["Q"], dtype=float) / TWOPI
    omega_e = float(np.asarray(q_result["omega_freq_per_year"], dtype=float)[1])
    Q_e = Q[1]

    colorbar_half_width = 0.15
    levels = np.linspace(
        omega_e - colorbar_half_width,
        omega_e + colorbar_half_width,
        17,
    )

    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    im = ax.contourf(
        xx,
        yy,
        Q_e,
        levels=levels,
        cmap="seismic",
        vmin=levels[0],
        vmax=levels[-1],
        extend="both",
    )
    ax.contour(
        xx,
        yy,
        Q_e,
        levels=levels,
        colors="0.25",
        linewidths=0.15,
        alpha=1.0,
    )

    q_e_min = float(np.nanmin(Q_e))
    q_e_max = float(np.nanmax(Q_e))
    if q_e_min <= 0.0 <= q_e_max:
        ax.contour(xx, yy, Q_e, levels=[0.0], colors="0.15", linewidths=0.25, alpha=1.0)

    for i, initial_phase in enumerate(ODE_INITIAL_PHASES):
        phi_a_ode, phi_e_ode = simulate_q_ode(
            q_result,
            initial_phase=initial_phase,
            step_count=ODE_STEPS,
        )
        plot_wrapped_trajectory(
            ax,
            phi_a_ode,
            phi_e_ode,
            color=ODE_LINE_COLORS[i % len(ODE_LINE_COLORS)],
            label="_nolegend_",
        )

    phase_ticks = np.linspace(0.0, TWOPI, 5)
    # 軌道の計算・描画座標はradのまま、目盛りは2piで割って表示する。
    phase_ticklabels = ["0", "0.25", "0.5", "0.75", "1"]
    ax.set_xlim(0.0, TWOPI)
    ax.set_ylim(0.0, TWOPI)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(MONTH_TICKS)
    ax.set_yticks(phase_ticks)
    ax.set_xticklabels(MONTH_TICKLABELS)
    ax.set_yticklabels(phase_ticklabels)
    ax.tick_params(axis="x", labelsize=8, pad=2)
    ax.set_xlabel(r"$\frac{\theta_a}{2\pi}$ (Month)")
    ax.set_ylabel(r"$\frac{\theta_e}{2\pi}$")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(im, cax=cax, label=r"$\frac{Q_e(\theta_a,\theta_e)}{2\pi}$ [cycle / year]")
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    cbar.ax.axhline(omega_e, color="black", linewidth=1.2)

    save_figure(fig, savedir, png_output=png_output, eps_output=eps_output)
    plt.close(fig)

if __name__ == "__main__":
    main(savedir=f"Q_e_ODE_amp{AMP_val}_M{Q_M_val}_{DATAKEY}_ODEoverlay",
         png_output=True, eps_output=True)
