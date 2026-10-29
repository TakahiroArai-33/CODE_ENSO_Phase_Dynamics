#!/usr/bin/env python3
# coding: utf-8
"""Langevin軌道によるEl Niño / La Niñaの初通過時間のポンチ絵。

プレゼン素材/ランジュバンシミュレーション_経過時間測定のポンチ絵 の
推定Qe版を元に、Fig6の4月・9月の単独図をカラーバー付きで出力する。
計算はrad・year、背景のQeのみcycles/year（README.mdを参照）。
"""

from __future__ import annotations

from plot_config import FIGURE_OUTPUT_DIR

import argparse
import calendar
import os
from pathlib import Path


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator

from Common import get_data


PAPERFIGURE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = FIGURE_OUTPUT_DIR / "素材" / "LangevinPassage"
DATAKEY = "ERSSTv5_2026"
AMP_VAL = 30
Q_M_VAL = 3
START_MONTHS = (4, 9)
N_TRAJECTORIES = 20
SDE_DT_YEARS = 1.0 / 120.0
SDE_MAX_YEARS = 10.0
TWOPI = 2.0 * np.pi

# (ファイル名, 表示名, 開始位相 [rad], 到達位相 [rad, unwrap], seed)
EVENTS = (
    ("ElNino", "El Niño", 13.0 * np.pi / 8.0, TWOPI + 3.0 * np.pi / 8.0, 2300),
    ("LaNina", "La Niña", 3.0 * np.pi / 8.0, 9.0 * np.pi / 8.0, 2400),
)
MONTH_TICKLABELS = ["J.", "F.", "M.", "A.", "M.", "J.", "J.", "A.", "S.", "O.", "N.", "D.", "J."]


def load_model():
    dataset_dir = get_data(key=DATAKEY)["all"]["dataset_dir"]
    path = dataset_dir / "SRC_py" / "Interim_data_Q" / f"Q_result_fixed_M{Q_M_VAL}_amp{AMP_VAL}.npz"
    with np.load(path) as result:
        theta_a = np.asarray(result["phi1"], dtype=float)
        theta_e = np.asarray(result["phi2"], dtype=float)
        q_rad = np.asarray(result["Q"], dtype=float)  # rad/year
        covariance = np.asarray(result["noise_cov"], dtype=float)  # E=HH^T, rad^2/year
        omega_e = float(result["omega_freq_per_year"][1])  # cycles/year

    axes = (theta_a[0, :], theta_e[:, 0])
    drift_a = RegularGridInterpolator(axes, q_rad[0].T)
    drift_e = RegularGridInterpolator(axes, q_rad[1].T)
    covariance = 0.5 * (covariance + covariance.T)
    noise_cholesky = np.linalg.cholesky(covariance + 1.0e-14 * np.eye(2))
    return theta_a, theta_e, q_rad[1] / TWOPI, omega_e, drift_a, drift_e, noise_cholesky


def simulate_sde_until_passage(
    initial_phase, target_theta_e, drift_a, drift_e, noise_cholesky, seed,
    dt_years=SDE_DT_YEARS, max_years=SDE_MAX_YEARS,
):
    """Euler–Maruyama法で最初の上向き到達まで計算し、軌道と到達時間[year]を返す。"""
    rng = np.random.default_rng(seed)
    trajectory = [np.asarray(initial_phase, dtype=float)]
    for step in range(int(round(max_years / dt_years))):
        previous = trajectory[-1]
        point = tuple(np.mod(previous, TWOPI))
        drift = np.array([
            drift_a(point).item(),
            drift_e(point).item(),
        ])
        noise = np.sqrt(dt_years) * (rng.standard_normal(2) @ noise_cholesky.T)
        current = previous + drift * dt_years + noise
        # modulo後の位相では判定しない。El Niñoでは2piを越えた後に到達する。
        if current[1] >= target_theta_e:
            fraction = (target_theta_e - previous[1]) / (current[1] - previous[1])
            fraction = np.clip(fraction, 0.0, 1.0)
            current = previous + fraction * (current - previous)
            current[1] = target_theta_e
            trajectory.append(current)
            return np.asarray(trajectory), (step + fraction) * dt_years
        trajectory.append(current)
    raise RuntimeError(f"SDE did not reach within {max_years:g} years (seed={seed}).")


def wrap_trajectory_for_plot(trajectory):
    """周期境界をまたぐ線をNaNで切る（数値計算ではunwrapした位相を保持）。"""
    wrapped = np.mod(trajectory, TWOPI)
    jumps = np.flatnonzero(np.any(np.abs(np.diff(wrapped, axis=0)) > np.pi, axis=1)) + 1
    return np.insert(wrapped, jumps, np.nan, axis=0)


def seed_for_month(base_seed, month):
    # 元notebookの3月・9月の乱数系列をそのまま使う。
    return base_seed + {3: 0, 9: 1000}.get(month, 10000 + month * 1000)


def plot_panel(fig, ax, model, trajectories, title, start_e, end_e, show_colorbar=True):
    theta_a, theta_e, q_e, omega_e, *_ = model
    levels = np.linspace(omega_e - 0.15, omega_e + 0.15, 17)
    # EPSでも半透明の軌道を保つため、背景と軌道のみ300dpiでラスタ化する。
    # 軸・文字・開始点・到達点はベクトルとして保存する。
    ax.set_rasterization_zorder(3)
    im = ax.contourf(
        theta_a, theta_e, q_e, levels=levels, cmap="seismic",
        vmin=levels[0], vmax=levels[-1], extend="both", zorder=0,
    )
    ax.contour(theta_a, theta_e, q_e, levels=levels, colors="0.25", linewidths=0.15, zorder=1)
    for phase in (start_e, end_e):
        ax.axhline(phase % TWOPI, color="black", linewidth=1.0, zorder=7)
    for trajectory in trajectories:
        wrapped = wrap_trajectory_for_plot(trajectory)
        ax.plot(wrapped[:, 0], wrapped[:, 1], color="black", linewidth=0.55, alpha=0.60, zorder=2)
    start = trajectories[0][0] % TWOPI
    ends = np.array([trajectory[-1] for trajectory in trajectories]) % TWOPI
    ax.scatter(*start, s=30, color="yellow", edgecolor="black", linewidth=0.7, zorder=10)
    ax.scatter(ends[:, 0], ends[:, 1], s=16, color="limegreen", edgecolor="black", linewidth=0.4, zorder=10)
    ax.set(
        xlim=(0, TWOPI), ylim=(0, TWOPI),
        xlabel=r"$\frac{\theta_a}{2\pi}$ (Month)", ylabel=r"$\frac{\theta_e}{2\pi}$",
        title=title,
    )
    ax.set_aspect("equal", adjustable="box")
    # SDE計算・描画座標はradのまま、軸の表示だけcyclesに揃える。
    phase_ticks = np.linspace(0, TWOPI, 5)
    phase_ticklabels = ["0", "0.25", "0.5", "0.75", "1"]
    ax.set_xticks(np.linspace(0, TWOPI, 13), MONTH_TICKLABELS)
    ax.tick_params(axis="x", labelsize=8, pad=2)
    ax.set_yticks(phase_ticks, phase_ticklabels)
    if show_colorbar:
        cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.08)
        cbar = fig.colorbar(im, cax=cax, label=r"$\frac{Q_e(\theta_a,\theta_e)}{2\pi}$ [cycle / year]")
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        cbar.ax.axhline(omega_e, color="black", linewidth=1.2)


def save_figure(fig, output_dir, stem):
    fig.tight_layout()
    for extension in ("png", "eps"):
        path = output_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"saved figure: {path}")
    plt.close(fig)


def main(months=START_MONTHS, n_trajectories=N_TRAJECTORIES, output_dir=OUTPUT_DIR):
    if not months or any(month not in range(1, 13) for month in months):
        raise ValueError("months must contain month numbers from 1 to 12.")
    if len(set(months)) != len(months):
        raise ValueError("months must be unique.")
    if not 1 <= n_trajectories <= 1000:
        raise ValueError("n_trajectories must be between 1 and 1000 (independent seeds).")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "mathtext.fontset": "stix", "font.family": "Times New Roman", "font.size": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10, "figure.facecolor": "white",
    })
    model = load_model()
    *_, drift_a, drift_e, noise_cholesky = model

    for event, label, start_e, end_e, base_seed in EVENTS:
        stem = f"Langevin_passage_{event}_amp{AMP_VAL}_M{Q_M_VAL}_{DATAKEY}_estimated"
        for month in months:
            paths_and_times = [
                simulate_sde_until_passage(
                    [(month - 1) * TWOPI / 12.0, start_e], end_e,
                    drift_a, drift_e, noise_cholesky,
                    seed_for_month(base_seed, month) + index,
                )
                for index in range(n_trajectories)
            ]
            trajectories = [item[0] for item in paths_and_times]
            times_months = np.array([item[1] for item in paths_and_times]) * 12.0
            title = f"{label}: {calendar.month_name[month]} start"
            print(f"{event} estimated month={month:02d}: n={n_trajectories}, "
                  f"sample passage time mean={times_months.mean():.3f} months, "
                  f"range={times_months.min():.3f}..{times_months.max():.3f} months")
            fig, ax = plt.subplots(figsize=(3.5, 2.9))
            plot_panel(fig, ax, model, trajectories, title, start_e, end_e, show_colorbar=True)
            save_figure(fig, output_dir, f"{stem}_month{month:02d}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--months", nargs="+", type=int, default=START_MONTHS, help="開始月（既定: 4 9）")
    parser.add_argument("--trajectories", type=int, default=N_TRAJECTORIES, help="各パネルの軌道数（既定: 20）")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    main(args.months, args.trajectories, args.output_dir)
