
from plot_config import FIGURE_OUTPUT_DIR
#!/usr/bin/env python3
# coding: utf-8

# # $Q_e$ の通過時間計算
#
# `SRC/plotting/Plot_Qe.py` と同じ方法で推定結果を読み込み、まず $Q_e(\theta_a,\theta_e)$ を表示します。
# 通過時間の計算は、この読み込み・表示結果を確認した後のセルに追加します。

# ## Setup and parameters
#
# 既定値は `ERSSTv5_2026`、`M=3`、振幅フィルター30%です。`AMP_val=30` は既存の `Plot_Qe.py` の既定値に合わせています。

# In[1]:


import os
import sys
from pathlib import Path
from typing import get_args


#* shellからどのcurrent directoryで実行しても、このscriptと同じ場所のCommon.pyを読みます。
PAPERFIGURE_DIR = Path(__file__).resolve().parent

if str(PAPERFIGURE_DIR) not in sys.path:
    sys.path.insert(0, str(PAPERFIGURE_DIR))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator

from Common import DatasetKey, get_data

plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 10
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["figure.facecolor"] = "white"

DATAKEY: DatasetKey = "ERSSTv5_2026"
AMP_val = int(30)
Q_M_val = int(3)
VALID_DATAKEYS = tuple(get_args(DatasetKey))
if DATAKEY not in VALID_DATAKEYS:
    raise ValueError(f"Unknown DATAKEY {DATAKEY!r}; expected one of {VALID_DATAKEYS}")


# ## Load $Q_e$
#
# `Common.get_data()` で対象データセットの場所を取得し、`Plot_Qe.py` と同じ `Q_result_fixed_M{M}_amp{amp}.npz` を読み込みます。保存されている $Q$ は rad/year のため、$2\pi$ で割って yr$^{-1}$ に変換します。

# In[2]:


def amp_rank(amp_val):
    return int(round(amp_val * 100.0 if amp_val <= 1.0 else amp_val))


def main(
        DATAKEY=DATAKEY,
        AMP_val=AMP_val,
        Q_M_val=Q_M_val,
        SDE_TRIAL_COUNT=5000,
        savedir=FIGURE_OUTPUT_DIR,
        output_variants=("SDE",),
):
    if DATAKEY not in VALID_DATAKEYS:
        raise ValueError(
            f"Unknown DATAKEY {DATAKEY!r}; expected one of {VALID_DATAKEYS}"
        )

    output_dir = Path(savedir)
    if not output_dir.is_absolute():
        output_dir = FIGURE_OUTPUT_DIR / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_stem = f"Q_e_passage_ElNino_amp{AMP_val}_M{Q_M_val}_{DATAKEY}"

    def save_figure(fig, suffix):
        if suffix not in output_variants:
            plt.close(fig)
            return
        output_base = output_dir / f"{output_stem}_{suffix}"
        eps_path = output_base.with_suffix(".eps")
        png_path = output_base.with_suffix(".png")
        fig.savefig(eps_path, format="eps", bbox_inches="tight")
        fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"saved figure: {eps_path}")
        print(f"saved figure: {png_path}")


    datas = get_data(key=DATAKEY)
    amp_rank_val = amp_rank(AMP_val)
    q_result_path = (
        datas["all"]["dataset_dir"]
        / "SRC_py"
        / "Interim_data_Q"
        / f"Q_result_fixed_M{Q_M_val}_amp{amp_rank_val}.npz"
    )

    if not q_result_path.exists():
        raise FileNotFoundError(q_result_path)

    with np.load(q_result_path) as q_result:
        theta_a_grid = np.asarray(q_result["phi1"], dtype=float)
        theta_e_grid = np.asarray(q_result["phi2"], dtype=float)
        Q_rad_per_year = np.asarray(q_result["Q"], dtype=float)
        noise_cov = np.asarray(q_result["noise_cov"], dtype=float)
        Q = Q_rad_per_year / (2.0 * np.pi)
        omega_e = float(np.asarray(q_result["omega_freq_per_year"], dtype=float)[1])

    Q_a = Q[0]
    Q_e = Q[1]

    if theta_a_grid.shape != theta_e_grid.shape or Q_e.shape != theta_a_grid.shape:
        raise ValueError(
            f"shape mismatch: theta_a={theta_a_grid.shape}, "
            f"theta_e={theta_e_grid.shape}, Q_e={Q_e.shape}"
        )

    print({
        "DATAKEY": DATAKEY,
        "amp_rank": f"{amp_rank_val}%",
        "M": Q_M_val,
        "result_path": q_result_path,
        "grid_shape": Q_e.shape,
        "Q_e_min [yr^-1]": float(np.nanmin(Q_e)),
        "Q_e_max [yr^-1]": float(np.nanmax(Q_e)),
        "omega_e [yr^-1]": omega_e,
    })


    # ## Display $Q_e(\theta_a,\theta_e)$
    #
    # 等高線、カラーバー範囲、位相 bin の境界・ラベル、annual phase の月表示を `Plot_Qe.py` と同じ形式にします。

    # In[3]:


    TWOPI = 2.0 * np.pi
    PHASE_BIN_COUNT = 8
    MONTH_TICKS = np.linspace(0.0, TWOPI, 13)
    MONTH_TICKLABELS = ["J.", "F.", "M.", "A.", "M.", "J.", "J.", "A.", "S.", "O.", "N.", "D.", "J."]


    def add_phase_bin_labels(ax, phase_bin_count=PHASE_BIN_COUNT):
        width = TWOPI / phase_bin_count
        text_style = {"fontsize": 7, "ha": "left", "color": "0.25", "zorder": 5}

        for boundary in width * (np.arange(phase_bin_count) + 0.5):
            ax.axhline(boundary, color="0.35", linestyle="--", linewidth=0.55, zorder=4)

        for i in range(phase_bin_count):
            y = i * width
            if i == 0:
                ax.text(0.08, 0.08, "j=1", va="bottom", **text_style)
                ax.text(0.08, TWOPI - 0.08, "j=1", va="top", **text_style)
            else:
                ax.text(0.08, y, f"j={i + 1}", va="center", **text_style)


    def set_phase_axis(ax):
        phase_ticks = np.linspace(0.0, TWOPI, 5)
        phase_ticklabels = [r"$0$", r"$0.5\pi$", r"$\pi$", r"$1.5\pi$", r"$2\pi$"]
        ax.set_xlim(0.0, TWOPI)
        ax.set_ylim(0.0, TWOPI)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks(MONTH_TICKS)
        ax.set_xticklabels(MONTH_TICKLABELS)
        ax.tick_params(axis="x", labelsize=8, pad=2)
        ax.set_xlabel(r"$\theta_a$ (Month)")
        ax.set_yticks(phase_ticks)
        ax.set_yticklabels(phase_ticklabels)
        ax.set_ylabel(r"$\theta_e$")


    colorbar_half_width = 0.15
    levels = np.linspace(omega_e - colorbar_half_width, omega_e + colorbar_half_width, 17)

    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    im = ax.contourf(
        theta_a_grid, theta_e_grid, Q_e,
        levels=levels, cmap="seismic",
        vmin=levels[0], vmax=levels[-1], extend="both",
    )
    ax.contour(
        theta_a_grid, theta_e_grid, Q_e,
        levels=levels, colors="0.25", linewidths=0.15, alpha=1.0,
    )

    q_e_min = float(np.nanmin(Q_e))
    q_e_max = float(np.nanmax(Q_e))
    if q_e_min <= 0.0 <= q_e_max:
        ax.contour(
            theta_a_grid, theta_e_grid, Q_e,
            levels=[0.0], colors="0.15", linewidths=0.25, alpha=1.0,
        )

    add_phase_bin_labels(ax)
    set_phase_axis(ax)
    ax.set_title(r"$Q_e(\theta_a,\theta_e)$")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(im, cax=cax, label=r"[yr$^{-1}$]")
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    cbar.ax.axhline(omega_e, color="black", linewidth=1.2)

    fig.tight_layout()
    save_figure(fig, "Qe")


    # ## Passage time from the left edge of $j=8$ to the right edge of $j=2$
    #
    # `Plot_Qe.py` の phase bin は幅 $2\pi/8$ で、$j=1$ の中心が $\theta_e=0$ です。したがって、増加方向に進む今回の区間は
    #
    # $$\theta_{e,\mathrm{start}}=\frac{13\pi}{8} \quad (j=8\text{ の左端}),$$
    #
    # $$\theta_{e,\mathrm{end}}=2\pi+\frac{3\pi}{8} \quad (j=2\text{ の右端})$$
    #
    # です。各月の $\theta_a$ は初期値としてだけ指定し、その後は $Q_a$ と $Q_e$ の両方を使って時間発展させます。
    #
    # $$\frac{d\theta_a}{dt}=Q_a(\theta_a,\theta_e),\qquad \frac{d\theta_e}{dt}=Q_e(\theta_a,\theta_e).$$
    #
    # 初期条件は $(\theta_a(0),\theta_e(0))=(\text{各月の位相},13\pi/8)$、終了条件は非ラップの $\theta_e\geq2\pi+3\pi/8$ です。場を補間するときだけ両位相を $[0,2\pi)$ に戻します。`Plot_Qe_ODEOverlay.py` と同じEuler法、$dt=1/120$ yearを使い、最後のstep内の到達時刻を線形補間します。

    # In[4]:


    import pandas as pd

    MONTH_NAMES = ["Jan.", "Feb.", "Mar.", "Apr.", "May", "Jun.",
                   "Jul.", "Aug.", "Sep.", "Oct.", "Nov.", "Dec."]
    MONTH_PHASES = np.arange(12, dtype=float) * TWOPI / 12.0

    J8_LEFT_EDGE = 13.0 * np.pi / 8.0
    J2_RIGHT_EDGE = 3.0 * np.pi / 8.0
    TARGET_THETA_E = TWOPI + J2_RIGHT_EDGE
    SIM_DT_YEARS = 1.0 / 120.0
    SIM_MAX_YEARS = 10.0
    SIM_MAX_STEPS = int(round(SIM_MAX_YEARS / SIM_DT_YEARS))

    theta_a_axis = np.asarray(theta_a_grid[0, :], dtype=float)
    theta_e_axis = np.asarray(theta_e_grid[:, 0], dtype=float)

    f_a = RegularGridInterpolator(
        (theta_a_axis, theta_e_axis),
        Q_rad_per_year[0].T,
        method="linear",
        bounds_error=True,
    )
    f_e = RegularGridInterpolator(
        (theta_a_axis, theta_e_axis),
        Q_rad_per_year[1].T,
        method="linear",
        bounds_error=True,
    )

    def phase_drift(theta):
        """Return [Q_a, Q_e] in rad/year at one unwrapped phase point."""
        theta_a, theta_e = np.mod(np.asarray(theta, dtype=float), TWOPI)
        point = (theta_a, theta_e)
        return np.array([
            f_a(point).item(),
            f_e(point).item(),
        ])


    def simulate_ode_first_passage(theta_a_initial):
        """Integrate both Q_a and Q_e and return first-passage time [year]."""
        theta = np.array([theta_a_initial, J8_LEFT_EDGE], dtype=float)

        for step in range(SIM_MAX_STEPS):
            previous_theta = theta.copy()
            theta = theta + phase_drift(theta) * SIM_DT_YEARS

            if theta[1] >= TARGET_THETA_E:
                fraction = (
                    (TARGET_THETA_E - previous_theta[1])
                    / (theta[1] - previous_theta[1])
                )
                fraction = float(np.clip(fraction, 0.0, 1.0))
                return (step + fraction) * SIM_DT_YEARS

        raise RuntimeError(f"ODE did not reach within {SIM_MAX_YEARS:g} years.")


    passage_times_years = []
    for theta_a in MONTH_PHASES:
        passage_times_years.append(simulate_ode_first_passage(theta_a))

    passage_times_years = np.asarray(passage_times_years)
    passage_times_months = 12.0 * passage_times_years

    passage_time_table = pd.DataFrame({
        "month": MONTH_NAMES,
        "initial theta_a [rad]": MONTH_PHASES,
        "passage time [year]": passage_times_years,
        "passage time [month]": passage_times_months,
    })

    print(passage_time_table.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


    # ## ODE passage time for each initial month

    # In[5]:


    fig, ax = plt.subplots(figsize=(3.4, 2.0))
    month_index = np.arange(12)

    ax.plot(
        month_index, passage_times_months,
        color="royalblue", linewidth=1.2, marker="o", markersize=4,
    )
    ax.set_xticks(month_index)
    ax.set_xticklabels([month[0] + "." for month in MONTH_NAMES])
    ax.set_xlabel("Month")
    ax.set_ylabel("passage time [month]")
    ax.set_title("El Niño")
    ax.grid(axis="y", linestyle="-.", color="0.82", alpha=1.0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))

    # for x, y in zip(month_index, passage_times_months):
    #     ax.annotate(f"{y:.2f}", (x, y), xytext=(0, 6), textcoords="offset points",
    #                 ha="center", va="bottom", fontsize=7)

    fig.tight_layout()
    save_figure(fig, "ODE")


    # ## SDE passage-time statistics (5000 trials per month)
    #
    # ODE版と同じく $Q_a$ と $Q_e$ の両方を使い、2次元SDE
    #
    # $$d\boldsymbol{\theta}=\boldsymbol{Q}(\boldsymbol{\theta})dt+L\,d\boldsymbol{W},\qquad LL^{\mathsf{T}}=E$$
    #
    # をEuler--Maruyama法で進めます。$E$ は推定結果に保存された2次元ノイズ共分散行列 `noise_cov` で、Cholesky因子 $L$ により $\theta_a$ と $\theta_e$ の相関ノイズを生成します。初期条件・終了条件・$dt=1/120$ yearはODE版と同じで、各初期月5000試行の初通過時間を計算します。

    # In[6]:


    SDE_SEED = 23

    noise_cov = 0.5 * (noise_cov + noise_cov.T)
    noise_cholesky = np.linalg.cholesky(noise_cov + 1e-14 * np.eye(2))


    def simulate_sde_first_passage_times(theta_a_initial, trial_count, seed):
        """Integrate both phases with covariance E; return times [month]."""
        rng = np.random.default_rng(seed)
        first_passage_months = np.full(int(trial_count), np.nan, dtype=float)
        noise_scale = np.sqrt(SIM_DT_YEARS)

        for trial in range(int(trial_count)):
            theta = np.array([theta_a_initial, J8_LEFT_EDGE], dtype=float)

            for step in range(SIM_MAX_STEPS):
                previous_theta = theta.copy()
                correlated_noise = noise_scale * (rng.standard_normal(2) @ noise_cholesky.T)
                theta = theta + phase_drift(theta) * SIM_DT_YEARS + correlated_noise

                if theta[1] >= TARGET_THETA_E:
                    fraction = (
                        (TARGET_THETA_E - previous_theta[1])
                        / (theta[1] - previous_theta[1])
                    )
                    fraction = float(np.clip(fraction, 0.0, 1.0))
                    first_passage_months[trial] = (step + fraction) * SIM_DT_YEARS * 12.0
                    break

        return first_passage_months


    sde_passage_times_months = np.empty((12, SDE_TRIAL_COUNT), dtype=float)
    for month_index_value, theta_a in enumerate(MONTH_PHASES):
        sde_passage_times_months[month_index_value] = simulate_sde_first_passage_times(
            theta_a_initial=theta_a, trial_count=SDE_TRIAL_COUNT,
            seed=SDE_SEED + month_index_value,
        )

    sde_reached_count = np.sum(np.isfinite(sde_passage_times_months), axis=1)
    if np.any(sde_reached_count != SDE_TRIAL_COUNT):
        raise RuntimeError(
            f"Some trials did not reach within {SIM_MAX_YEARS:g} years: "
            f"{sde_reached_count.tolist()}"
        )

    sde_mean_months = np.nanmean(sde_passage_times_months, axis=1)
    sde_std_months = np.nanstd(sde_passage_times_months, axis=1, ddof=1)
    sde_median_months = np.nanmedian(sde_passage_times_months, axis=1)
    sde_p025_months = np.nanpercentile(sde_passage_times_months, 2.5, axis=1)
    sde_p975_months = np.nanpercentile(sde_passage_times_months, 97.5, axis=1)

    sde_statistics_table = pd.DataFrame({
        "month": MONTH_NAMES,
        "trials": SDE_TRIAL_COUNT,
        "reached": sde_reached_count,
        "mean [month]": sde_mean_months,
        "std [month]": sde_std_months,
        "median [month]": sde_median_months,
        "p2.5 [month]": sde_p025_months,
        "p97.5 [month]": sde_p975_months,
    })

    print({
        "trials per month": SDE_TRIAL_COUNT,
        "dt [year]": SIM_DT_YEARS,
        "noise covariance E [rad^2/year]": noise_cov,
        "seed base": SDE_SEED,
    })
    print(sde_statistics_table.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


    # ## SDE passage time for each initial month

    # In[7]:


    fig, ax = plt.subplots(figsize=(3.4, 2.0))
    month_index = np.arange(12)

    ax.fill_between(
        month_index, sde_p025_months, sde_p975_months,
        color="#DCE5FA", alpha=1.0, linewidth=0.0,
    )
    ax.plot(
        month_index, sde_mean_months,
        color="royalblue", linewidth=1.2, marker="o", markersize=4,
    )
    ax.set_xticks(month_index)
    ax.set_xticklabels([month[0] + "." for month in MONTH_NAMES])
    ax.set_xlabel("Month")
    ax.set_ylabel("passage time [month]")
    ax.set_title("El Niño")
    ax.grid(axis="y", linestyle="-.", color="0.82", alpha=1.0)
    if DATAKEY == "ERSSTv5_2026" and amp_rank_val == 30 and Q_M_val == 3:
        ymin, ymax = ax.get_ylim()
        ax.set_yticks([tick for tick in (15, 17, 19, 21, 23) if ymin <= tick <= ymax])
        ax.set_ylim(ymin, ymax)
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    fig.tight_layout()
    save_figure(fig, "SDE")





if __name__ == "__main__":
    trials = int(os.environ.get("SDE_TRIAL_COUNT", "5000"))
    if trials < 1:
        raise ValueError("SDE_TRIAL_COUNT must be positive.")
    main(SDE_TRIAL_COUNT=trials, savedir=FIGURE_OUTPUT_DIR)
