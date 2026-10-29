
from plot_config import FIGURE_OUTPUT_DIR
#!/usr/bin/env python3
# coding: utf-8

# %% 




import datetime
from pathlib import Path
from typing import get_args
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Common import DatasetKey, get_data

plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams['font.family'] = 'Times New Roman' # font familyの設定
plt.rcParams['font.size'] = 10 # 全体のフォントサイズが変更されます。
plt.rcParams['xtick.labelsize'] = 10 # 軸だけ変更されます.
plt.rcParams['ytick.labelsize'] = 10 # 軸だけ変更されます.
plt.rcParams['figure.facecolor']= 'white'




DATAKEY = "ERSSTv5_2026"
AMP_val = int(30) #* 振幅フィルター30%の条件で推定した結果を読み込みます
Q_M_val = int(3) #* M=3の条件で推定した結果を読み込みます

TWOPI = 2.0 * np.pi
PHASE_BIN_COUNT = 8
MONTH_TICKS = np.linspace(0.0, TWOPI, 13)
MONTH_TICKLABELS = ["J.", "F.", "M.", "A.", "M.", "J.", "J.", "A.", "S.", "O.", "N.", "D.", "J."]


def add_phase_bin_labels(ax, phase_bin_count=PHASE_BIN_COUNT):
    width = TWOPI / phase_bin_count
    x = 0.08
    text_style = {
        "fontsize": 7,
        "ha": "left",
        "color": "0.25",
        "zorder": 5,
    }

    for boundary in width * (np.arange(phase_bin_count) + 0.5):
        ax.axhline(
            boundary,
            color="0.35",
            linestyle="--",
            linewidth=0.55,
            zorder=4,
        )

    for i in range(phase_bin_count):
        y = i * width
        if i == 0:
            ax.text(x, 0.08, "j=1", va="bottom", **text_style)
            ax.text(x, TWOPI - 0.08, "j=1", va="top", **text_style)
        else:
            ax.text(x, y, f"j={i + 1}", va="center", **text_style)


def set_phase_axis(ax, theta_a_month_tick_key=False):
    phase_ticks = np.linspace(0.0, TWOPI, 5)
    # 描画座標はradのまま、目盛りと軸ラベルをcyclesで表示する。
    phase_ticklabels = ["0", "0.25", "0.5", "0.75", "1"]
    ax.set_xlim(0.0, TWOPI)
    ax.set_ylim(0.0, TWOPI)
    ax.set_aspect("equal", adjustable="box")

    if theta_a_month_tick_key:
        ax.set_xticks(MONTH_TICKS)
        ax.set_xticklabels(MONTH_TICKLABELS)
        ax.tick_params(axis="x", labelsize=8, pad=2)
        ax.set_xlabel(r"$\frac{\theta_a}{2\pi}$ (Month)")
    else:
        ax.set_xticks(phase_ticks)
        ax.set_xticklabels(phase_ticklabels)
        ax.set_xlabel(r"$\frac{\theta_a}{2\pi}$")

    ax.set_yticks(phase_ticks)
    ax.set_yticklabels(phase_ticklabels)
    ax.set_ylabel(r"$\frac{\theta_e}{2\pi}$")


def main(
        DATAKEY=DATAKEY,
        AMP_val=AMP_val,
        Q_M_val=Q_M_val,
        savedir = "./",
        png_output = True,
        eps_output = False,
        p_label_key = False,
        theta_a_month_tick_key = False,
        amplitude_threshold_title = False,
):

    #* datas key: "annual_enso", "all"
    #* datas["annual_enso"] key: "date", "nino3.4", "Z_annual", "Z_enso",
    #*     "Z_annual_conj", "Z_enso_conj", "raw_Z_annual", "raw_Z_enso",
    #*     "raw_Z_annual_conj", "raw_Z_enso_conj", "raw_Z_annual_left",
    #*     "raw_Z_enso_left", "raw_Z_annual_conj_left", "raw_Z_enso_conj_left",
    #*     "lambda_annual", "lambda_enso", "lambda_annual_conj",
    #*     "lambda_enso_conj", "gamma_annual", "gamma_enso",
    #*     "gamma_annual_conj", "gamma_enso_conj", "mode_info", "phase_origin"
    #* datas["all"] key: "key", "dataset_dir", "myfolder", "date", "Z",
    #*     "Z_left", "nino34", "Lambda", "mu", "mode_info", "phase_origin"
    #* -----
    datas = get_data(key=DATAKEY)

   
    amp_rank = int(round(AMP_val * 100.0 if AMP_val <= 1.0 else AMP_val))
    q_result_path = (
        datas["all"]["dataset_dir"]
        / "SRC_py"
        / "Interim_data_Q"
        / f"Q_result_fixed_M{Q_M_val}_amp{amp_rank}.npz"
    )

    with np.load(q_result_path) as q_result:
        xx = np.asarray(q_result["phi1"], dtype=float)
        yy = np.asarray(q_result["phi2"], dtype=float)
        #* Q は rad/year で保存されているので、ここで 2π で割って yr^{-1} に変換する。
        Q = np.asarray(q_result["Q"], dtype=float) / (2.0 * np.pi)
        omega_e = float(np.asarray(q_result["omega_freq_per_year"], dtype=float)[1])

    Q_a = Q[0]
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

    if p_label_key:
        add_phase_bin_labels(ax)

    set_phase_axis(ax, theta_a_month_tick_key=theta_a_month_tick_key)
    if amplitude_threshold_title:
        ax.set_title(
            rf"$\tilde{{R}}_{{\mathrm{{e}}}} = \tilde{{R}}_{{{amp_rank}\,\mathrm{{th}}}}$"
        )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(im, cax=cax, label=r"$\frac{Q_e(\theta_a,\theta_e)}{2\pi}$ [cycle / year]")
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    cbar.ax.axhline(omega_e, color="black", linewidth=1.2)

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
    if eps_output:
        fig.savefig(eps_path, format="eps", bbox_inches="tight")

    plt.close(fig)

    if png_output:
        print(f"saved figure: {png_path}")
    if eps_output:
        print(f"saved figure: {eps_path}")




# %% 


# %% 

if __name__ == "__main__":
    main(AMP_val=30, savedir=f"Q_e_amp30_M3_{DATAKEY}",
         png_output=True, eps_output=True, p_label_key=True, theta_a_month_tick_key=True)
    for amp in (10, 20, 30, 40, 50):
        main(AMP_val=amp, savedir=f"素材/Q_e/Q_e_amp{amp}_M3_{DATAKEY}_threshold_title",
             png_output=True, eps_output=True, p_label_key=True, theta_a_month_tick_key=True,
             amplitude_threshold_title=True)
