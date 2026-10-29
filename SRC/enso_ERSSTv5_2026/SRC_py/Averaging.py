#!/usr/bin/env python3
# coding: utf-8

"""Average estimated Q_a(phi_a, phi_e) and Q_e(phi_a, phi_e).

このスクリプトは、推定済みの ``Interim_data_Q/Q_result_fixed_M*_amp*.npz``
を読み込み、指定した m:n ごとに

    y1: Q_a を同期軌道方向に平均したもの
    y2: Q_e を同期軌道方向に平均したもの
    y3: d/dt(phi_a - (n/m) phi_e)

を ``Interim_Averaging`` に保存する。

``y1``, ``y2``, ``y3`` は DirectAveraging の結果と同じく yr^-1 単位で保存する。
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from tempfile import gettempdir

os.environ.setdefault("MPLCONFIGDIR", str(Path(gettempdir()) / "ersstv5-matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(gettempdir()) / "ersstv5-cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
from scipy.interpolate import RegularGridInterpolator


SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / "Interim_data_Q"
INTERIM_DIR = SCRIPT_DIR / "Interim_Averaging"

TWOPI = 2.0 * np.pi
DATA_NAME_PATTERN = re.compile(r"^Q_result(?:_fixed)?_M(?P<M>\d+)_amp(?P<amp>\d+)\.npz$")
MN_PATTERN = re.compile(r"^(?P<m>\d+):(?P<n>\d+)$")
ALLOWED_MN = ("1:2", "1:3", "1:4", "1:5", "2:5")


def parse_mn(value: str) -> tuple[int, int]:
    match = MN_PATTERN.fullmatch(value)
    if match is None:
        raise argparse.ArgumentTypeError("m:n の形式で指定してください。例: 1:2")
    ratio = (int(match["m"]), int(match["n"]))
    if f"{ratio[0]}:{ratio[1]}" not in ALLOWED_MN:
        raise argparse.ArgumentTypeError(f"m:n は {', '.join(ALLOWED_MN)} のいずれかを指定してください。")
    return ratio


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Average Q_result_fixed_M*_amp*.npz.")
    parser.add_argument(
        "data_name",
        type=Path,
        nargs="?",
        help="単発実行する Q_result npz。省略時は input-dir 内の全ファイルを処理する。",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=INPUT_DIR,
        help="Q_result_fixed_M*_amp*.npz を読むディレクトリ。デフォルト: Interim_data_Q",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=INTERIM_DIR,
        help="平均化結果を保存するディレクトリ。デフォルト: Interim_Averaging",
    )
    parser.add_argument(
        "--ratio",
        "--ratios",
        dest="ratios",
        type=parse_mn,
        nargs="+",
        default=[parse_mn(value) for value in ALLOWED_MN],
        metavar="m:n",
        help="同期状態における振動数比 m:n。省略時は全候補を実行する。",
    )
    parser.add_argument("--orders", type=int, nargs="+", default=None, metavar="M")
    parser.add_argument("--amps", type=int, nargs="+", default=None, metavar="AMP")
    return parser.parse_args()


def parse_data_name(data_name: Path) -> tuple[int, int]:
    match = DATA_NAME_PATTERN.fullmatch(data_name.name)
    if match is None:
        raise ValueError(f"invalid data name: {data_name.name}")
    return int(match["M"]), int(match["amp"])


def ratio_label(ratio: tuple[int, int]) -> str:
    return f"{ratio[0]:d}:{ratio[1]:d}"


def resolve_input_path(data_name: Path, input_dir: Path, M: int) -> Path:
    if data_name.is_absolute() and data_name.exists():
        return data_name

    candidates = [
        input_dir / data_name,
        input_dir / f"M{M:d}" / data_name,
        data_name,
    ]
    for path in candidates:
        if path.exists():
            return path

    tried = "\n  ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Q result file not found. Tried:\n  {tried}")


def collect_data_names(input_dir: Path, orders: set[int] | None, amps: set[int] | None) -> list[Path]:
    data_names = []
    for path in input_dir.glob("Q_result*_M*_amp*.npz"):
        match = DATA_NAME_PATTERN.fullmatch(path.name)
        if match is None:
            continue
        order = int(match["M"])
        amp = int(match["amp"])
        if orders is not None and order not in orders:
            continue
        if amps is not None and amp not in amps:
            continue
        data_names.append(path)
    return sorted(data_names, key=lambda name: parse_data_name(name))


def make_phase_interpolator(
    phi_a_grid: np.ndarray,
    phi_e_grid: np.ndarray,
    values: np.ndarray,
):
    # values[i_e, j_a] = Q(phi_a_grid[j_a], phi_e_grid[i_e])
    interpolator = RegularGridInterpolator((phi_a_grid, phi_e_grid), values.T)

    def evaluate(phi_a, phi_e):
        points = np.column_stack((np.ravel(np.mod(phi_a, TWOPI)), np.ravel(np.mod(phi_e, TWOPI))))
        return interpolator(points).reshape(np.shape(phi_a))

    return evaluate


def noise_cov_values(data: np.lib.npyio.NpzFile, m_e: int, m_a: int) -> tuple[float, float, float, float]:
    if "noise_cov" not in data.files:
        return 0.0, 0.0, 0.0, 0.0

    noise_cov = np.asarray(data["noise_cov"], dtype=float)
    ratio = m_a / m_e
    Eaa = noise_cov[0, 0] / TWOPI**2
    Eee = noise_cov[1, 1] / TWOPI**2
    Eae = noise_cov[0, 1] / TWOPI**2
    E_delta = Eaa - 2.0 * ratio * Eae + ratio**2 * Eee
    return Eaa, Eee, Eae, E_delta


def average_q(data: np.lib.npyio.NpzFile, ratio_mn: tuple[int, int]) -> dict[str, np.ndarray | float | int]:
    m_e, m_a = ratio_mn
    ratio = m_a / m_e

    Q = np.asarray(data["Q"], dtype=float)
    Qa = Q[0]
    Qe = Q[1]
    phi_a_grid = np.asarray(data["phi1_grid"], dtype=float)
    phi_e_grid = np.asarray(data["phi2_grid"], dtype=float)

    func_Qa = make_phase_interpolator(phi_a_grid, phi_e_grid, Qa)
    func_Qe = make_phase_interpolator(phi_a_grid, phi_e_grid, Qe)

    psi_grid = np.linspace(0.0, TWOPI, phi_a_grid.size, endpoint=True)
    theta_num = (phi_e_grid.size - 1) * m_e + 1
    theta_grid = np.linspace(0.0, TWOPI * m_e, theta_num, endpoint=True)

    gamma_a_from_e = np.zeros_like(psi_grid)
    gamma_a_to_e = np.zeros_like(psi_grid)

    for i, psi in enumerate(psi_grid):
        phi_a_for_average = psi + ratio * theta_grid
        phi_e_for_average = theta_grid
        gamma_a_from_e[i] = np.trapz(func_Qa(phi_a_for_average, phi_e_for_average), theta_grid) / (TWOPI * m_e)
        gamma_a_to_e[i] = np.trapz(func_Qe(phi_a_for_average, phi_e_for_average), theta_grid) / (TWOPI * m_e)

    # DirectAveraging の bx と同じ向きにするため、Q_e 側は -psi の関数として保存する。
    gamma_a_to_e_rev = np.concatenate(([gamma_a_to_e[0]], gamma_a_to_e[-2:0:-1], [gamma_a_to_e[0]]))
    gamma_delta = gamma_a_from_e - ratio * gamma_a_to_e

    Eaa, Eee, Eae, E_delta = noise_cov_values(data, m_e=m_e, m_a=m_a)
    return {
        "x": psi_grid,
        "y1": gamma_a_from_e / TWOPI,
        "y2": gamma_a_to_e_rev / TWOPI,
        "y3": gamma_delta / TWOPI,
        "y1_err": np.zeros_like(psi_grid),
        "y2_err": np.zeros_like(psi_grid),
        "y3_err": np.zeros_like(psi_grid),
        "Gamma_a_from_e": gamma_a_from_e,
        "Gamma_a_to_e": gamma_a_to_e,
        "Gamma_a_to_e_rev": gamma_a_to_e_rev,
        "Gamma_delta": gamma_delta,
        "Eaa": Eaa,
        "Eee": Eee,
        "Eae": Eae,
        "E_delta": E_delta,
        "m_e": m_e,
        "m_a": m_a,
        "ratio": np.array([m_e, m_a], dtype=int),
    }


def save_averaged_result(data_path: Path, output_dir: Path, ratio_mn: tuple[int, int]) -> Path:
    M, amp = parse_data_name(data_path)
    with np.load(data_path) as data:
        result = average_q(data, ratio_mn)
        result["M"] = np.array(M, dtype=int)
        result["amp_rank"] = np.array(amp, dtype=float)
        if "omega" in data.files:
            result["omega"] = data["omega"]
            result["omega_freq_per_year"] = data["omega"] / TWOPI
        elif "omega_freq_per_year" in data.files:
            result["omega_freq_per_year"] = data["omega_freq_per_year"]

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{data_path.stem}_ratio{ratio_label(ratio_mn)}_averaged_result.npz"
    np.savez_compressed(output_path, **result)
    print(f"saved averaged result: {output_path}")
    return output_path


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    orders = set(args.orders) if args.orders is not None else None
    amps = set(args.amps) if args.amps is not None else None

    if args.data_name is None:
        data_names = collect_data_names(input_dir, orders=orders, amps=amps)
    else:
        data_names = [args.data_name]

    if not data_names:
        print(f"no Q_result files found in {input_dir}")
        return 1

    saved = []
    for data_name in data_names:
        M, amp = parse_data_name(data_name)
        data_path = resolve_input_path(data_name, input_dir=input_dir, M=M)
        for ratio_mn in args.ratios:
            print(f"M={M:d}, amp={amp:d}, ratio={ratio_label(ratio_mn)}")
            saved.append(save_averaged_result(data_path, output_dir=output_dir, ratio_mn=ratio_mn))

    # Langevin シミュレーションとヒストグラム用の保存は、今回の用途では使わないので実行しない。
    # 必要になった場合は古い履歴から simulate_langevin 周りを戻す。
    print(f"done: {len(saved)} files saved in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
