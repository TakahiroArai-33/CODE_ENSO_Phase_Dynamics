#!/usr/bin/env python3
# coding: utf-8

"""Corrected 2D Q(phi_a, phi_e) estimator.

修正点:
- ``Module_Bayes_Routine.estimate_q`` は ``修正版`` 内にコピーした実装を使う。
- 入力位相だけを修正版に差し替える。
- ENSO 位相は ``PHI_e_SHIFT`` で Nino3.4 ピーク平均位相を 0 に合わせる。
- 出力先は ``SRC_py/Interim_data_Q``。
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "Interim_data_Q"
TWOPI = 2.0 * 3.141592653589793


def parse_args() -> argparse.Namespace:
    """Parse command-line options for the corrected two-dimensional Q estimate."""
    parser = argparse.ArgumentParser(description="Estimate Q(phi_a, phi_e) with corrected phase convention.")
    parser.add_argument("--env", default=None, help="conda environment name, e.g. virtual_env")
    parser.add_argument("--M-list", nargs="+", type=int, default=[3, 5])
    parser.add_argument("--amp-ranks", nargs="+", type=float, default=list(range(0, 81, 10)))
    parser.add_argument("--lambda-0", type=float, default=1e-6)
    parser.add_argument("--num-iterations", type=int, default=100)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--grid-num", type=int, default=2**9)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def rerun_in_conda_env(args: argparse.Namespace) -> None:
    """Re-execute this script inside the requested conda environment if needed."""
    if not args.env:
        return
    if os.environ.get("CONDA_DEFAULT_ENV") == args.env:
        return
    if os.environ.get("Q2D_FIXED_ENV_REEXEC") == args.env:
        return
    env = os.environ.copy()
    env["Q2D_FIXED_ENV_REEXEC"] = args.env
    cmd = ["conda", "run", "-n", args.env, "python", str(Path(__file__).resolve()), *sys.argv[1:]]
    raise SystemExit(subprocess.run(cmd, env=env).returncode)


def evaluate_phase_equation(result, num: int):
    """
    推定済みの位相方程式を、規則的な ``phi_a``/``phi_e`` 格子上で評価する。

    推定器は任意の位相点で Q を返せるが、図にしたり保存したりするには
    2 次元格子上の値が必要になる。この関数はそのために、0 から 2*pi まで
    の格子を作り、``Q_a`` と ``Q_e`` を全格子点で計算する。
    """
    import numpy as np

    phi1_grid = np.linspace(0.0, TWOPI, num=num, endpoint=True)
    phi2_grid = np.linspace(0.0, TWOPI, num=num, endpoint=True)
    phi1_mesh, phi2_mesh = np.meshgrid(phi1_grid, phi2_grid, indexing="xy")
    q1 = result.phi1.predict(phi1_mesh.ravel(), phi2_mesh.ravel()).reshape(phi1_mesh.shape)
    q2 = result.phi2.predict(phi1_mesh.ravel(), phi2_mesh.ravel()).reshape(phi2_mesh.shape)
    return phi1_grid, phi2_grid, phi1_mesh, phi2_mesh, q1, q2


def save_result(args, M, amp_rank, result, phase_origin):
    """
    Q 推定結果を `.npz` ファイルに保存する。

    保存するのは、格子上の ``Q_a``/``Q_e`` だけではない。推定に実際に
    使った中点位相、位相速度、amp filter の mask、ノイズ共分散、平均
    drift、位相基準の情報も一緒に保存する。後で図を作るとき、推定時と
    同じ座標・同じデータを確認できるようにするためである。
    """
    import numpy as np

    phi1_grid, phi2_grid, phi1_mesh, phi2_mesh, q1, q2 = evaluate_phase_equation(result, args.grid_num)
    Q = np.stack((q1, q2), axis=0)
    omega = np.array([np.real(result.phi1.coef[0]), np.real(result.phi2.coef[0])], dtype=float)
    mode_indices = phase_origin.get("mode_indices", {})
    output = args.output_dir / f"Q_result_fixed_M{M}_amp{int(amp_rank) if float(amp_rank).is_integer() else amp_rank:g}.npz"
    np.savez_compressed(
        output,
        phi1=phi1_mesh,
        phi2=phi2_mesh,
        Q=Q,
        phi1_grid=phi1_grid,
        phi2_grid=phi2_grid,
        data_phi1=result.data.phi1_ast,
        data_phi2=result.data.phi2_ast,
        data_dot_phi1=result.data.dot_phi1,
        data_dot_phi2=result.data.dot_phi2,
        keep=result.data.keep,
        noise_cov=result.noise_cov,
        omega=omega,
        omega_freq_per_year=omega / TWOPI,
        M=np.array(M, dtype=int),
        amp_rank=np.array(amp_rank, dtype=float),
        lambda_0=np.array(args.lambda_0, dtype=float),
        phi_a_shift_rad=np.array(phase_origin["phi_a_shift_rad"], dtype=float),
        phi_e_shift_rad=np.array(phase_origin["phi_e_shift_rad"], dtype=float),
        annual_mode_index=np.array(mode_indices.get("annual", -1), dtype=int),
        enso_mode_index=np.array(mode_indices.get("enso", -1), dtype=int),
        enso_peak_mask=np.asarray(phase_origin["enso_peak_mask"], dtype=bool),
        enso_peak_mean_phase_rad=np.array(phase_origin["enso_peak_mean_phase_rad"], dtype=float),
        phase_convention=np.array(phase_origin["phase_convention"]),
    )
    print(f"saved: {output}")
    return output


def main() -> int:
    """
    指定された M と amp の組み合わせすべてについて Q 推定を実行する。

    まず ``common_fixed.get_data()`` で修正版の複素時系列を取り出し、
    ``get_Rangle()`` で annual/ENSO の位相と振幅を作る。その後、
    ``estimate_q`` に渡す。
    """
    args = parse_args()
    rerun_in_conda_env(args)

    import numpy as np

    import common_fixed as common

    from Module_Bayes_Routine import estimate_q

    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = common.get_data()
    phi_a, _amp_a = common.get_Rangle(data["annual"])
    phi_e, amp_e = common.get_Rangle(data["enso"])
    phase_origin = data["phase_origin"]
    dt_year = 1.0 / 12.0

    outputs = []
    for M in args.M_list:
        for amp_rank in args.amp_ranks:
            print(f"=== Q fixed: M={M}, amp_rank={amp_rank:g} ===")
            result = estimate_q(
                phi1=phi_a,
                phi2=phi_e,
                dt=dt_year,
                amp=amp_e,
                amp_rank=amp_rank,
                M1=M,
                M2=M,
                lambda_0=args.lambda_0,
                num_iterations=args.num_iterations,
                rtol=args.rtol,
            )
            output = save_result(args, M, amp_rank, result, phase_origin)
            outputs.append(output)
            omega_freq = np.real([result.phi1.coef[0], result.phi2.coef[0]]) / TWOPI
            print(f"data points: {result.data.phi1_ast.size}, omega [freq/year]: {omega_freq}")
    print(f"done: {len(outputs)} files saved in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
