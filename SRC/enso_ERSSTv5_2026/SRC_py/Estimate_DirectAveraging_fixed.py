#!/usr/bin/env python3
# coding: utf-8

"""Corrected DirectAveraging estimator.

修正点:
- 推定ロジックは ``SRC_py/DirectAvaraging`` にコピーした実装を使う。
- 入力位相だけを修正版に差し替える。
- ENSO 位相は ``PHI_e_SHIFT`` で Nino3.4 ピーク平均位相を 0 に合わせる。
- 出力先は ``SRC_py/Interim_data_DirectAveraging``。
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
LOCAL_DIRECT_SCRIPT = SCRIPT_DIR / "DirectAvaraging" / "Module_DirectAveraging_Core.py"
OUTPUT_DIR = SCRIPT_DIR / "Interim_data_DirectAveraging"
DEFAULT_RATIOS = ((1, 2), (1, 3), (1, 4), (1, 5), (2, 5))
DEFAULT_AMPS = tuple(range(0, 81, 10))
DEFAULT_ORDERS = (2, 3, None)


def parse_ratio(value: str) -> tuple[int, int]:
    """Parse a ratio string such as ``1:3`` into integer phase coefficients."""
    m_text, n_text = value.split(":", maxsplit=1)
    return int(m_text), int(n_text)


def parse_order(value: str) -> int | None:
    """Parse a fixed Fourier order or the ``auto`` model-selection sentinel."""
    if value.lower() in {"auto", "none", "select"}:
        return None
    return int(value)


def parse_args() -> argparse.Namespace:
    """Parse command-line options for corrected DirectAveraging estimation."""
    parser = argparse.ArgumentParser(description="DirectAveraging with corrected phase convention.")
    parser.add_argument("--env", default=None, help="conda environment name, e.g. virtual_env")
    parser.add_argument("--amp", type=float, default=None, help="single-run amplitude percentile")
    parser.add_argument(
        "--amps",
        "--amp-ranks",
        dest="amps",
        nargs="+",
        type=float,
        default=DEFAULT_AMPS,
        help="batch amplitude percentile ranks",
    )
    parser.add_argument("--ratio", type=parse_ratio, default=None, help="single-run ratio, e.g. 1:2")
    parser.add_argument("--ratios", nargs="+", type=parse_ratio, default=DEFAULT_RATIOS)
    parser.add_argument("--M", "--order", dest="order", type=parse_order, default=None)
    parser.add_argument("--orders", nargs="+", type=parse_order, default=DEFAULT_ORDERS)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--log10-lambda-min", type=int, default=0)
    parser.add_argument("--log10-lambda-max", type=int, default=8)
    parser.add_argument("--alpha-0", type=float, default=2.0)
    parser.add_argument("--beta-0", type=float, default=2.0)
    parser.add_argument("--langevin-steps", type=int, default=5000 * 12)
    parser.add_argument("--langevin-step-factor", type=int, default=None)
    parser.add_argument("--langevin-seed", type=int, default=11)
    return parser.parse_args()


def rerun_in_conda_env(args: argparse.Namespace) -> None:
    """Re-execute this script inside the requested conda environment if needed."""
    if not args.env:
        return
    if os.environ.get("CONDA_DEFAULT_ENV") == args.env:
        return
    if os.environ.get("DIRECT_AVERAGING_FIXED_ENV_REEXEC") == args.env:
        return
    env = os.environ.copy()
    env["DIRECT_AVERAGING_FIXED_ENV_REEXEC"] = args.env
    cmd = ["conda", "run", "-n", args.env, "python", str(Path(__file__).resolve()), *sys.argv[1:]]
    raise SystemExit(subprocess.run(cmd, env=env).returncode)


def load_direct_module():
    """Load the local DirectAvaraging estimation core and its runtime imports."""
    spec = importlib.util.spec_from_file_location("local_direct_averaging", LOCAL_DIRECT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {LOCAL_DIRECT_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.load_runtime_dependencies()
    return module


def make_direct_args(base_args: argparse.Namespace, amp: float, ratio: tuple[int, int], order: int | None) -> argparse.Namespace:
    """Build the argument namespace expected by the DirectAvaraging core module."""
    return argparse.Namespace(
        env=None,
        amp=amp,
        amps=base_args.amps,
        order=order,
        orders=base_args.orders,
        ratio=ratio,
        ratios=base_args.ratios,
        log10_lambda_min=base_args.log10_lambda_min,
        log10_lambda_max=base_args.log10_lambda_max,
        alpha_0=base_args.alpha_0,
        beta_0=base_args.beta_0,
        output_dir=base_args.output_dir,
        output=base_args.output,
        langevin_steps=base_args.langevin_steps,
        langevin_step_factor=base_args.langevin_step_factor,
        langevin_seed=base_args.langevin_seed,
    )


def estimate_one(direct, common, base_args, phi_a, phi_e, amp_e, phase_origin, amp, ratio, order) -> Path:
    """
    amp, ratio, M の 1 組について DirectAveraging 推定を実行する。

    修正版で作った ``phi_a``/``phi_e`` と振幅を、DirectAveraging の推定
    コアに渡す。振幅 filter 後に連続して残っている区間だけを使い、
    推定結果とその条件を `.npz` に保存する。
    """
    if not (0.0 <= float(amp) <= 100.0):
        raise ValueError("amp must be between 0 and 100.")
    run_args = make_direct_args(base_args, amp=amp, ratio=ratio, order=order)
    data = direct.make_estimation_data(
        phi1=phi_a,
        phi2=phi_e,
        amp=amp_e,
        amp_rank=amp,
    )
    result, phase_origin_out, data = direct.estimate_direct(
        run_args,
        phase_origin,
        data,
        initial_phase=common.np.array([phi_a[0], phi_e[0]], dtype=float),
        data_length=phi_a.size,
    )
    output = direct.save_result(
        run_args,
        result,
        phase_origin_out,
        data,
    )
    print(f"saved: {output} selected_M=({result['M1']}, {result['M2']})")
    return output


def main() -> int:
    """
    DirectAveraging 推定を単発または batch で実行する。

    ここでも ``common_fixed.get_data()`` で複素時系列を取り出し、
    ``get_Rangle()`` で annual/ENSO の位相と振幅を作る。Q 推定や
    phase composite と同じ修正版の位相基準になる。
    """
    args = parse_args()
    rerun_in_conda_env(args)

    import common_fixed as common

    direct = load_direct_module()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = common.get_data()
    phi_a, _amp_a = common.get_Rangle(data["annual"])
    phi_e, amp_e = common.get_Rangle(data["enso"])
    phase_origin = data["phase_origin"]

    single = args.amp is not None or args.ratio is not None
    if single:
        if args.amp is None or args.ratio is None:
            raise ValueError("single run requires both --amp and --ratio")
        estimate_one(direct, common, args, phi_a, phi_e, amp_e, phase_origin, args.amp, args.ratio, args.order)
        return 0

    if args.output is not None:
        raise ValueError("--output is only for a single run; use --output-dir for batch runs.")

    outputs = []
    orders = [args.order] if args.order is not None else list(args.orders)
    total = len(args.amps) * len(args.ratios) * len(orders)
    count = 0
    for amp in args.amps:
        for ratio in args.ratios:
            for order in orders:
                count += 1
                label = "auto" if order is None else str(order)
                print(f"[{count}/{total}] amp={amp:g}, ratio={ratio[0]}:{ratio[1]}, M={label}")
                outputs.append(estimate_one(direct, common, args, phi_a, phi_e, amp_e, phase_origin, amp, ratio, order))
    print(f"batch done: {len(outputs)} files saved in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
