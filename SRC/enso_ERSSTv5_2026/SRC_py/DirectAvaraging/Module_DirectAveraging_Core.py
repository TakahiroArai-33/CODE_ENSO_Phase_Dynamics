#!/usr/bin/env python3
# coding: utf-8

"""Direct averaging estimation core for the corrected ENSO phase convention.

This module intentionally contains no phase-data loader and no batch CLI.
Phase preparation is owned by ``../common_fixed.py`` and orchestration is owned by
``../Estimate_DirectAveraging_fixed.py``.  The functions here are the estimation
and result-saving pieces copied from the original DirectAvaraging workflow.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
TWOPI = 2.0 * math.pi
DT_YEAR = 1.0 / 12.0
DEFAULT_M_MIN = 0
DEFAULT_M_MAX = 10

np = None
data_arrangement = None
routine_minibatch_update = None
simulate_langevin_averaged = None


def load_runtime_dependencies() -> None:
    """Import heavy runtime dependencies after the caller has set up ``sys.path``."""
    global np, data_arrangement, routine_minibatch_update, simulate_langevin_averaged

    import numpy as _np

    sys.path.insert(0, str(SCRIPT_DIR))
    from Module_Langevin import simulate_langevin_averaged as _simulate_langevin_averaged
    from Module_routine import data_arrangement as _data_arrangement
    from Module_routine import routine_minibatch_update as _routine_minibatch_update

    np = _np
    data_arrangement = _data_arrangement
    routine_minibatch_update = _routine_minibatch_update
    simulate_langevin_averaged = _simulate_langevin_averaged


def number_label(value: float) -> str:
    """Format a numeric option for stable output-file names."""
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}".replace(".", "p")


def order_label(order: int | None) -> str:
    """Return the output-file label for fixed-order or automatic model selection."""
    return "Mauto" if order is None else f"M{order}"


def order_range(order: int | None) -> tuple[int, int]:
    """Convert an optional fixed order into the estimator's ``M_min, M_max`` range."""
    if order is None:
        return DEFAULT_M_MIN, DEFAULT_M_MAX
    return order, order


@dataclass
class EstimationData:
    """Data used by the DirectAveraging estimator after amplitude filtering."""

    phi1_list: list[np.ndarray]
    phi2_list: list[np.ndarray]
    index_list: list[np.ndarray]
    keep: np.ndarray

    @property
    def keep_pair_indices(self) -> np.ndarray:
        """Return indices of adjacent phase increments kept for estimation."""
        return np.flatnonzero(self.keep)


def make_estimation_data(
    phi1,
    phi2,
    amp=None,
    amp_rank: float = 0.0,
) -> EstimationData:
    """
    Prepare DirectAveraging data with the same amplitude-filter style as Q.

    The keep mask is defined on adjacent phase increments. DirectAveraging still
    needs continuous phase segments, so consecutive kept increments are converted
    back to point-index segments before estimation.
    """
    phi1 = np.asarray(phi1, dtype=float).reshape(-1)
    phi2 = np.asarray(phi2, dtype=float).reshape(-1)
    if phi1.size != phi2.size:
        raise ValueError("phi1 and phi2 must have the same length.")
    if phi1.size < 2:
        raise ValueError("phi1 and phi2 must contain at least two points.")

    keep = np.ones(phi1.size - 1, dtype=bool)
    if amp is not None and amp_rank > 0:
        amp = np.asarray(amp, dtype=float).reshape(-1)
        if amp.size != phi1.size:
            raise ValueError("amp must have the same length as phi1 and phi2.")
        amp_ast = 0.5 * (amp[:-1] + amp[1:])
        keep = amp_ast >= np.percentile(amp_ast, float(amp_rank))

    pair_indices = np.flatnonzero(keep)
    if pair_indices.size == 0:
        raise ValueError("No adjacent phase increment remains after amplitude filtering.")

    split_points = np.where(np.diff(pair_indices) != 1)[0] + 1
    phi1_list: list[np.ndarray] = []
    phi2_list: list[np.ndarray] = []
    index_list: list[np.ndarray] = []
    for segment_pairs in np.split(pair_indices, split_points):
        segment_indices = np.arange(segment_pairs[0], segment_pairs[-1] + 2)
        phi1_list.append(phi1[segment_indices])
        phi2_list.append(phi2[segment_indices])
        index_list.append(segment_indices)

    return EstimationData(
        phi1_list=phi1_list,
        phi2_list=phi2_list,
        index_list=index_list,
        keep=keep,
    )


def ratio_to_phase_coeffs(ratio: tuple[int, int]) -> tuple[float, float]:
    """Convert an integer synchronization ratio into phase-argument coefficients."""
    m, n = ratio
    return 1.0, n / m


def evidence_rows(evidence_part: dict[str, object], key: str) -> np.ndarray:
    """Convert evidence dictionaries to a numeric array for saving in NPZ files."""
    rows = evidence_part[key]
    return np.array([
        [float(row["M"]), float(row["lambda_0"]), float(row["log10_evidence"])]
        for row in rows
    ])


def simulate_langevin_direct(
    result: dict[str, object],
    phase_coeffs: tuple[float, float],
    initial_phase: np.ndarray,
    step_count: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """Simulate the averaged DirectAveraging phase equation from a fitted result."""
    x_grid = np.asarray(result["x"], dtype=float)
    y1 = np.asarray(result["y1"], dtype=float)
    y2 = np.asarray(result["y2"], dtype=float)
    p1, p2 = phase_coeffs
    step_count = max(2, int(step_count))
    t_langevin = DT_YEAR * np.arange(step_count)
    noise_cholesky = np.diag([float(result["D1"]), float(result["D2"])])
    rng = np.random.default_rng(seed)

    def func_gamma_a(psi):
        """Evaluate the fitted annual coupling function on wrapped phase difference."""
        return np.interp(np.mod(psi, TWOPI), x_grid, y1)

    def func_gamma_e(psi):
        """Evaluate the fitted ENSO coupling function on wrapped phase difference."""
        return np.interp(np.mod(psi, TWOPI), x_grid, y2)

    if not np.isclose(p1, 1.0):
        raise ValueError("simulate_langevin_averaged assumes phase_coeffs[0] == 1.0.")

    phi_a_langevin, phi_e_langevin = simulate_langevin_averaged(
        func_gamma_a,
        func_gamma_e,
        t_langevin,
        initial_phase,
        ratio=p2,
        noise_cholesky=noise_cholesky,
        rng=rng,
    )

    return {
        "t_langevin": t_langevin,
        "phi_a_langevin": phi_a_langevin,
        "phi_e_langevin": phi_e_langevin,
        "delta_phi_langevin": np.mod(
            p1 * phi_a_langevin - p2 * phi_e_langevin,
            TWOPI,
        ),
    }


def estimate_direct(
    args: argparse.Namespace,
    phase_origin: dict[str, object],
    data: EstimationData,
    initial_phase: np.ndarray,
    data_length: int,
) -> tuple[dict[str, object], dict[str, object], EstimationData]:
    """Fit DirectAveraging models from prepared estimation data."""
    p1, p2 = ratio_to_phase_coeffs(args.ratio)
    M_min, M_max = order_range(args.order)
    result_model_1, result_model_2, model_evidence = routine_minibatch_update(
        data.phi1_list,
        data.phi2_list,
        dt=DT_YEAR,
        p1=p1,
        p2=p2,
        M_min=M_min,
        M_max=M_max,
        log10_lambda_0_min=args.log10_lambda_min,
        log10_lambda_0_max=args.log10_lambda_max,
        alpha_0=args.alpha_0,
        beta_0=args.beta_0,
        return_evidence=True,
    )
    result = data_arrangement(result_model_1, result_model_2, DT_YEAR)
    result["lambda_1"] = float(result_model_1.lambda_0)
    result["lambda_2"] = float(result_model_2.lambda_0)
    result["model_evidence"] = model_evidence
    result["phase_origin"] = phase_origin
    result["keep_pair_indices"] = data.keep_pair_indices
    result.update(
        simulate_langevin_direct(
            result,
            phase_coeffs=(p1, p2),
            initial_phase=initial_phase,
            step_count=args.langevin_steps if args.langevin_step_factor is None else args.langevin_step_factor * data_length,
            seed=args.langevin_seed,
        )
    )
    return result, phase_origin, data


def save_result(
    args: argparse.Namespace,
    result: dict[str, object],
    phase_origin: dict[str, object],
    data: EstimationData,
) -> Path:
    """Save a fitted DirectAveraging result and diagnostics to a compressed NPZ."""
    ratio_m, ratio_n = args.ratio
    ratio_key = f"{ratio_m}:{ratio_n}"
    amp_label = number_label(args.amp)
    M_min, M_max = order_range(args.order)
    order_name = order_label(args.order)
    output_path = args.output
    if output_path is None:
        output_path = args.output_dir / f"DirectAvaraging_result_{order_name}_amp{amp_label}_ratio{ratio_key}.npz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    p1, p2 = ratio_to_phase_coeffs(args.ratio)
    data_phi_a = np.concatenate(data.phi1_list)
    data_phi_e = np.concatenate(data.phi2_list)
    valid_indices = np.concatenate(data.index_list)
    keep_pair_indices = data.keep_pair_indices
    model_evidence = result["model_evidence"]
    mode_indices = phase_origin.get("mode_indices", {})
    y1 = np.asarray(result["y1"], dtype=float) / TWOPI
    y2 = np.asarray(result["y2"], dtype=float) / TWOPI
    y3 = np.asarray(result["y3"], dtype=float) / TWOPI
    y1_err = np.asarray(result["y1_err"], dtype=float) / TWOPI
    y2_err = np.asarray(result["y2_err"], dtype=float) / TWOPI
    y3_err = np.asarray(result["y3_err"], dtype=float) / TWOPI

    np.savez_compressed(
        output_path,
        x=np.asarray(result["x"]),
        y1=y1,
        y2=y2,
        y3=y3,
        y1_err=y1_err,
        y2_err=y2_err,
        y3_err=y3_err,
        y_unit=np.array("freq/year"),
        M=np.array(-1 if args.order is None else args.order, dtype=int),
        M_mode=np.array("auto" if args.order is None else "fixed"),
        M_min=np.array(M_min, dtype=int),
        M_max=np.array(M_max, dtype=int),
        M1=np.array(result["M1"], dtype=int),
        M2=np.array(result["M2"], dtype=int),
        amp_rank=np.array(args.amp, dtype=float),
        ratio=np.array([ratio_m, ratio_n], dtype=int),
        ratio_key=np.array(ratio_key),
        phase_argument_coeffs=np.array([p1, p2], dtype=float),
        D1=np.array(result["D1"], dtype=float),
        D2=np.array(result["D2"], dtype=float),
        D3=np.array(result["D3"], dtype=float),
        omega1=np.array(result["omega1"], dtype=float),
        omega2=np.array(result["omega2"], dtype=float),
        lambda_1=np.array(result["lambda_1"], dtype=float),
        lambda_2=np.array(result["lambda_2"], dtype=float),
        log10_lambda_range=np.array([args.log10_lambda_min, args.log10_lambda_max], dtype=int),
        dt=np.array(DT_YEAR, dtype=float),
        data_phi_a=np.mod(data_phi_a, TWOPI),
        data_phi_e=np.mod(data_phi_e, TWOPI),
        data_delta_phi=np.mod(p1 * data_phi_a - p2 * data_phi_e, TWOPI),
        delta_phi_data=np.mod(p1 * data_phi_a - p2 * data_phi_e, TWOPI),
        t_langevin=np.asarray(result["t_langevin"], dtype=float),
        phi_a_langevin=np.mod(np.asarray(result["phi_a_langevin"], dtype=float), TWOPI),
        phi_e_langevin=np.mod(np.asarray(result["phi_e_langevin"], dtype=float), TWOPI),
        delta_phi_langevin=np.asarray(result["delta_phi_langevin"], dtype=float),
        delta_phi_averaged_langevin=np.asarray(result["delta_phi_langevin"], dtype=float),
        distribution_bins=np.linspace(0.0, TWOPI, num=13),
        langevin_steps=np.array(args.langevin_steps, dtype=int),
        langevin_step_factor=np.array(-1 if args.langevin_step_factor is None else args.langevin_step_factor, dtype=int),
        langevin_seed=np.array(args.langevin_seed, dtype=int),
        keep=data.keep,
        valid_indices=valid_indices,
        keep_pair_indices=keep_pair_indices,
        segment_lengths=np.array([phi.size for phi in data.phi1_list], dtype=int),
        phi_a_shift_rad=np.array(phase_origin["phi_a_shift_rad"], dtype=float),
        phi_e_shift_rad=np.array(phase_origin["phi_e_shift_rad"], dtype=float),
        annual_mode_index=np.array(mode_indices.get("annual", -1), dtype=int),
        enso_mode_index=np.array(mode_indices.get("enso", -1), dtype=int),
        enso_peak_percentile=np.array(phase_origin["enso_peak_percentile"], dtype=float),
        enso_peak_mask=np.asarray(phase_origin["enso_peak_mask"], dtype=bool),
        phase_convention=np.array(phase_origin["phase_convention"]),
        evidence_model1_all=evidence_rows(model_evidence["model1"], "all"),
        evidence_model1_best_by_M=evidence_rows(model_evidence["model1"], "best_by_M"),
        evidence_model2_all=evidence_rows(model_evidence["model2"], "all"),
        evidence_model2_best_by_M=evidence_rows(model_evidence["model2"], "best_by_M"),
    )
    return output_path


if __name__ == "__main__":
    raise SystemExit("Use ../Estimate_DirectAveraging_fixed.py. This module only contains estimation-core functions.")
