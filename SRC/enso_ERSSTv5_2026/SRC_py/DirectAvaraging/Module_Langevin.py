#!/usr/bin/env python3
# coding: utf-8

"""Minimal Langevin simulation utilities."""

from __future__ import annotations

import numpy as np


def simulate_langevin(
    func_a,
    func_e,
    t_span,
    initial_state,
    noise_cholesky=None,
    rng=None,
):
    """
    Euler-Maruyamaで位相方程式をシミュレーションする。

    Returns
    -------
    tuple[ndarray, ndarray]
        phi_a, phi_e.
    """
    t_span = np.asarray(t_span, dtype=float).reshape(-1)
    if t_span.size < 1:
        raise ValueError("t_span must contain at least one time.")

    phi = np.empty((t_span.size, 2), dtype=float)
    phi[0] = np.asarray(initial_state, dtype=float).reshape(2)

    H = np.zeros((2, 2), dtype=float) if noise_cholesky is None else np.asarray(noise_cholesky, dtype=float)
    if rng is None:
        rng = np.random.default_rng()

    for i in range(t_span.size - 1):
        dt = t_span[i + 1] - t_span[i]
        y = phi[i]
        drift = np.array([func_a(y[0], y[1]), func_e(y[0], y[1])], dtype=float).reshape(2)
        noise = H @ (np.sqrt(dt) * rng.standard_normal(2))
        phi[i + 1] = y + drift * dt + noise

    return phi[:, 0], phi[:, 1]


def simulate_langevin_averaged(
    func_gamma_a,
    func_gamma_e,
    t_span,
    initial_state,
    ratio,
    noise_cholesky=None,
    rng=None,
):
    """
    平均化した位相方程式をEuler-Maruyamaでシミュレーションする。

    psi_a = phi_a - ratio * phi_e,
    psi_e = ratio * phi_e - phi_a として、
        dot{phi_a} = Gamma_a(psi_a) + noise_a
        dot{phi_e} = Gamma_e(psi_e) + noise_e
    を計算する。
    """
    t_span = np.asarray(t_span, dtype=float).reshape(-1)
    if t_span.size < 1:
        raise ValueError("t_span must contain at least one time.")

    phi = np.empty((t_span.size, 2), dtype=float)
    phi[0] = np.asarray(initial_state, dtype=float).reshape(2)

    H = np.zeros((2, 2), dtype=float) if noise_cholesky is None else np.asarray(noise_cholesky, dtype=float)
    if rng is None:
        rng = np.random.default_rng()

    dt = np.diff(t_span)
    noise = (np.sqrt(dt)[:, None] * rng.standard_normal((dt.size, 2))) @ H.T

    for i, dt_i in enumerate(dt):
        phi_a, phi_e = phi[i]
        psi_a = phi_a - ratio * phi_e
        psi_e = ratio * phi_e - phi_a
        drift_a = np.asarray(func_gamma_a(psi_a), dtype=float).item()
        drift_e = np.asarray(func_gamma_e(psi_e), dtype=float).item()
        phi[i + 1, 0] = phi_a + drift_a * dt_i + noise[i, 0]
        phi[i + 1, 1] = phi_e + drift_e * dt_i + noise[i, 1]

    return phi[:, 0], phi[:, 1]
