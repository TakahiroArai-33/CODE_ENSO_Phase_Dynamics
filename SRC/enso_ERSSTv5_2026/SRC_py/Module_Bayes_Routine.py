#!/usr/bin/env python3
# coding: utf-8

"""
Minimal routine for estimating Q(phi1, phi2) from phase time series.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Module_Bayes_Estimation import BayesEstimate, estimate_components


TWOPI = 2.0 * np.pi


@dataclass
class EstimationData:
    """
    Data class for storing estimation data.
    """
    phi1_ast: np.ndarray
    phi2_ast: np.ndarray
    dot_phi1: np.ndarray
    dot_phi2: np.ndarray
    keep: np.ndarray


@dataclass
class RoutineResult:
    """
    estimate_q の返り値をまとめるデータクラス。

    data には推定に実際に使われた時刻データを保存し、
    phi1, phi2 にはそれぞれの位相方程式の推定結果を保存する。

    EstimationDataクラス: make_estimation_data で作った、filtering 後の推定用データ。
    BayesEstimateクラス: estimate_components で作った、dot_phi1, dot_phi2 の位相方程式の推定結果。
    """
    data: EstimationData
    phi1: BayesEstimate
    phi2: BayesEstimate

    @property
    def noise_cov(self):
        """
        dot_chi = [dot_phi1, dot_phi2] の同時推定で得たノイズ共分散行列 E を返す。
        phi1, phi2 の BayesEstimate は同じ E を保持している。
        """
        return self.phi1.noise_cov

    def noise_cholesky(self):
        """
        ノイズ共分散行列 E の Cholesky 分解を返す。
        """
        return np.linalg.cholesky(self.noise_cov)

    def predict(self, phi1, phi2, return_std: bool = False):
        """
        任意の (phi1, phi2) に対して、
        推定された Q_1(phi1, phi2), Q_2(phi1, phi2)を評価する。

        実際の評価処理は Module_Bayes_Estimation.BayesEstimate.predict に任せる。
        ここでは phi1 成分と phi2 成分の結果をまとめて返すだけ。
        """
        return (
            self.phi1.predict(phi1, phi2, return_std=return_std),
            self.phi2.predict(phi1, phi2, return_std=return_std),
        )


def make_estimation_data(
    phi1,
    phi2,
    dt: float,
    amp=None,
    amp_rank: float = 0.0,
) -> EstimationData:
    """
    phi1, phi2: 1D array-like, shape (N,)
        Phase time series (in radians).
    dt: float
        Time step of the time series.
    つぎの引数はオプションで、データ選択(filtering)のために使用されます。
    amp: 1D array-like, shape (N,), optional
        Amplitude time series (for data selection).
    amp_rank: float, optional
        Percentile rank (0-100) for amplitude-based data selection. Only points with amplitude
    """
    phi1 = np.asarray(phi1, dtype=float).reshape(-1)
    phi2 = np.asarray(phi2, dtype=float).reshape(-1)
    if phi1.size != phi2.size:
        raise ValueError("phi1 and phi2 must have the same length.")
    if phi1.size < 2:
        raise ValueError("phi1 and phi2 must contain at least two points.")

    #* 中点
    phi1_ast = np.mod(0.5 * (phi1[:-1] + phi1[1:]), TWOPI)
    phi2_ast = np.mod(0.5 * (phi2[:-1] + phi2[1:]), TWOPI)
    #* 微分値
    dot_phi1 = np.diff(phi1) / float(dt)
    dot_phi2 = np.diff(phi2) / float(dt)

    keep = np.ones(phi1_ast.size, dtype=bool)
    if amp is not None and amp_rank > 0:
        amp = np.asarray(amp, dtype=float).reshape(-1)
        if amp.size != phi1.size:
            raise ValueError("amp must have the same length as phi1 and phi2.")
        #* 中点の振幅 (近似値)
        amp_ast = 0.5 * (amp[:-1] + amp[1:])
        #* (filtering) 振幅時系列に対して、amp_rankパーセンタイル以上の点を選択
        keep = amp_ast >= np.percentile(amp_ast, float(amp_rank))

    return EstimationData(
        phi1_ast=phi1_ast[keep],
        phi2_ast=phi2_ast[keep],
        dot_phi1=dot_phi1[keep],
        dot_phi2=dot_phi2[keep],
        keep=keep,
    )


def estimate_q(
    phi1,
    phi2,
    dt: float,
    amp=None,
    amp_rank: float = 0.0,
    M1: int = 10,
    M2: int = 10,
    lambda_0: float = 1e-12,
    num_iterations: int = 100,
    rtol: float = 1e-8,
) -> RoutineResult:
    """
    2 本の位相時系列から ``Q(phi1, phi2)`` の 2 成分を推定する。

    この関数の入力 ``phi1`` と ``phi2`` は時刻ごとの位相である。内部では
    隣り合う時刻の中点位相を説明変数にし、差分から求めた位相速度を
    目的変数にする。つまり、
    ``dot_phi1 = Q1(phi1_mid, phi2_mid)`` と
    ``dot_phi2 = Q2(phi1_mid, phi2_mid)`` をデータから推定している。
    """
    
    #* データの準備
    data = make_estimation_data(
        phi1=phi1,
        phi2=phi2,
        dt=dt,
        amp=amp,
        amp_rank=amp_rank,
    )
    #* dot{phi_1} = Q_1(phi_1, phi_2) + noise
    #* dot{phi_2} = Q_2(phi_1, phi_2) + noise
    result_phi1, result_phi2 = estimate_components(
        phi1_ast=data.phi1_ast,
        phi2_ast=data.phi2_ast,
        dot_phi1=data.dot_phi1,
        dot_phi2=data.dot_phi2,
        M1=M1,
        M2=M2,
        lambda_0=lambda_0,
        dt=dt,
        num_iterations=num_iterations,
        rtol=rtol,
    )
    return RoutineResult(data=data, phi1=result_phi1, phi2=result_phi2)
