#!/usr/bin/env python3
# coding: utf-8

"""
Minimal CPU estimator for Q(phi1, phi2).

Model:
    dot_phi = Q(phi1, phi2)

Q is represented by a two-dimensional complex Fourier basis:
    exp(i (m phi1 + n phi2))
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy import linalg as LA


TWOPI = 2.0 * np.pi


@dataclass
class BayesEstimate:
    """
    1つの位相方程式成分 Q_j(phi1, phi2) のベイズ推定結果。

    target:
        推定対象の成分名。"phi1" または "phi2"。
    modes:
        複素フーリエ基底 exp(i(m1 phi1 + m2 phi2)) の (m1, m2) 配列。
        modes[0] = (0, 0) なので coef[0] が定数項に対応する。
    coef:
        Q_j の複素フーリエ係数。Q_j の評価値は実部を物理量として使う。
    cov:
        coef の事後共分散行列。
    noise_cov:
        同時推定で得たノイズ共分散行列 E。
    lambda_0:
        事前精度行列 Xi_prior の L2 正則化係数。
    """
    target: str
    modes: np.ndarray
    coef: np.ndarray
    cov: np.ndarray
    noise_cov: np.ndarray
    lambda_0: float

    def predict(self, phi1, phi2, return_std: bool = False):
        """
        推定済みのパラメータから Q_j(phi1, phi2) を計算する。
        複素基底を使うため、内部では X @ coef は複素数になり得る。
        位相方程式の値としては、その実部 mean = Re[X @ coef] を返す。
        return_std=True の場合:
            事後共分散 cov から sqrt(Re[X^H cov X]) を標準偏差として返す。
            これは係数推定の不確かさであり、ノイズ共分散 E そのものではない。
            E は self.noise_cov に保存されている。
        """
        X = design_matrix(phi1, phi2, self.modes)
        mean = np.real(X @ self.coef)
        if not return_std:
            return mean
        #* var は Q_j(phi1, phi2) の係数事後分散。ノイズ共分散 E は self.noise_cov。
        var = np.einsum("ij,jk,ik->i", np.conj(X), self.cov, X)
        return mean, np.sqrt(np.maximum(np.real(var), 0.0))


def make_modes(M1: int = 10, M2: int = 10) -> np.ndarray:
    """ 
    (m1, m2) の組み合わせを作成する。
    return shape: (num_modes, 2)
        mode[:, 0] = m1のarray
        mode[:, 1] = m2のarray

    注意:
        modes[0] = (0, 0) として定数項を先頭に置く。
        これにより、推定後の coef[0] が exp(i*0)=1 に掛かる定数項になる。
    """
    modes = [(0, 0)]
    for m in range(-M1, M1 + 1):
        for n in range(-M2, M2 + 1):
            if m == 0 and n == 0:
                continue
            modes.append((m, n))
    return np.asarray(modes, dtype=int)


def design_matrix(phi1, phi2, modes: np.ndarray) -> np.ndarray:
    """
    位相方程式を次の複素2次元フーリエ級数展開で表現するための計画行列を作成する。
    sum_{m1, m2} a_{m1,m2} exp(i (m1 phi1 + m2 phi2))
    """
    phi1 = np.mod(np.asarray(phi1, dtype=float).reshape(-1), TWOPI)
    phi2 = np.mod(np.asarray(phi2, dtype=float).reshape(-1), TWOPI)
    if phi1.size != phi2.size:
        raise ValueError("phi1 and phi2 must have the same length.")

    #* modes[0] = (0, 0) なので、返り値の第0列は exp(i*0)=1 の定数基底。
    #* exp(i (m1 phi1 + m2 phi2)) = exp(i m1 phi1) * exp(i m2 phi2) なので expの中身phaseを外積で計算。
    phase = np.outer(phi1, modes[:, 0]) + np.outer(phi2, modes[:, 1])
    return np.exp(1j * phase)


def estimate_components(
    phi1_ast,
    phi2_ast,
    dot_phi1,
    dot_phi2,
    M1: int = 10,
    M2: int = 10,
    lambda_0: float = 1e-6,
    dt: float = 1.0,
    num_iterations: int = 100,
    rtol: float = 1e-7,
) -> tuple[BayesEstimate, BayesEstimate]:
    """
    位相方程式成分 Q_j(phi1, phi2) のベイズ推定を行う。
     (フィルタリングされたデータが使用可能)
    phi1_ast, phi2_ast:
        観測された位相データの配列。
    dot_phi1, dot_phi2:
        観測された位相変化率の配列。
    M1, M2:
        フーリエ基底の最大次数。
    lambda_0:
        事前精度行列 Xi_prior の L2 正則化係数。
        (事前分布の分散の設定。平均値は0としている。)
    dt:
        時系列データをサンプルした時間刻み。
    num_iterations:
        推定の反復回数。
    rtol:
        収束判定の相対許容誤差。
    return:
        BayesEstimate オブジェクトのタプル (phi1 成分, phi2 成分)。
    """
    phi1_ast = np.asarray(phi1_ast, dtype=float).reshape(-1)
    phi2_ast = np.asarray(phi2_ast, dtype=float).reshape(-1)
    dot_phi1 = np.asarray(dot_phi1, dtype=float).reshape(-1)
    dot_phi2 = np.asarray(dot_phi2, dtype=float).reshape(-1)
    if not (phi1_ast.size == phi2_ast.size == dot_phi1.size == dot_phi2.size):
        raise ValueError("phi1_ast, phi2_ast, dot_phi1, and dot_phi2 must have the same length.")
    Delta_t = float(dt)
    if Delta_t <= 0.0:
        raise ValueError("dt must be positive.")

    modes = make_modes(M1=M1, M2=M2)
    g_array = design_matrix(phi1_ast, phi2_ast, modes)  # size(N, R)

    #* 今回は chi = (phi1, phi2) として、dot_chi = [dot_phi1, dot_phi2] を同時推定する。
    dot_chi = np.array([dot_phi1, dot_phi2], dtype=float)

    #* R: 1成分あたりの基底関数の総数。a = [a_phi1, a_phi2] のサイズは 2R。
    R = g_array.shape[1]
    N = dot_phi1.size

    #* G_normal: size(2, 2R, N)
    #* G_normal[:, :, n] = [[g_n^T, 0], [0, g_n^T]]
    G_normal = np.zeros((2, 2 * R, N), dtype=np.complex128)
    G_normal[0, :R, :] = g_array.T
    G_normal[1, R:, :] = g_array.T

    #* G_partial: size(2, 2R, N)
    #* phi1成分は partial_phi1 G、phi2成分は partial_phi2 G。
    #* d exp(i(m1 phi1 + m2 phi2))/d phi_j = i m_j exp(i(m1 phi1 + m2 phi2))
    G_partial = np.zeros_like(G_normal)
    G_partial[0, :R, :] = (1j * modes[:, 0])[:, None] * g_array.T
    G_partial[1, R:, :] = (1j * modes[:, 1])[:, None] * g_array.T

    datas = {
        "G_normal": G_normal,
        "G_partial": G_partial,
        "dot_chi": dot_chi,
    }

    #* 属性として保持しているテンソルを取り出す
    G_normal = datas["G_normal"]  # size(2, 2R, N)
    G_partial = datas["G_partial"]  # size(2, 2R, N)
    dot_chi = datas["dot_chi"]  # size(2, N)

    #* 事前精度行列 Xi_prior。
    Xi_prior = (float(lambda_0)) * np.eye(2 * R, dtype=np.complex128)

    #* --- Tool: Hermitian (conjugate) transpose ---
    G_normal_H = np.conj(np.transpose(G_normal, (1, 0, 2)))  # size(2R, 2, N)
    G_partial_H = np.conj(np.transpose(G_partial, (1, 0, 2)))  # size(2R, 2, N)

    #* 事前平均は反復中も0に固定し、反復値 a_old と区別する。
    a_prior = np.zeros(2 * R, dtype=np.complex128)
    a_old = a_prior.copy()

    for _ in range(int(num_iterations)):
        #* --- Residuals for E update ---
        #* G_normal: size(2, 2R, N), a_old: size(2R,)
        Ga = np.tensordot(G_normal, a_old, axes=([1], [0]))  # size(2, N)
        residual = (dot_chi - Ga).real

        #* --- E update ---
        #* residual は実数としたので、式におけるdaggerは転置で良い。
        E_new_complex = (Delta_t / N) * (residual @ residual.T)
        E_new = ((E_new_complex + E_new_complex.T) / 2).real.astype(np.float64)
        Einv = LA.inv(E_new)

        #* --- Xi update ---
        #* tensordot で E_inv を各 G_n に同時適用する。
        Einv_G = np.tensordot(Einv, G_normal, axes=([1], [0]))  # size(2, 2R, N)
        Xi_new = Xi_prior + Delta_t * np.einsum("kin,iqn->kq", G_normal_H, Einv_G, optimize=True)

        #* 数値誤差で非エルミート化を抑える
        Xi_new = (Xi_new + Xi_new.conj().T) / 2

        #* --- r calculate ---
        e_vec = np.ones(2, dtype=float)
        E_inv_dotchi = Einv @ dot_chi  # size(2, N)
        sum_GH_Einv_dotchi = np.einsum("kin,in->k", G_normal_H, E_inv_dotchi, optimize=True)
        sum_partial = np.einsum("kin,i->k", G_partial_H, e_vec, optimize=True)
        r = Xi_prior @ a_prior + Delta_t * sum_GH_Einv_dotchi - (Delta_t / 2.0) * sum_partial

        #* --- a update (Xi @ a = r ↔ a = Xi^{-1} @ r)---
        a_new = LA.solve(Xi_new, r)

        #* 係数の相対変化が十分小さければ収束とみなす。
        denom = max(LA.norm(a_old), 1e-8)
        if LA.norm(a_new - a_old) / denom < float(rtol):
            a_old = a_new
            break
        a_old = a_new

    #* 最終係数で E を評価する。
    a = a_old
    Ga = np.tensordot(G_normal, a, axes=([1], [0]))  # size(2, N)
    residual = (dot_chi - Ga).real
    E_new_complex = (Delta_t / N) * (residual @ residual.T)
    E_new = ((E_new_complex + E_new_complex.T) / 2).real.astype(np.float64)

    #* 係数共分散は、最終Eで評価した精度行列の逆行列。
    Einv = LA.inv(E_new)
    Einv_G = np.tensordot(Einv, G_normal, axes=([1], [0]))  # size(2, 2R, N)
    Xi_new = Xi_prior + Delta_t * np.einsum("kin,iqn->kq", G_normal_H, Einv_G, optimize=True)
    Xi_new = (Xi_new + Xi_new.conj().T) / 2
    Sigma = LA.inv(Xi_new)

    #* --- 推定後の係数の取り出し ---
    #* a は [Q1の係数R個, Q2の係数R個] の順に並んでいる。
    #* それぞれの先頭 coef[0] が modes[0]=(0,0) に対応する定数項。
    result_phi1 = BayesEstimate(
        target="phi1",
        modes=modes,
        #* a = [a_phi1, a_phi2] なので、前半 R 個が dot_phi1 = Q1(phi1, phi2) の係数。
        #* modes[0] = (0, 0) のため、coef[0] が Q1 の定数項。
        coef=a[:R],
        cov=Sigma[:R, :R],
        noise_cov=E_new,
        lambda_0=float(lambda_0),
    )
    result_phi2 = BayesEstimate(
        target="phi2",
        modes=modes,
        #* a = [a_phi1, a_phi2] なので、後半 R 個が dot_phi2 = Q2(phi1, phi2) の係数。
        #* modes[0] = (0, 0) のため、coef[0] が Q2 の定数項。
        coef=a[R:],
        cov=Sigma[R:, R:],
        noise_cov=E_new,
        lambda_0=float(lambda_0),
    )
    return result_phi1, result_phi2
