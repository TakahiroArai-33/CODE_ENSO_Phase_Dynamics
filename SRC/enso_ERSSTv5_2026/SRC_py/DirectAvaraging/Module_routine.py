import numpy as np
import numpy.linalg as LA
import math
from copy import deepcopy, copy
from Module_Estimator import Estimator


log10_lambda_0_min = 1
log10_lambda_0_max = 8
M_min = int(0)
M_max = int(10)
alpha_0 = 2.0
beta_0 = 2.0
lambda_0 = 2.0



def routine(phi1, phi2, dt=1, 
            p1=int(1), p2=int(1),
            M_min=M_min, M_max=M_max,
            log10_lambda_0_min=log10_lambda_0_min, 
            log10_lambda_0_max=log10_lambda_0_max, 
            alpha_0=alpha_0, beta_0=beta_0,
            return_evidence=False
            ):
    """Fit two coupled phase models from one continuous pair of phase series."""

    print("initial parameter setting")
    print("M_min={:d}, M_max={:d}".format(M_min, M_max))
    print("log10_lambda_0_min={:d}, log10_lambda_0_max={:d}".\
          format(int(log10_lambda_0_min), int(log10_lambda_0_max)))
    print("alpha_0={:0.1e}, beta_0={:0.1e}".format(alpha_0, beta_0))
    print("!! setting is for the investigation for p1:p2 sync.:  {:g}*phi1 and {:g}*phi2".format(p1, p2))
    print("")


    # 推定に用いるデータを作成
    phase_series = np.array([phi1, 
                            phi2]).T

    L = phi1.size # データサイズ
    N = 2 # 振動子数


    # 平均速度(初期条件として与える) [rad/day]
    v_1_mean = (phi1[-1]-phi1[0]) / (dt*(L-1))
    v_2_mean = (phi2[-1]-phi2[0]) / (dt*(L-1))


    ### 推定: receiver=1
    Estimator_1 = Estimator(receiver=0, N=2, p=p1, q=p2,
                                            M_min=M_min, M_max=M_max,
                                            log10_lambda_0_min=log10_lambda_0_min, 
                                            log10_lambda_0_max=log10_lambda_0_max, 
                                            alpha_0=alpha_0, beta_0=beta_0)
    # set intial condition (initial condition of omega)
    Estimator_1.set_init_all(v_1_mean)
    # update
    Estimator_1.update_all(phase_series, delta_t=dt)
    # model section
    result_model_1, _ = Estimator_1.model_selection()


    ### 推定: receiver=2
    Estimator_2 = Estimator(receiver=1, N=2, p=p2, q=p1,
                                            M_min=M_min, M_max=M_max,
                                            log10_lambda_0_min=log10_lambda_0_min, 
                                            log10_lambda_0_max=log10_lambda_0_max,  
                                            alpha_0=alpha_0, beta_0=beta_0)
    # set intial condition (initial condition of omega)
    Estimator_2.set_init_all(v_2_mean)
    # update
    Estimator_2.update_all(phase_series, delta_t=dt)
    # model section
    result_model_2, _ = Estimator_2.model_selection()

    if return_evidence:
        evidence = evidence_arrangement(Estimator_1, Estimator_2)
        return result_model_1, result_model_2, evidence
    return result_model_1, result_model_2


def _evidence_by_M(estimator, sender):
    """
    Return model evidence candidates and the best lambda_0 for each M.
    """
    model_list, log10_evidence_list = estimator.all_evidence()

    rows = []
    for model, log10_evidence in zip(model_list, log10_evidence_list):
        rows.append({
            "M": int(model.return_Gamma_order(sender=sender)),
            "lambda_0": float(model.lambda_0),
            "log10_evidence": float(log10_evidence),
        })

    best_by_M = []
    for M in sorted({row["M"] for row in rows}):
        candidates = [row for row in rows if row["M"] == M]
        best_by_M.append(max(candidates, key=lambda row: row["log10_evidence"]))

    return {"all": rows, "best_by_M": best_by_M}


def evidence_arrangement(estimator_1, estimator_2):
    """
    Arrange all candidate model evidence for the two fitted estimators.
    """
    return {
        "model1": _evidence_by_M(estimator_1, sender=1),
        "model2": _evidence_by_M(estimator_2, sender=0),
    }



def routine_minibatch_update(phi1_list:list, phi2_list:list, dt=1, 
            p1=int(1), p2=int(1),
            M_min=M_min, M_max=M_max,
            log10_lambda_0_min=log10_lambda_0_min, 
            log10_lambda_0_max=log10_lambda_0_max, 
            alpha_0=alpha_0, beta_0=beta_0,
            return_evidence=False
            ):
    """Fit coupled phase models by accumulating several filtered data segments."""

    # print("initial parameter setting")
    # print("M_min={:d}, M_max={:d}".format(M_min, M_max))
    # print("log10_lambda_0_min={:d}, log10_lambda_0_max={:d}".\
    #       format(int(log10_lambda_0_min), int(log10_lambda_0_max)))
    # print("alpha_0={:0.1e}, beta_0={:0.1e}".format(alpha_0, beta_0))
    # print("!! setting is for the investigation for p1:p2 sync.:  {:g}*phi1 and {:g}*phi2".format(p1, p2))
    # print("")




    N = 2 # 振動子数


    # 平均速度(初期条件として与える) [rad/day]
    v1, v2 = [], []
    for phi1, phi2 in zip(phi1_list, phi2_list):
        v1.append(np.diff(phi1)/dt)
        v2.append(np.diff(phi2)/dt)
    v_1_mean = np.mean(np.concatenate(v1))
    v_2_mean = np.mean(np.concatenate(v2))


    ### 推定: receiver=1
    Estimator_1 = Estimator(receiver=0, N=2, p=p1, q=p2,
                                            M_min=M_min, M_max=M_max,
                                            log10_lambda_0_min=log10_lambda_0_min, 
                                            log10_lambda_0_max=log10_lambda_0_max, 
                                            alpha_0=alpha_0, beta_0=beta_0)
    # set intial condition (initial condition of omega)
    Estimator_1.set_init_all(v_1_mean)
    # update
    for phi1, phi2 in zip(phi1_list, phi2_list):
        # 推定に用いるデータを作成
        if phi1.size > 1: 
            phase_series = np.array([phi1,  phi2]).T
            Estimator_1.update_all(phase_series, delta_t=dt)
    # model section
    result_model_1, _ = Estimator_1.model_selection()


    ### 推定: receiver=2
    Estimator_2 = Estimator(receiver=1, N=2, p=p2, q=p1,
                                            M_min=M_min, M_max=M_max,
                                            log10_lambda_0_min=log10_lambda_0_min, 
                                            log10_lambda_0_max=log10_lambda_0_max,  
                                            alpha_0=alpha_0, beta_0=beta_0)
    # set intial condition (initial condition of omega)
    Estimator_2.set_init_all(v_2_mean)
    # update
    for phi1, phi2 in zip(phi1_list, phi2_list):
        # 推定に用いるデータを作成
        if phi1.size > 1: 
            phase_series = np.array([phi1,  phi2]).T
            Estimator_2.update_all(phase_series, delta_t=dt)
    # model section
    result_model_2, _ = Estimator_2.model_selection()


    if return_evidence:
        evidence = evidence_arrangement(Estimator_1, Estimator_2)
        return result_model_1, result_model_2, evidence
    return result_model_1, result_model_2




def data_arrangement(result_model_1, result_model_2, dt):
    """Collect selected-model parameters into arrays used by plotting and saving."""

    ##### データ整理
    Gamma_1_mean, Gamma_1_std, phi  = result_model_1.return_Gamma(sender=1)
    Omega_1_mean, Omega_1_std = result_model_1.return_omega()
    D_1_mean, D_1_std = result_model_1.return_D(delta_t=dt)
    sqrt2D_1_mean = np.sqrt(2.0*D_1_mean)  # 実際のノイズのSTD
    M_1 = result_model_1.return_Gamma_order(sender=1)
    #
    Gamma_2_mean, Gamma_2_std, _   = result_model_2.return_Gamma(sender=0)
    Omega_2_mean, Omega_2_std = result_model_2.return_omega()
    D_2_mean, D_2_std = result_model_2.return_D(delta_t=dt) # 2Dでノイズの分散
    sqrt2D_2_mean = np.sqrt(2.0*D_2_mean)  # 実際のノイズのSTD
    M_2 = result_model_2.return_Gamma_order(sender=0)
    #

    #* f(phi1 - phi2)の形式
    x = phi
    y1 = (Gamma_1_mean + Omega_1_mean)[::-1]
    y1_err = (np.sqrt(Gamma_1_std**2 + Omega_1_std**2))[::-1]
    M1 = M_1
    
    #* f(phi2 - phi1)の形式
    y2 = (Gamma_2_mean + Omega_2_mean)[::-1]
    y2_err = (np.sqrt(Gamma_2_std**2 + Omega_2_std**2))[::-1]
    M2 = M_2


    #* f(phi1 - phi2)の形式
    p = result_model_1.p
    q = result_model_1.q
    y3 = p*y1 - q*y2[::-1]
    y3_err = np.sqrt( (p*y1_err)**2 + (q*y2_err[::-1])**2 )
    sqrt2D_3_mean = np.sqrt( (p*sqrt2D_1_mean)**2 + (q*sqrt2D_2_mean)**2  )


    #### 選択されたモデルのパラメータ
    # print("選択されたモデルのパラメータ")
    # print("model1: lambda_0 = {:0.2e},  M = {:d}".format(result_model_1.lambda_0, M_1))
    # print("model2: lambda_0 = {:0.2e},  M = {:d}".format(result_model_2.lambda_0, M_2))
    # print("")

  #### 
    # print("推定された位相方程式のノイズ強度(sqrt(2D))")
    # print("model1: {:0.4e}".format(sqrt2D_1_mean))
    # print("model2: {:0.4e}".format(sqrt2D_2_mean))
    # print("odd (m:n): {:0.4e}".format(sqrt2D_3_mean))
    # print("")

    #### 
    # print("推定された位相方程式の振動数(sqrt(2D))")
    # print("model1: {:0.4e}".format(Omega_1_mean))
    # print("model2: {:0.4e}".format(Omega_2_mean))
    # print("")



    return {"x":x, 
                "y1":y1, "y1_err":y1_err,
                "y2":y2, "y2_err":y2_err,
                "y3":y3, "y3_err":y3_err,
                "M1":M1, "M2":M2,
                "D1": sqrt2D_1_mean,
                "D2": sqrt2D_2_mean,
                "D3": sqrt2D_3_mean,
                "omega1":Omega_1_mean,
                "omega2":Omega_2_mean,
                "pq": (result_model_1.p, result_model_1.q)
                }

    
