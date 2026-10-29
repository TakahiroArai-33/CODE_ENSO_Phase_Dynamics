import numpy as np
import numpy.linalg as LA
import math
from copy import deepcopy, copy

# 位相振動子から, 振動数とPRCを求めるモジュール 

log10_lambda_0_min = 1
log10_lambda_0_max = 8
M_min = int(0)
M_max = int(10)
alpha_0 = 2.0
beta_0 = 2.0
lambda_0 = 2.0
TWOPI = 2.0 * np.pi







class Model():
    """
    receiver = 0 ~ N-1
    """
    def __init__(self, receiver, M_arr, p=1, q=1,
                     alpha_0=alpha_0, beta_0=beta_0, lambda_0=lambda_0):
        """Initialize one Bayesian phase-oscillator model for a receiver."""
        self.receiver = receiver
        self.M_arr = M_arr.astype(np.int32) # M_ijの列
        self.M_arr[receiver] = int(0)  # M_ii=0を強制する
        self.N = self.M_arr.size
        self.Nc = int(1 + 2*self.M_arr.sum())  # パラメータmu,Sigmaの次数
        self.grid = 2**10+1 # 推定結果を出力する際のphi, Delta_phiのグリット数
        self.p = p
        self.q = q

        def return_Sigma_0():
            """Build an order-independent isotropic prior covariance matrix."""
            return np.identity(self.Nc) / lambda_0
        
        def return_mu_0():
            """Build the zero prior mean vector for all model parameters."""
            return np.zeros([self.Nc,])

        self.T = int(0) # 推定に用いたデータのサンプリング数の累積をカウント
        # hyper parameter (パラメータの初期化にも使う)
        self.alpha_0 = alpha_0
        self.beta_0 = beta_0
        self.lambda_0 = lambda_0
        self.mu_0 = return_mu_0()
        self.Sigma_0 = return_Sigma_0()
        # set_param
        self.alpha = np.copy(self.alpha_0)
        self.beta = np.copy(self.beta_0)
        self.mu = np.copy(self.mu_0)
        self.Sigma = np.copy(self.Sigma_0)


    def make_F(self, phase):
        """ return design matrix """
        F = np.zeros([phase.shape[0]-1, self.Nc])
        # 位相速度 np.diff(phase) と同じ区間の中点位相で結合項を評価する。
        phase_ast = 0.5 * (phase[:-1, :] + phase[1:, :])
        # receiverの位相
        receiver = self.receiver
        receiver_phase = phase_ast[:, receiver]
        # omegaの更新に対応する領域
        F[:, 0] = 1
        # Gammaの更新に対応する領域
        col = int(1)
        for sender in range(self.N):
            # senderの位相
            sender_phase = phase_ast[:, sender]
            # 位相差 (q*\phi_s-p*\phi_r)
            Delta_phi = np.mod(self.q*sender_phase - self.p*receiver_phase, TWOPI)
            # 結合関数の次数
            M_ij = self.M_arr[sender] # dimention of Gamma_{ij} i:receiver, j:sender
            if sender == self.receiver:
                pass
            else:
                #　各orderごとに計算
                for m_1 in range(1, M_ij+1):
                    F[:, col] = np.cos(m_1 * Delta_phi)
                    F[:, col + 1] = np.sin(m_1 * Delta_phi)
                    col += 2
        return F


    def make_d(self, phase, delta_t):
        """ 学習データの従属変数(位相速度)を返す. """
        return np.diff(phase[:, self.receiver]) / delta_t


    def update(self, phase, delta_t):
        """ 
        update instance 
        input:
            phase: size[T(time), N(num of osci)]
            pert: size[T(time),]
            delta_t
        """
        F, d = self.make_F(phase), self.make_d(phase, delta_t)
        # set param (old)
        Sigma_old = np.copy(self.Sigma)
        mu_old = np.copy(self.mu)
        alpha_old = copy(self.alpha)
        beta_old = copy(self.beta)
        # update by Bayes theorem
        inv_Sigma_old = LA.inv(Sigma_old)
        T = d.size # Tは今回追加したデータ数, self.Tは累積
        self.Sigma = LA.inv(inv_Sigma_old  + F.T @ F)
        self.mu = self.Sigma @ (F.T @ d + inv_Sigma_old  @ mu_old)
        self.alpha = alpha_old + (T / 2.0)
        self.beta = beta_old + 0.5 * (d @ d + mu_old @ inv_Sigma_old @ mu_old - self.mu @ LA.inv(self.Sigma) @ self.mu)
        self.T += T
    
    
    def return_Gamma(self, sender):
        """
        return mean and std of coupling function between receiver and sender
        """
        def make_G(phi, M):
            """Build the Fourier basis vector for one coupling-function phase."""
            G = np.zeros([2*M, ])
            G[0::2] = np.cos(np.arange(1, M+1)*phi)
            G[1::2] = np.sin(np.arange(1, M+1)*phi)
            return G
        
        # preparation        
        Gamma_mean, Gamma_var = np.zeros([self.grid, ]), np.zeros([self.grid, ])
        phi = np.linspace(0, 2.0*np.pi, num=self.grid)

        if sender != self.receiver:
            M = self.M_arr[sender]
            marker = int(1+2*self.M_arr[:sender].sum())
            # variance and mean of Gamma
            sub_mu = self.mu[marker:(marker+2*M)]
            sub_Sigma = self.Sigma[marker:(marker+2*M), marker:(marker+2*M)]
            for i_1 in range(phi.size):
                G = make_G(phi[i_1], M)
                Gamma_mean[i_1] = G @ sub_mu
                Gamma_var[i_1] = G @ sub_Sigma @ G
            Gamma_var = (self.beta / (self.alpha - 1.0)) * Gamma_var
        elif sender == self.receiver:
            pass
        return Gamma_mean, np.sqrt(Gamma_var), phi
    
    
    def return_Gamma_order(self, sender):
        """Return the selected Fourier order for coupling from ``sender``."""
        return self.M_arr[sender]
    

    def return_omega(self):
        """
        return mean and variance of constant term
        """
        omega_mean = self.mu[0]
        omega_var =  (self.beta / (self.alpha - 1.0)) * self.Sigma[0, 0]
        return omega_mean, np.sqrt(omega_var)

    
    def return_D(self, delta_t):
        """
        return strength of noise
        """
        D_hat_mean = (self.beta / (self.alpha - 1.0))
        D_hat_var = self.beta**2 / ((self.alpha - 1.0)**2 * (self.alpha - 2.0))
        D_mean, D_var = D_hat_mean * (delta_t / 2.0), D_hat_var * (delta_t / 2.0)**2
        return D_mean, np.sqrt(D_var)


    def evidence(self):
        """Return the log10 marginal likelihood used for model selection."""
        Sigma = self.Sigma
        alpha = self.alpha
        beta = self.beta
        Sigma_0 = self.Sigma_0
        alpha_0 = self.alpha_0
        beta_0 = self.beta_0
        T = self.T

        # 共分散行列は半正定値
        log10det_Sigma = LA.slogdet(Sigma)[1] * np.log10(np.e) # log10[Sigma]
        log10det_Sigma_0 = LA.slogdet(Sigma_0)[1] * np.log10(np.e) # log10[Sigma_0]            
        f_1 = - (T / 2.0) * np.log10(2.0 * np.pi)
        f_2 = log10det_Sigma / 2.0
        f_3 = - log10det_Sigma_0 / 2.0
        f_4 = alpha_0 * np.log10(beta_0)
        f_5 = - alpha * np.log10(beta)
        f_6 = math.lgamma(alpha) / np.log(10.0)
        f_7 = - math.lgamma(alpha_0) / np.log(10.0)
        return f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7



class Estimator():
    """Manage candidate phase models and select the model with highest evidence."""

    def __init__(self, receiver, N, p=int(1), q=int(1),
                  M_min=M_min, M_max=M_max,
                  log10_lambda_0_min=log10_lambda_0_min,
                  log10_lambda_0_max=log10_lambda_0_max,
                  alpha_0=alpha_0, beta_0=beta_0):
        """
        モデル選択の候補となるインスタンスのリストを生成. 
        各インスタンスは異なる結合関数の次数やlambda_0に対応している.
        p,qはm:n結合に関するパラメータで
        Gamma(q\phi_s-p\phi_r)となる.
        """
        def update_M_arr(M_arr):
            """Advance the candidate Fourier-order vector across all senders."""
            M_arr[0] += 1
            for i_1 in range(M_arr.size - 1):
                if M_arr[i_1] > M_max:
                    M_arr[i_1] = M_min
                    M_arr[i_1+1] += 1
            if M_arr[-1] > M_max:
                return None
            else:
                return M_arr.astype(int)
            
        # クラス Model のインスタンス列
        self.model_list = []
        # 結合関数の次数の列 {M_ij}_j ただし, M_ii は除く 
        M_arr_excluding_receiver = np.array([M_min for _ in range(N-1)], dtype=np.int32)
        # iteration (Gamma)
        while True:
            # iteration (lambda) 
            for log10_lambda_0 in range(int(log10_lambda_0_min), int(log10_lambda_0_max+1.0)):
                lambda_0 = 10.0**log10_lambda_0
                # Modelインスタンス作成,　インスタンス列に追加
                M_arr = np.insert(M_arr_excluding_receiver, 
                                  receiver, int(0)) # self-coupling (次数0) を含める
                self.model_list.append(Model(receiver, M_arr, p, q,
                                              alpha_0, beta_0, lambda_0))
            # update M_arr
            M_arr_excluding_receiver = update_M_arr(M_arr_excluding_receiver)
            if M_arr_excluding_receiver == None:
                break
            
            
    def set_init_all(self, value):
        """
        モデル全てに初期条件を与える, 
        パラメータ self.mu, self.mu_0 の第0要素(omegaに対応)に指定の値を与える. (平均角速度など)
        """
        for model in self.model_list:
            model.mu_0[0] = value
            model.mu[0] = value


    def update_all(self, phase, delta_t):
        """
        モデルをすべて更新したものを返す. 
        input: phase, phase_data
        """
        for model in self.model_list:
            model.update(phase, delta_t)
            

    def model_selection(self):
        """
        モデル選択: エビデンス近似によるモデル選択
        """
        L_list = np.array([model.evidence() for model in self.model_list])
        L = L_list.max()
        model = self.model_list[np.argmax(L_list)] 
        return model, L


    def all_evidence(self):
        """
        エビデンス近似によるモデル選択のために周辺尤度(エビデンス)の値をインスタンスごとにを返す
        """

        # L_listは各instanceのevidence, instance_list(input)からL_listのargmaxをとれば, model selectionとなる.
        L_list = np.array([model.evidence() for model in self.model_list])
        return self.model_list, L_list


     
