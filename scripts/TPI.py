import numpy as np
import scipy.optimize as opt
import firm
import households as hh
import aggregates as agg


def solve_tp(r_path_init, BQ_path_init, params):
    '''
    Solves for the time path equlibrium using TPI
    '''
    (beta, sigma, alpha, A, delta, theta, T, xi, b_sp1_pre, r_ss, BQ_ss, b_sp1_ss) = params
    tpi_dist = 7.0
    tpi_tol = 1e-8
    tpi_iter = 0
    tpi_max_iter = 300
    r_path = np.append(r_path_init, np.ones(S) * r_ss)
    BQ_path = np.append(BQ_path_init, np.ones(S) * BQ_ss)
    while (tpi_dist > tpi_tol) & (tpi_iter < tpi_max_iter):
        w_path = firm.get_w(r_path, alpha, A, delta)
        # Solve HH problem
        b_sp1_mat = np.zeros((T + S, S))
        euler_errors_mat = np.zeros((T + S, S))
        # solve upper right elements before the first full lifetime
        UpMaskb = np.triu(np.ones((T + S, S)), 1)
        for t in range(T + S):
            foc_args = (beta, sigma, r_path[t:t+S], w_path[t:t+S],    0.0, l_tilde, chi, theta, rho_s)
            b_sp1_guess = b_sp1_ss
            result = opt.root(hh.FOCs, b_sp1_guess, args=foc_args)
            b_sp1_mat[t:t+S, :] = (UpMaskb * result.x +
                                   b_sp1_mat[t:t+S, :])
            euler_errors_mat[t:t+S, :] = (UpMaskb * result.fun + euler_errors_mat[t:t+S, :])
        # solve all full lifetimes
        DiagMaskb = np.eye(S, dtype=bool)
        for t in range(T + S):
            foc_args = (beta, sigma, r_path[t:t+S], w_path[t:t+S],    0.0, l_tilde, chi, theta, rho_s)
            b_sp1_guess = b_sp1_ss
            result = opt.root(hh.FOCs, b_sp1_guess, args=foc_args)
            b_sp1_mat[t:t+S, :] = (DiagMaskb * result.x +
                                   b_sp1_mat[t:t+S, :])
            euler_errors_mat[t:t+S, :] = (DiagMaskb * result.fun + euler_errors_mat[t:t+S, :])
        # create a b_s_mat
        b_s_mat = np.zeros((T + S, S))
        b_s_mat[]
