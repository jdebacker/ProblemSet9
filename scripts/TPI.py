import numpy as np
import scipy.optimize as opt
import firm
import households as hh
import aggregates as agg


def solve_tp(r_path_init, w_path_init, params):
    '''
    Solves for the time path equlibrium using TPI
    '''
    (beta, sigma, alpha, A, delta, theta, T, xi, b_sp1_pre, r_ss, w_ss,
     b_sp1_ss) = params
    tpi_dist = 7.0
    tpi_tol = 1e-8
    tpi_iter = 0
    tpi_max_iter = 300
    r_path = np.append(r_path_init, np.ones(S-1) * r_ss)
    w_path = np.append(w_path_init, np.ones(S) * b_ss)
    while (tpi_dist > tpi_tol) & (tpi_iter < tpi_max_iter):
        w_path = firm.get_w(L, K, alpha, A, delta)
        # Solve HH problem
        b_sp1_mat = np.zeros((T + S - 1, S))
        euler_errors_mat = np.zeros((T + S - 1, S))
        
