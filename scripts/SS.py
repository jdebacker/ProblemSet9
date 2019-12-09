import numpy as np
import scipy.optimize as opt
import os
import household as hh
import firms
import aggregates as aggr



def get_SS(params, graph=False):
    (KL_init, beta, sigma, chi_n, l_tilde, b, upsilon, S, alpha, A,
        delta) = params
    dist = 10
    mindist = 1e-08
    maxiter = 500
    ss_iter = 0
    xi = 0.2

    r_params = (alpha, A, delta)
    w_params = (alpha, A)

    while dist > mindist and ss_iter < maxiter:
        ss_iter += 1
        K, L = KL_init
        r = firms.get_r(K, L, r_params)
        w = firms.get_w(K, L, w_params)
        c1_guess = 0.5
        c1_params = (r, w, beta, sigma, chi_n, l_tilde, b, upsilon, S)
        results_c1 = opt.root(hh.get_b_sp1, c1_guess, params=(c1_params))
        c1 = results_c1.x
        c = hh.get_recurs_c(c1, r, beta, sigma, S)
        n = hh.get_n_s(c, w, sigma, chi_n, l_tilde, b,
                          upsilon)
        b = hh.get_recurs_b(c, n, r, w)
        K_new = aggr.get_K(b[:-1])
        L_new = aggr.get_L(n)
        KL_new = np.array([K_new, L_new])
        dist = ((KL_new - KL_init) ** 2).sum()
        KL_init = xi * KL_new + (1 - xi) * KL_init
        print('iter:', ss_iter, ' dist: ', dist)

    K_ss = K_new
    L_ss = L_new
    r_ss = r
    w_ss = w
    c_ss = c
    n_ss = n
    b_ss = b
    
    Y_params = (alpha, A)
    Y_ss = aggr.get_Y(K_ss, L_ss, Y_params)
    C_ss = aggr.get_C(c_ss)

    ss_output = {'c_ss': c_ss, 'n_ss': n_ss, 'b_ss': b_ss, 'K_ss': K_ss,
                 'L_ss': L_ss, 'r_ss': r_ss, 'w_ss': w_ss, 'Y_ss': Y_ss,
                 'C_ss': C_ss}

    return ss_output