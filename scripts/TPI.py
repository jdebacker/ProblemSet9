import numpy as np
import scipy.optimize as opt
import aggregates as agg
import households as hh
import firm


def solve_tp(g_n_path, omega_S_preTP, rho_s, imm_rates_path, params):
    '''
    Solves for the time path equilibrium using TPI
    '''
    # Missing some elements of params
    b_ss, r_11, T, S = params
    dist = 8.0
    mindist = 1e-08
    maxiter = 300
    tpi_iter = 0
    xi = 0.2
    while dist > mindist and tpi_iter < maxiter:

        # Define paths
        b_11 = 1.1 * b_ss
        BQ_params = (g_n_path[0], omega_S_preTP, rho_s)
        K_params = (g_n_path[0], omega_S_preTP, imm_rates_path[0, :])
        BQ_11 = agg.get_BQ(b_11, r_11, BQ_params)
        K_11 = agg.get_K(b_11, K_params)
        BQpath_init = np.zeros(T + S - 1)
        BQpath_init[:T] = np.linspace(BQ_11, BQ_ss, T)
        BQpath_init[T:] = BQ_ss
        r_11 = firm.get_r(K_11, L_ss, r_params)
        r_path = firm.get_r(L_ss, r_params)
        w_path = firm.get_w(L_ss, w_params)
        bmat = np.zeros((S - 1, T + S - 1))
        bmat[:, 0] = b_1

        # Solve for households
        for p in range(2, S):
            b_guess = np.diagonal(bmat[S - p:, :p - 1])
            b_init = bmat[S - p - 1, 0]
            b_params = (b_init, n[-p:], r_path[:p], w_path[:p],
                        BQpath_init[:p], rho_s[-p:], beta, sigma)
            results_bp = opt.root(hh.FOCs, b_guess, args=(b_params))
            b_solve_p = results_bp.x
            DiagMaskbp = np.eye(p - 1, dtype=bool)
            bmat[S - p:, 1:p] = DiagMaskbp * b_solve_p + bmat[S - p:, 1:p]

        for t in range(1, T + 1):
            b_guess = np.diagonal(bmat[:, t - 1:t + S - 2])
            b_init = 0.0
            b_params = (b_init, n, r_path[t - 1:t + S - 1],
                        w_path[t - 1:t + S - 1], BQpath_init[t - 1:t + S - 1],
                        rho_s, beta, sigma)
            results_bt = opt.root(hh.FOCs, b_guess, args=(b_params))
            b_solve_t = results_bt.x
            DiagMaskbt = np.eye(S - 1, dtype=bool)
            bmat[:, t:t + S - 1] = (DiagMaskbt * b_solve_t +
                                    bmat[:, t:t + S - 1])

        new_Kpath = np.zeros(T)
        new_Kpath[0] = K_1
        new_Kpath[1:] = \
            (1 / (omega_path_S[:T - 1, :-1]) *
                bmat[:, 1:T].T +
                imm_rates_path[:T - 1, 1:] *
                omega_path_S[:T - 1, 1:] * bmat[:, 1:T].T).sum(axis=1)
        new_BQpath = np.zeros(T)
        new_BQpath[0] = BQ_1
        new_BQpath[1:] = \
            ((1 + r_path[1:T]) / (rho_s[:-1]) *
                omega_path_S[:T - 1, :-1] * bmat[:, 1:T].T).sum(axis=1)

        dist = ((BQ_init - new_BQ) ** 2).sum()
        BQpath_init[:T] = xi * new_BQpath[:T] + (1 - xi) * BQpath_init[:T]
        # update iteration counter
        tpi_iter += 1

    if tpi_iter < maxiter:
        print('The time path solved! ->', ' iter:', tpi_iter, ', dist: ', dist)
    else:
        print('The time path did not solve.')

    return [new_Kpath, new_BQpath]
