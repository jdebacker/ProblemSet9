# Import packages
import numpy as np
import scipy.optimize as opt


def HHFOCs(b, *params):
    b_init, n, rpath, wpath, BQpath, rho_s, beta, sigma, g_y = params
    b_s = np.append(b_init, b)
    b_sp1 = np.append(b, 0.0)
    c = (1 + rpath) * b_s + wpath * n + BQpath - np.exp(g_y) * b_sp1
    MU_c = c ** (-sigma)
    errors = (MU_c[:-1] - np.exp(-sigma * g_y) * beta * (1 + rpath[1:]) *
              (1 - rho_s[:-1]) * MU_c[1:])
    # rpath = length p
    # n = length p
    # wpath = length p
    # p = number of periods left in life time
    # b = length p-1
    return errors

# TPI solution
b1 = 1.1 * b_ss
K_params = (g_n_path[0], omega_S_preTP, imm_rates_path[0, :])
K1 = get_K(b1, K_params)
Kpath_init = np.zeros(T + S - 1)
Kpath_init[:T] = np.linspace(K1, K_ss, T)
Kpath_init[T:] = K_ss
BQ_params = (g_n_path[0], omega_S_preTP, rho_s)
BQ1 = get_BQ(b1, r1, BQ_params)
BQpath_init = np.zeros(T + S - 1)
BQpath_init[:T] = np.linspace(BQ1, BQ_ss, T)
BQpath_init[T:] = BQ_ss
r1 = get_r(K1, L_ss, r_params)

dist = 7.0
mindist = 1e-08
maxiter = 300
tpi_iter = 0
xi = 0.2

while dist > mindist and tpi_iter < maxiter:
    tpi_iter += 1
    # Get r and w paths
    rpath = get_r(Kpath_init, L_ss, r_params)
    wpath = get_w(Kpath_init, L_ss, w_params)
    bmat = np.zeros((S - 1, T + S - 1))
    bmat[:, 0] = b1
    # Solve for households' problems
    for p in range(2, S):
        b_guess = np.diagonal(bmat[S - p:, :p - 1])
        b_init = bmat[S - p - 1, 0]
        b_params = (b_init, n[-p:], rpath[:p], wpath[:p],
                  BQpath_init[:p], rho_s[-p:], beta, sigma, g_y)
        results_bp = opt.root(HHFOCs, b_guess, args=(b_params))
        b_sol_p = results_bp.x
        DiagMaskbp = np.eye(p - 1, dtype=bool)
        bmat[S - p:, 1:p] = DiagMaskbp * b_sol_p + bmat[S - p:, 1:p]

    for t in range(1, T + 1):  # Go from periods 1 to S+T-1
        b_guess = np.diagonal(bmat[:, t - 1:t + S - 2])
        b_init = 0.0
        b_params = (b_init, n, rpath[t - 1:t + S - 1],
                  wpath[t - 1:t + S - 1], BQpath_init[t - 1:t + S - 1],
                  rho_s, beta, sigma, g_y)
        results_bt = opt.root(HHFOCs, b_guess, args=(b_params))
        b_sol_t = results_bt.x
        DiagMaskbt = np.eye(S - 1, dtype=bool)
        bmat[:, t:t + S - 1] = (DiagMaskbt * b_sol_t +
                                bmat[:, t:t + S - 1])

    Kpath_new = np.zeros(T)
    BQpath_new = np.zeros(T)
    Kpath_new[0] = K1
    BQpath_new[0] = BQ1
    Kpath_new[1:] = \
        (1 / (1 + g_n_path[1:T])) * (omega_path_S[:T - 1, :-1] *
         bmat[:, 1:T].T + imm_rates_path[:T - 1, 1:] *
         omega_path_S[:T - 1, 1:] * bmat[:, 1:T].T).sum(axis=1)
    BQpath_new[1:] = \
        ((1 + rpath[1:T]) / (1 + g_n_path[1:T])) * (rho_s[:-1] *
         omega_path_S[:T - 1, :-1] *
         bmat[:, 1:T].T).sum(axis=1)
    KBQ_init = np.append(Kpath_init[:T], BQpath_init[:T])
    KBQ_new = np.append(Kpath_new[:T], BQpath_new[:T])
    dist = ((KBQ_init - KBQ_new) ** 2).sum()
    Kpath_init[:T] = xi * Kpath_new[:T] + (1 - xi) * Kpath_init[:T]
    BQpath_init[:T] = xi * BQpath_new[:T] + (1 - xi) * BQpath_init[:T]
    
    print('iter:', tpi_iter, ' dist: ', dist)
