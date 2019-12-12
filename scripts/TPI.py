import numpy as np
import scipy.optimize as opt


while dist > mindist and tpi_iter < maxiter:
    
    # define paths
    b_11 = 1.1 * b_ss
    BQ_params = (g_n_path[0], omega_S_preTP, rho_s)
    K_params = (g_n_path[0], omega_S_preTP, imm_rates_path[0, :])
    BQ_11 = get_BQ(b_11, r_11, BQ_params)
    K_11 = get_K(b_11, K_params)
    BQpath_init = np.zeros(T + S - 1)
    Kpath_init = np.zeros(T + S - 1)
    BQpath_init[:T] = np.linspace(BQ_11, BQ_ss, T)
    Kpath_init[:T] = np.linspace(K_11, K_ss, T)
    BQpath_init[T:] = BQ_ss
    Kpath_init[T:] = K_ss
    r_11 = get_r(K_11, L_ss, r_params)
    rpath = get_r(Kpath_init, L_ss, r_params)
    wpath = get_w(Kpath_init, L_ss, w_params)
    bmat = np.zeros((S - 1, T + S - 1))
    bmat[:, 0] = b_1
    
    tpi_iter += 1
    dist = 8.0
    mindist = 1e-08
    maxiter = 300
    tpi_iter = 0
    xi = 0.2

    # Solve for households
    for p in range(2, S):
        b_guess = np.diagonal(bmat[S - p:, :p - 1])
        b_init = bmat[S - p - 1, 0]
        b_params = (b_init, n[-p:], rpath[:p], wpath[:p],
                  BQpath_init[:p], rho_s[-p:], beta, sigma)
        results_bp = opt.root(HHFOCs, b_guess, args=(b_params))
        b_solve_p = results_bp.x
        DiagMaskbp = np.eye(p - 1, dtype=bool)
        bmat[S - p:, 1:p] = DiagMaskbp * b_solve_p + bmat[S - p:, 1:p]

    for t in range(1, T + 1): 
        b_guess = np.diagonal(bmat[:, t - 1:t + S - 2])
        b_init = 0.0
        b_params = (b_init, n, rpath[t - 1:t + S - 1],
                  wpath[t - 1:t + S - 1], BQpath_init[t - 1:t + S - 1],
                  rho_s, beta, sigma)
        results_bt = opt.root(HHFOCs, b_guess, args=(b_params))
        b_solve_t = results_bt.x
        DiagMaskbt = np.eye(S - 1, dtype=bool)
        bmat[:, t:t + S - 1] = (DiagMaskbt * b_solve_t +
                                bmat[:, t:t + S - 1])

    new_Kpath= np.zeros(T)
    new_Kpath[0] = K_1
    new_Kpath[1:] = \
        (1 / (omega_path_S[:T - 1, :-1] *
         bmat[:, 1:T].T + imm_rates_path[:T - 1, 1:] *
         omega_path_S[:T - 1, 1:] * bmat[:, 1:T].T).sum(axis=1)
    new_BQpath = np.zeros(T)
    new_BQpath[0] = BQ_1
    new_BQpath[1:] = \
        ((1 + rpath[1:T]) / (rho_s[:-1] *
         omega_path_S[:T - 1, :-1] *
         bmat[:, 1:T].T).sum(axis=1)
    
    KBQ_init = np.append(Kpath_init[:T], BQpath_init[:T])
    Kpath_init[:T] = xi * new_Kpath[:T] + (1 - xi) * Kpath_init[:T]
    new_KBQ = np.append(new_Kpath[:T], new_BQpath[:T])
    dist = ((KBQ_init - new_KBQ) ** 2).sum()
    BQpath_init[:T] = xi * new_BQpath[:T] + (1 - xi) * BQpath_init[:T]
    
    print('iter:', tpi_iter, ' dist: ', dist)
