from scipy import optimize as opt
import numpy as np
import firm
import households as hh
import aggregates as agg


def solve_ss(r_init, w_init, params):
    '''
    Solves for the steady-state equlibrium of the OG model
    '''
    beta, sigma, n, alpha, A, delta, xi = params
    ss_dist = 7.0
    ss_tol = 1e-8
    ss_iter = 0
    ss_max_iter = 300
    r = r_init
    w = w_init
    while (ss_dist > ss_tol) & (ss_iter < ss_max_iter):
        # solve HH problem
        foc_args = (beta, sigma, r, w, 0.0)
        n_s_guess = np.ones(S)
        b_sp1_guess = np.ones(S-1) * 0.5
        HH_guess = np.append(b_sp1_guess, n_s_guess)
        result = opt.root(hh.FOCs, HH_guess, args=foc_args)
        b_sp1 = result.x[0: S-1]
        n_s = result.x[S-1:]
        euler_errors = result.fun
        b_s = np.append(0.0, b_sp1)
        # use market clearing
        L = agg.get_L(n)
        K = agg.get_K(b_s)
        # find implied r
        r_prime = firm.get_r(L, K, alpha, A, delta)
        # find implied w
        w_prime = firm.get_w(L, K, G, alpha, A)
        # check distance
        ss_dist_r = np.absolute(r - r_prime)
        ss_dist_w = np.absolute(w - w_prime)
        print('Iteration = ', ss_iter, ', Distance r = ', ss_dist_r,
              ', Disrance w = ' ss_dist_w,
              ', r = ', r, ', w = ', w)
        # update r
        r = xi * r_prime + (1 - xi) * r
        # update w
        w = xi * w_prime + (1 - xi) * w
        # update iteration counter
        ss_iter += 1

    return r, w, b_sp1, euler_errors
