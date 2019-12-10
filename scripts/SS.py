from scipy import optimize as opt
import numpy as np
import firm
import aggregates as agg
import households as hh


def solve_ss(r_init, w_init, params):
    '''
    Solves for the steady-state equlibrium of the OG model
    '''
    beta, sigma, alpha, A, delta, xi, omega_SS, imm_rates_SS, S = params
    ss_dist = 7.0
    ss_tol = 1e-8
    ss_iter = 0
    ss_max_iter = 300
    r = r_init
    # w = w_init
    # Why do we need w as well? Are we not going to get it by providing r_init to firm.get_w()?
    # I think Jason said this too and I'm going to stick to what we did in class.
    while (ss_dist > ss_tol) & (ss_iter < ss_max_iter):
        # get w
        w = firm.get_w(r, alpha, A, delta)
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
        L = agg.get_L(n_s, omega_SS)
        K = agg.get_K(b_s, omega_SS, imm_rates_SS)
        # find implied r
        r_prime = firm.get_r(L, K, alpha, A, delta)
        # find implied w
        w_prime = firm.get_w(r, alpha, A, delta)
        # check distance
        ss_dist_r = np.absolute(r - r_prime)
        ss_dist_w = np.absolute(w - w_prime)
        print('Iteration = ', ss_iter, ', Distance r = ', ss_dist_r,
              ', Disrance w = ', ss_dist_w,
              ', r = ', r, ', w = ', w)
        # update r
        r = xi * r_prime + (1 - xi) * r
        # update w
        w = xi * w_prime + (1 - xi) * w
        # update iteration counter
        ss_iter += 1

    return r, w, b_sp1, euler_errors
