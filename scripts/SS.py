import numpy as np
import scipy.optimize as opt


def EuFunc(b, *params):
    (n, beta, sigma, alpha, A, delta, g_n_SS, omega_SS, imm_rates_SS, rho_s, g_y) = params
    b_s = np.append(0.0, b)
    b_sp1 = np.append(b, 0.0)
    c = (1 + r) * b_s + w * n - np.exp(g_y) * b_sp1 + BQ
    MU_c = c**(-sigma)
    errors = (MU_c[:-1] - (np.exp(-sigma * g_y) * (1 - rho_s[:-1]) *
                           beta * (1 + r) * MU_c[1:]))

    return errors

# SS solutions

b_init = np.ones(S - 1) * 0.1
b_params = (g_n_SS, omega_SS, n, beta, sigma, alpha, A, delta, imm_rates_SS, rho_s, g_y)
results_b = opt.root(EuFunc, b_init, args=(b_params))
print(results_b)
b_ss = results_b.x
print('SS savings: ', b_ss)
K_params = (g_n_SS, omega_SS, imm_rates_SS)
K_ss = get_K(b_ss, K_params)
L_ss = get_L(n, omega_SS)
print('K_ss and L_ss', np.array([K_ss, L_ss]))
r_params = (alpha, A, delta)
w_params = (alpha, A)
r_ss = get_r(K_ss, L_ss, r_params)
w_ss = get_w(K_ss, L_ss, w_params)
print('r and w in SS: ', np.array([r_ss, w_ss]))
b_s_ss = np.append(0.0, b_ss)
b_sp1_ss = np.append(b_ss, 0.0)
BQ_params = (g_n_SS, omega_SS, rho_s)
BQ_ss = get_BQ(b_ss, r_ss, BQ_params)
print('BQ SS: ', BQ_ss)
c_ss = (1 + r_ss) * b_s_ss + w_ss * n - np.exp(g_y) * b_sp1_ss + BQ_ss
print('Consumption SS: ', c_ss)
Y_params = (alpha, A)
Y_ss = get_Y(K_ss, L_ss, Y_params)
C_ss = get_C(c_ss, omega_SS)
I_params = (delta, g_n_SS, omega_SS, imm_rates_SS, g_y)
I_ss = get_I(K_ss, K_ss, b_ss, I_params)
RC_error = Y_ss - C_ss - I_ss
print('RC Error: ', RC_error)

