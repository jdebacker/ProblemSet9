# market clearning conditions


import numpy as np


def get_L(n, omega_SS):
    L = (omega_SS * n).sum()

    return L


def get_C(c, omega_SS):
    C = (omega_SS * c).sum()

    return C


def get_Y(K, L, params):
    alpha, A = params
    Y = A * (K ** alpha) * (L ** (1 - alpha))

    return Y


def get_I(K, Kp1, b, params):
    (delta, g_n_SS, omega_SS, imm_rates_SS, g_y) = params
    capital_flow = ((1 + g_n_SS) * np.exp(g_y) * (imm_rates_SS[1:] *
                    omega_SS[1:] * b).sum())

    I = ((1 + g_n_SS) * np.exp(g_y) * Kp1 - (1.0 - delta) * K -
         capital_flow)

    return I


def get_K(b, params):
    '''
    b = (S-1,) vector of savings
    '''
    (g_n_SS, omega_SS, imm_rates_SS) = params
    K = ((1 / (1 + g_n_SS)) * (omega_SS[:-1] * b + imm_rates_SS[1:] *
                               omega_SS[1:] * b).sum())

    return K



def get_BQ(b, r, params):
    (g_n_SS, omega_SS, rho_s) = params
    BQ = (((1 + r) / (1 + g_n_SS)) * (rho_s[:-1] * omega_SS[:-1] *
                                      b).sum())

    return BQ


def get_r(K, L, params):
    alpha, A, delta = params
    r = alpha * A * ((L / K) ** (1 - alpha)) - delta

    return r


def get_w(K, L, params):
    alpha, A = params
    w = (1 - alpha) * A * ((K / L) ** alpha)

    return w

