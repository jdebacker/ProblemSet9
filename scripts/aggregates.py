# market clearning conditions


import numpy as np


def get_L(n_s, omega_SS):
    '''
    Returns labor (a scalar)
    '''
    L = (omega_SS[1:] * n_s[1:]).sum() # indexing from 1 since the summation starts at E + 1, not E

    return L

def get_K(b_s, *params):
    '''
    b = vector of savings of length S - 1

    Returns capital (a scalar)
    '''

    (g_n_SS, omega_SS, imm_rates_SS) = params
    K = ((1 / (1 + g_n_SS)) * (omega_SS[:-1] * b_s + imm_rates_SS[1:] *
                               omega_SS[1:] * b_s).sum())

    return K


def get_C(c, omega_SS):
    C = (omega_SS[1:] * c[1:]).sum() # indexing from 1 since the summation starts at E + 1, not E

    return C

def get_I(K, K_sp1, delta, g_n_SS):
    # Need to think about how we're getting Kp1 - something analogous to what we did for b_sp_1 in class I presume
    I = ((1 + g_n_SS) * K_sp1) - ((1 - delta) * K)
    return I


def get_Y(b_s, C, I, omega_SS, imm_rates_SS):
    '''
    Returns total goods consumed (a scalar)
    '''

    Y = C + I - (imm_rates_SS[1:] * omega_SS[1:] * b_s).sum()

    return Y


def get_BQ(b_s, r, params):
    (omega_SS, rho_s) = params
    BQ = (1 + r) * (rho_s[:-1] * omega_SS[:-1] * b_s).sum()

    return BQ


