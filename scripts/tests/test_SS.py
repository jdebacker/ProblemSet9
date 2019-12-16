import numpy as np
import demographics
import SS


def test_get_pop_objs():
    """
    Test of the that omega_SS and the last period of omega_path_S are
    close to each other.
    """
    E = 20
    S = 80
    T = int(round(4.0 * S))
    start_year = 2018

    (omega, g_n_ss, omega_SS, surv_rate, rho, g_n_vector, imm_rates,
        omega_S_preTP) = demographics.get_pop_objs(E, S, T, 1, 100,
                                                   start_year, False)
    return omega_SS,imm_rates

omega_SS,imm_rates = test_get_pop_objs()

def test_SS():
    """
    Test whether the SS funtion works well and what the returns look like
    """
    r_init = 0.1
    beta = 0.8
    sigma = 1.5
    alpha = 0.3
    A = 1.0
    delta = 0.1
    xi = 0.1
    S = 80
    imm_rates_SS = imm_rates

    (r,w,b_sp1,euler_errors) = SS.solve_ss(r_init,params=(beta, sigma, alpha,
                                                          A, delta, xi, omega_SS,
                                                          imm_rates_SS, S))
    print("r = ",r)
    print("w = ",w)
    print("b_sp1 = ",b_sp1)
    print("euler_errors = ",euler_errors)

test_SS()
