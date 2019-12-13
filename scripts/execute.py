# import packages
import numpy as np
import matplotlib.pyplot as plt
import elliptical_u_est as ellip
import demographics as demog
import SS
# import TPI

# model parameters
S = int(40)
E = 20
T = 4 * S
annual_beta = 0.9
beta = annual_beta ** (40 / S)
sigma = 3.0
chi = 1.0 * np.ones(S)
l_tilde = 1.0   # per period endowment for agent/ maximum labor supply
b_s = 0.6
upsilon = 1.5
theta = 0.9  # frisch elasticity of labor
# get estimated b and upsilon from ellipitical utility function
b_ellip, upsilon = ellip.estimation(theta, l_tilde)

# Firm
A = 1.0
annual_delta = 0.05
delta = 1 - ((1 - annual_delta) ** (40 / S))
alpha = 0.5

# Population
min_age = 1
max_age = 100
start_year = 2013
pop_graphs = False
# Get the immigration rates for every period
imm_rates_SS = demog.get_imm_resid(E+S, min_age, max_age, pop_graphs)

# Get population objects
(omega_path_S, g_n_SS, omega_SS, surv_rates_S, mort_rates_S, g_n_path,
 imm_rates_mat, omega_S_preTP) = demog.get_pop_objs(E, S, T, min_age,
                                                    max_age, start_year,
                                                    pop_graphs)

# Get moratality rates
totpers = E + S
(rho_s, infmort_rate) = demog.get_mort(totpers, min_age, max_age, graph=False)

# Solve the SS
r_init = 1 / beta - 1
xi = 0.1
ss_params = (beta, sigma, alpha, A, delta, xi, l_tilde,
             chi, theta, omega_SS, imm_rates_SS, rho_s, S)
r_ss, b_sp1_ss, euler_errors_ss = SS.solve_ss(r_init, ss_params)

# Economic growth
g_y = 0.02
