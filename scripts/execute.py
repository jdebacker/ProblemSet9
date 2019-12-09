# import packages
import numpy as np
# import matplotlib.pyplot as plt
import SS
import TPI
import demographics as demog
import elliptical_u_est

# model parameters
S = int(40)
E = 20
T = 4 * S
annual_beta = 0.9
beta = annual_beta ** (40 / S)
sigma = 3.0
chi_n = 1.0 * np.ones(S)
l_tilde = 1.0   # per period endowment for agent/ maximum labor supply
b_s = 0.6
upsilon = 1.5
theta = 0.9  # frisch elasticity of labor
# get estimated b and upsilon from ellipitical utility function
b_ellip, upsilon = elliptical_u_est.estimation(theta, l_tilde)

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
(omega_path_S, imm_rates_path, rho_s, omega_SS, surv_rates_S, g_n_path,
    g_n_SS, omega_S_preTP) = demog.get_pop_objs(E, S, T, min_age, max_age,
                                                start_year, pop_graphs)
imm_rates_SS = imm_rates_path[-1, :]

# Economic growth
g_y = 0.02
