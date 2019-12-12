

import numpy as np
import elliptical_u_est as ellip

# household functions

def HHFOCs(b_sp1, *params):
    '''
    Refer to the equations on chapter 8
    '''
    BQpath, rho_s, beta, sigma, b_init, n, rpath, wpath = params
    # rpath = length p, n = length p, 
    # p = number of periods left in life time, 
    #wpath = length p, b_sp1 = length p-1
    
    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)
    c = (1 + rpath) * b_s + wpath * n + BQpath - b_sp1
    MU_c = c ** (-sigma)
    errors = (MU_c[:-1] - beta * (1 + rpath[1:]) *
              (1 - rho_s[:-1]) * MU_c[1:])
    
    return errors




def mu_labor(n_s, l_tilde, chi, theta):
    '''
    Computes marginal utility with respect to labor supply
    '''

    b_ellipse, nu = ellip.estimation(theta, l_tilde)
    # b_ellipse is the constant defined in Equation (4.9) - has nothing to do with savings
    mu_n = chi * (b_ellipse / l_tilde) * (n_s / l_tilde) ** (nu - 1) * (
                                 1 - (n_s / l_tilde) ** nu) ** ((1 - nu) / nu )

    return mu_n
