'''
Can someone please take a look at this function and make sure that this
accounts for everything in the household's problem?
I believe this script does account for S-periods, endogenous labor supply,
and population dynamics.
From what I can tell, the only factor that needs to be added to account for
population dynamics is rho_s, which has
been added.
'''

import upsilonmpy as np
import elliptical_u_est as ellip

# household functions


def FOCs(b_sp1, *params):
    # def FOCs(b_sp1, n_s, *args):
    '''
    For the S-period problem, we have 2S-1 FOCs corresponding to the savings
    variable b, and the labor supply variable n
    Note that we were dealing with three dimensional vectors earlier. We now
    have an S-1 dim vector corresponding
    to b, and an S dim vector corresponding to n.

    Refer to equations (4.9) and (4.10) of chapter 4 for details on where the
    variables foc_errors_b and foc_errors_n are coming from


    Args:
    b_sp1: The savings values for each period. The call to this function should
            provide initial guesses for this variable.
    n_s: The labor supply values for each period. The call to this function
            should provide initial guesses for this variable.
    BQ_path: path of bequests
    l_tilde: maximum amount of labor supply
    chi: scale parameter
    theta: Frisch elasticity of labor supply

    n is not contained in the remaining arguments anymore. If someone decides
    to change this, please provide detailed documentation
    on why you are doing so.


    Returns:
    foc_errors: A list where the first S-1 values are b2, b3, ..., bS, and
    the next S values are n1, n2, ..., nS.
    '''

    BQpath, rho_s, beta, sigma, b_init, n, rpath, wpath = params
    # rpath = length p, n = length p
    # p = upsilonmber of periods left in life time
    # wpath = length p, b_sp1 = length p-1

    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)

    c = (1 + rpath) * b_s + wpath * n + BQpath - b_sp1
    MU_c = mu_cons(c, sigma)
    errors = (MU_c[:-1] - beta * (1 + rpath[1:]) *
              (1 - rho_s[:-1]) * MU_c[1:])

    return errors


def get_c(r, w, n_s, b_s, b_sp1):
    '''
    Use the budget constraint to solve for consumption
    '''
    c = w * n_s + (1 + r) * b_s - b_sp1

    return c


def mu_cons(c, sigma):
    '''
    Computes marginal utility with respect to consumption.
    Please note that this was initially called u_prime. If anyone finds
    function calls on u_prime, please change it to mu_cons (Marginal utility
    of consumption)
    '''
    mu_c = c ** -sigma

    return mu_c


def mu_labor(n_s, l_tilde, chi, theta):
    '''
    Computes marginal utility with respect to labor supply
    '''

    b_ellipse, upsilon = ellip.estimation(theta, l_tilde)
    # b_ellipse is the constant defined in Equation (4.9) - has nothing to do
    # with savings
    mu_n = chi * (b_ellipse / l_tilde) * (n_s / l_tilde) ** (upsilon - 1)
    * (1 - (n_s / l_tilde) ** upsilon) ** ((1 - upsilon) / upsilon)

    return mu_n
