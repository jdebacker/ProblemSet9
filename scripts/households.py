'''
Can someone please take a look at this function and make sure that this
accounts for everything in the household's problem?
I believe this script does account for S-periods, endogenous labor
supply, and population dynamics.
From what I can tell, the only factor that needs to be added to account
for population dynamics is rho_s, which has been added.
'''

import numpy as np
import elliptical_u_est as ellip
import aggregates as agg

# household functions
# FOC for savings
def FOC_save(c, params):
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
    rho_s: risk that someone alive at age-s will die at the end of that period
        and not be alive for age s+1
    l_tilde: maximum amount of labor supply
    chi: scale parameter
    theta: Frisch elasticity of labor supply

    n is not contained in the remaining arguments anymore. If someone decides
    to change this, please provide detailed documentation
    on why you are doing so.


    Returns:
    foc_errors_b: A list where the first S-1 values are b2, b3, ..., bS
    '''

    rho_s, beta, sigma, r = params

    mu_c = mu_cons(c, sigma)
    lhs = mu_c
    rhs = beta * (1 + r) * (1 - rho_s) * mu_c
    foc_errors_b = lhs[:-1] - rhs[1:]

    return foc_errors_b


# FOC for labor
def FOC_labor(c, n_s, args):
    '''
    Args:
    n_s: The labor supply values for each period. The call to this function
            should provide initial guesses for this variable.
    BQ_path: path of bequests
    rho_s: risk that someone alive at age-s will die at the end of that period
        and not be alive for age s+1
    l_tilde: maximum amount of labor supply
    chi: scale parameter
    theta: Frisch elasticity of labor supply


    Returns:
    foc_errors_n: A list of n2, n3, ..., nS
    '''
    sigma, w, l_tilde, chi, b_ellipse, upsilon = args
    mu_c = mu_cons(c, sigma)
    mu_n = mu_labor(n_s, l_tilde, chi, b_ellipse, upsilon)

    # Next get the Euler equations defined by (4.9) - S of these
    lhs_euler_n = w * mu_c
    rhs_euler_n = mu_n
    foc_errors_n = lhs_euler_n - rhs_euler_n

    return foc_errors_n


# Solve FOC_save and FOC_labor together
def FOCs(b_sp1, n_s, *args):
    (b_init, BQ, rho_s, omega, g_n, beta, sigma, l_tilde, chi,
     b_ellipse, upsilon, r, w, method) = args
    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)
    if method == 'SS':
        # if method is SS, solve for BQ, else (case of TPI), use BQ
        # as guessed in the "outer loop"
        BQ_params = (omega, g_n, rho_s)
        BQ = agg.get_BQ(b_s, r, BQ_params, method)
    c = get_c(r, w, n_s, b_s, b_sp1, BQ)
    b_args = (rho_s, beta, sigma, r)
    b_errors = FOC_save(c, b_args)
    n_args = (sigma, w, l_tilde, chi, b_ellipse, upsilon)
    n_errors = FOC_labor(c, n_s, n_args)
    errors = np.append(b_errors, n_errors)

    return errors


def get_c(r, w, n_s, b_s, b_sp1, BQ):
    '''
    Use the budget constraint to solve for consumption
    '''
    c = w * n_s + (1 + r) * b_s + BQ - b_sp1

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


def mu_labor(n_s, l_tilde, chi, b_ellipse, upsilon):
    '''
    Computes marginal utility with respect to labor supply
    '''
    mu_n = (
        chi * (b_ellipse / l_tilde) * (n_s / l_tilde) ** (upsilon - 1)
        * (1 - (n_s / l_tilde) ** upsilon) ** ((1 - upsilon) / upsilon))

    return mu_n
