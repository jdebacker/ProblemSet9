import numpy as np
import elliptical_u_est as ellip

# household functions


def FOCs(b_sp1, n_s, *args):
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
    l_tilde: maximum amount of labor supply
    chi: scale parameter
    theta: Frisch elasticity of labor supply

    n is not contained in the remaining arguments anymore. If someone decides
    to change this, please provide detailed documentation
    on why you are doing so.


    Returns:
    foc_errors: A list where the first S-1 values are b2, b3, ..., bS, and the
    next S values are n1, n2, ..., nS.
    '''

    beta, sigma, r, w, b_init, l_tilde, chi, theta = args
    # When working on SS.py, note that b_sp1_guess is now of length S-1
    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)
    c = get_c(r, w, n_s, b_s, b_sp1)
    mu_c = mu_cons(c, sigma)
    mu_n = mu_labor(n_s, l_tilde, chi, theta)

    # First get the Euler equations defined by Equation (4.10) - S-1 of these
    lhs_euler_b = mu_c
    rhs_euler_b = beta * (1+r) * mu_c
    foc_errors_b = lhs_euler_b[:-1] - rhs_euler_b[1:]
    # The above line doesn't need to be changed since it works for a general S
    # as well (writing out the vectors helps to see this)

    # Next get the Euler equations defined by (4.9) - S of these
    lhs_euler_n = w * mu_c
    rhs_euler_n = mu_n
    foc_errors_n = lhs_euler_n - rhs_euler_n

    foc_errors = np.append(foc_errors_b, foc_errors_n)

    return foc_errors


def get_c(r, w, n, b_s, b_sp1):
    '''
    Use the budget constraint to solve for consumption
    '''
    c = w * n + (1 + r) * b_s - b_sp1

    return c


def mu_cons(c, sigma):
    '''
    Computes marginal utility with respect to consumption.

    Please note that this was initially called u_prime. If anyone finds
    function calls on u_prime, please change it to mu_cons
    Marginal utility of consumption
    '''
    mu_c = c ** -sigma

    return mu_c


def mu_labor(n_s, l_tilde, chi, theta):
    '''
    Computes marginal utility with respect to labor supply
    '''
    # b_ellipse is the constant defined in Equation (4.9) - has nothing to do
    # with savings
    b_ellipse, nu = ellip.estimation(theta, l_tilde)
    mu_n = chi * (b_ellipse / l_tilde) * (n_s / l_tilde) ** (nu - 1)
    * (1 - (n_s / l_tilde) ** nu) ** ((1 - nu) / nu)

    return mu_n
