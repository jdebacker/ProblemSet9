import numpy as np

# household functions


def FOCs(b_sp1, *args):
    '''
    Use the HH FOCs to solve for optimal saves, b2 and b3
    (1) u'(c1) = beta * (1 + r)* u'(c2) -> b2
    (2) u'(c2) = beta * (1 + r)* u'(c3) -> b3
    '''
    beta, sigma, r, w, n, b_init = args
    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)
    c = get_c(r, w, n, b_s, b_sp1)
    mu_c = u_prime(c, sigma)
    lhs_euler = mu_c
    rhs_euler = beta * (1+r) * mu_c
    foc_errors = lhs_euler[:-1] - rhs_euler[1:]
    
    return foc_errors



def get_c(r, w, n, b_s, b_sp1):
    '''
    Use the budget constraint to solve for consumption
    '''
    c = w * n + (1 + r) * b_s - b_sp1

    return c


def u_prime(c, sigma):
    '''
    Marginal utility of consumption
    '''
    mu_c = c ** -sigma

    return mu_c

def mu_labor(n, chi, *args, b_sp1):
    '''
    Use FOC's for disutilitity of labor to get labor supply
    
    Refer to eqn 4.9 of ch 4 for the calculation of the 
    marginal utility of labor calculation
    
    n = household labor supply
    chi = utility weight on labor supply
    
    This model solves for the wage rate (w) given the exogeneous nature
    of labor supply
    
    There will be S labor supply Euler equations below
    
    ''' 
    
    theta = 0.9 # frisch elasticity of labor
    l_tilde = 1.0 # Maximum labor supply
    
    b, nu = elliptical_u_est.estimator(theta, l_tilde)
    beta, sigma, r, w, n, b_init = args
    b_s = np.append(b_init, b_sp1)
    b_sp1 = np.append(b_sp1, 0.0)
    c = get_c(r, w, n, b_s, b_sp1)
    mu_c = u_prime(c, sigma)
    
    lhs = w * mu_c
    rhs = chi * (b / l_tilde) * (n / l_tilde) ** (nu - 1) * (1 - (n / l_tilde) ** nu) ** ((1 - nu) / nu )
    foc_errors = lhs - rhs
    
    return foc_errors
