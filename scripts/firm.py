# Since we're not considering g_y, the expressions for r_t and w_t remain unchanged from what we did in class.

# firm functions

def get_r(L, K, alpha, A, delta):
    '''
    The interest rate implied by the firm FOC for the choice of capital
    Referring to equation (4.13)
    '''
    r = alpha * A * (L / K) ** (1 - alpha) - delta

    return r


def get_w(L, K, alpha, A, delta):
    '''
<<<<<<< HEAD
    The wage rate implied by the firm FOC for the choice of labor
    Referring to Equation (4.14)
=======
    The wage rate implied by the interest rate
>>>>>>> upstream/master
    '''
    w = ((1 - alpha) * A * ((r + delta) / (alpha * A)) **
         (alpha / (alpha - 1)))

    return w
