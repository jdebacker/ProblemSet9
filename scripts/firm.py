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
    The wage rate implied by the firm FOC for the choice of labor
    Referring to Equation (4.14)
    '''
    w = (1 - alpha) * A * (K/L)**alpha

    return w
