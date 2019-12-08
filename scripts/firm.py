# firm functions

import math

def get_r(L, K, alpha, A, g_y, delta):
    '''
    The interest rate implied by the firm FOC for the choice of capital
    '''
    r = alpha * A * (math.exp(g_y * t))**(1-alpha) * (L/K)**(1-alpha) - delta

    return r


def get_w(r, alpha, A, delta):
    '''
    The wage rate implied by the interest rate
    '''
    w = ((1 - alpha) * A * (K/L)**(alpha) * math.exp(g_y * t)

    return w
