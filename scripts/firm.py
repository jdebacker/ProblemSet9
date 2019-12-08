# firm functions

def get_r(L, K, G, alpha, A, delta):
    '''
    The interest rate implied by the firm FOC for the choice of capital
    G is the efficiency of labor with constant growth rate of g_y
    I am going to define G in TPI as:
          G = exp(g_y * t)
          Or
          G[t+1] = G[t] + G[t] * g_y with a for loop for t in range of S
    '''
    r = alpha * A * G**(1-alpha) * (L/K)**(1-alpha) - delta

    return r


def get_w(L, K, G, alpha, A):
    '''
    The wage rate implied by the firm FOC for the choice of labor
    G is the efficiency of labor with constant growth rate of g_y
    NO NEED for delta here. Consider that in using this function.
    '''
    w = (1 - alpha) * A *  G**(1-alpha) * (K/L)**(alpha) 

    return w
