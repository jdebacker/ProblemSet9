# firm functions
def get_Y (L, K, alpha, A, g_y, t):
    '''
    Updated based on eqn 8.14 
    '''
    import math
    
    Y = A * (K ** alpha) * (math.exp(g_y * t) * L) ** (1 - alpha)

    return Y
    
def get_r(Y, K, alpha, delta):
    '''
    Updated based on eqn 8.16 
    '''
    r = alpha * (Y / K) - delta

    return r


def get_w(Y, alpha, L):
    '''
    Updated based on eqn 8.17
    '''
    w = (1 - alpha) * (Y / L)

    return w
