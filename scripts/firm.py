# firm functions


def get_r(K, L, params):
    alpha, A, delta = params
    r = alpha * A * ((L / K) ** (1 - alpha)) - delta

    return r


def get_w(K, L, params):
    alpha, A = params
    w = (1 - alpha) * A * ((K / L) ** alpha)

    return w