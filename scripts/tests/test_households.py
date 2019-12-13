import numpy as np
import households


def test_FOC_labor():
    '''
    Test of households.FOC_labor()
    '''
    n_s = 0.4
    c = 0.7575
    w = 1.2
    sigma = 2.0
    l_tilde = 1
    chi = 0.5
    b_ellipse = 0.47
    upsilon = 1.2
    n_args = (sigma, w, l_tilde, chi, b_ellipse, upsilon)
    expected_error = 1.881986044
    test_error = households.FOC_labor(c, n_s, n_args)

    assert np.allclose(expected_error, test_error)


def test_FOC_save():
    '''
    Test of households.FOC_save()
    '''
    c = np.array([0.7575, 0.68])
    r = 0.05
    rho_s = 0.0
    beta = 0.96
    sigma = 2.0
    params = (rho_s, beta, sigma, r)
    expected_error = -0.437182264

    test_error = households.FOC_save(c, params)

    assert np.allclose(expected_error, test_error)


def test_get_c():
    '''
    Test of households.get_c()
    '''
    r = 0.05
    w = 1.2
    n_s = 0.4
    b_s = 0.55
    b_sp1 = 0.4
    BQ = 0.1
    expected_c = 0.7575

    test_c = households.get_c(r, w, n_s, b_s, b_sp1, BQ)

    assert np.allclose(expected_c, test_c)


def test_mu_cons():
    '''
    Test of households.mu_cons()
    '''
    c = 0.5
    sigma = 1.0
    expected_mu = 2.0

    test_mu = households.mu_cons(c, sigma)

    assert np.allclose(expected_mu, test_mu)


def test_mu_labor():
    '''
    Test households.mu_labor()
    '''
    n_s = 0.4
    l_tilde = 1
    chi = 0.5
    b_ellipse = 0.47
    upsilon = 1.2
    expected_mu = 0.209312195

    test_mu = households.mu_labor(n_s, l_tilde, chi, b_ellipse, upsilon)

    assert np.allclose(expected_mu, test_mu)


def test_FOCs():
    '''
    Test of households.FOCs()
    '''
    b_sp1 = 0.4
    n_s = np.array([0.4, 0.4])
    b_init = 0.55
    BQ = 0.1
    rho_s = 0.0
    beta = 0.96
    sigma = 2.0
    l_tilde = 1.0
    chi = 0.5
    b_ellipse = 0.47
    upsilon = 1.2
    r = 0.05
    w = 1.2
    omega_SS = np.ones(n_s.shape[0])
    g_n = 0.0
    method = 'TPI'
    foc_args = (b_init, BQ, rho_s, omega_SS, g_n, beta, sigma, l_tilde, chi,
                b_ellipse, upsilon, r, w, method)
    expected_errors = np.array([0.734748532, 1.881986044, 0.990687805])

    test_error = households.FOCs(b_sp1, n_s, *foc_args)

    assert np.allclose(expected_errors, test_error)
