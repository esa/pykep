import numpy as np
from numba import jit


@jit
def gravity_spherical_harmonic(x, *args):
    """Calculate the gravitational acceleration due to the spherical harmonics gravity model supplied.

    :param x: (N x 3) array of Cartesian coordinates of satellite in frame defined by gravity model.
    :type x: numpy.ndarray.
    :return: (N x 3) array of the gravitational acceleration.
    :rtype: numpy.ndarray.
    """
    acc = np.zeros(x.shape)
    for i in range(len(x[0])):
        acc[i] = _gottlieb(x[i], *args)

    return acc


@jit
def _gottlieb(x, r_planet, mu, c, s, n_max, m_max):
    """Calculate the gravitational acceleration on one set of coordinates due to
    the spherical harmonics gravity model supplied.

    :param x: Cartesian coordinates of satellite in frame defined by gravity model.
    :type x: numpy.ndarray.
    :param r_planet: Equatorial radius of planet.
    :type r_planet: float.
    :param mu: Gravitational parameter of planet.
    :type mu: float.
    :param c: Two-dimensional normalised C coefficient array: C[degree, order] = Cnm
    :type c: numpy.ndarray
    :param s: Two-dimensional normalised S coefficient array: S[degree, order] = Snm
    :type s: numpy.ndarray.
    :param n_max: Degree up to which to calculate the gravitational acceleration.
    :type n_max: int.
    :param m_max: Order up to which to calculate the gravitational acceleration. Cannot be higher than n_max.
    :type m_max: int.
    :return: Vector of the gravitational acceleration in a cartesian coordinate frame defined by the gravity model.
    :rtype: numpy.ndarray.

    This model was taken from a report by NASA:
    https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160011252.pdf
    This is the normalised gottlieb algorithm, as coded in MATLAB in the report and transferred to Python.
    """
    norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 = _calculate_normalisation_parameters(n_max)

    r = np.linalg.norm(x)

    r_inverted = 1 / r

    x_r = x[0] * r_inverted
    y_r = x[1] * r_inverted
    z_r = x[2] * r_inverted

    ep = z_r

    rp_r = r_planet * r_inverted
    rp_rn = rp_r

    mu_r2 = mu * r_inverted * r_inverted

    p = np.zeros((max(2, n_max + 1), max(2, n_max + 2)))

    p[0, 0] = 1
    p[1, 0] = np.sqrt(3) * ep
    p[1, 1] = np.sqrt(3)

    # sectorial ALFs
    for n in range(2, n_max + 1):
        p[n, n] = norm11[n] * p[n - 1, n - 1] * (2 * n - 1)

    ctil = np.zeros(max(3, n_max + 1))
    stil = np.zeros(max(3, n_max + 1))

    ctil[0] = 1
    ctil[1] = x_r
    stil[1] = y_r

    sumh = 0
    sumgm = 1
    sumj = 0
    sumk = 0

    for n in range(2, n_max + 1):
        rp_rn *= rp_r

        n2m1 = 2 * n - 1
        nm1 = n - 1
        nm2 = n - 2
        np1 = n + 1

        # tesseral ALFs
        p[n, nm1] = normn1[n, nm1] * ep * p[n, n]

        # zonal ALFs
        p[n, 0] = (n2m1 * ep * norm1[n] * p[nm1, 0] - nm1 * norm2[n] * p[nm2, 0]) / n
        p[n, 1] = (n2m1 * ep * norm1m[n, 1] * p[nm1, 1] - n * norm2m[n, 1] * p[nm2, 1]) / nm1

        sumhn = normn10[n] * p[n, 1] * c[n, 0]
        sumgmn = p[n, 0] * c[n, 0] * np1

        if m_max > 0:
            for m in range(2, n - 1):
                p[n, m] = (n2m1 * ep * norm1m[n, m] * p[nm1, m] -
                           (nm1 + m) * norm2m[n, m] * p[nm2, m]) / (n - m)

            sumjn = 0
            sumkn = 0

            ctil[n] = ctil[1] * ctil[nm1] - stil[1] * stil[nm1]
            stil[n] = stil[1] * ctil[nm1] + ctil[1] * stil[nm1]

            for m in range(1, min(n, m_max) + 1):
                mm1 = m - 1
                mp1 = m + 1
                mxpnm = m * p[n, m]

                bnmtil = c[n, m] * ctil[m] + s[n, m] * stil[m]

                sumhn += normn1[n, m] * p[n, mp1] * bnmtil
                sumgmn += (n + m + 1) * p[n, m] * bnmtil

                bnmtm1 = c[n, m] * ctil[mm1] + s[n, m] * stil[mm1]
                anmtm1 = c[n, m] * stil[mm1] - s[n, m] * ctil[mm1]

                sumjn += mxpnm * bnmtm1
                sumkn -= mxpnm * anmtm1

            sumj += rp_rn * sumjn
            sumk += rp_rn * sumkn

        sumh += rp_rn * sumhn
        sumgm += rp_rn * sumgmn

    _lambda = sumgm + ep * sumh

    acceleration = np.zeros(3)

    acceleration[0] = - mu_r2 * (_lambda * x_r - sumj)
    acceleration[1] = - mu_r2 * (_lambda * y_r - sumk)
    acceleration[2] = - mu_r2 * (_lambda * z_r - sumh)

    return acceleration


@jit
def _calculate_normalisation_parameters(n_max):
    """Calculate the parameters as defined in the report by NASA.
    (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160011252.pdf)

        lambda_n-1(n)       = norm1
        lambda_n-2(n)       = norm2
        lambda_n-1_n-1(n)   = norm11
        lambda_n-1_m(n, m)  = norm1m
        lambda_n-2_m(n, m)  = norm2m
        lambda_n_m+1(n, m)  = normn1
        lamdba_n_m+1(n, 0)  = normn10

    :param n_max: Maximum degree.
    :type n_max: int.
    :return: arrays of normalisation parameters.
    """
    m_max = n_max
    n_max += 1

    norm1 = np.zeros(n_max)
    norm2 = np.zeros(n_max)
    norm11 = np.zeros(n_max)
    normn10 = np.zeros(n_max)

    norm1m = np.zeros((n_max, m_max))
    norm2m = np.zeros((n_max, m_max))
    normn1 = np.zeros((n_max, m_max + 1))

    for n in range(2, n_max):
        norm1[n] = np.sqrt((2 * n + 1) / (2 * n - 1))
        norm2[n] = np.sqrt((2 * n + 1) / (2 * n - 3))
        norm11[n] = np.sqrt((2 * n + 1) / (2 * n)) / (2 * n - 1)
        normn10[n] = np.sqrt((n + 1) * n / 2)

        for m in range(1, n):  
            norm1m[n, m] = np.sqrt((n - m) * (2 * n + 1) / ((n + m) * (2 * n - 1)))  
            norm2m[n, m] = np.sqrt((n - m) * (n - m - 1) * (2 * n + 1) / ((n + m) * (n + m - 1) * (2 * n - 3)))  
            normn1[n, m] = np.sqrt((n + m + 1) * (n - m))  

    return norm1, norm2, norm11, normn10, norm1m, norm2m, normn1