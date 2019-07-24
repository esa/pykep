import numpy as np
from numba import jit


#@jit
def gravity_spherical_harmonic(x, r_planet, mu, c, s, n_max, m_max):
    """
    Calculate the gravitational acceleration due to the spherical harmonics gravity model supplied.

    Args:
        - x (``array-like``): (N x 3) array of Cartesian coordinates of satellite in frame defined by gravity model.
        - r_planet (``float``): Equatorial radius of central body.
        - mu (``float``): Gravitational parameter of central body.
        - c (``array-like``): Two-dimensional normalised C coefficient array: C[degree, order] = C_(n, m).
        - s (``array-like``): Two-dimensional normalised S coefficient array: S[degree, order] = S_(n, m).
        - n_max (``int``): Degree up to which to calculate the gravitational acceleration.
        - m_max (``int``): Order up to which to calculate the gravitational acceleration. Cannot be higher than n_max.
    
    Returns:
        - acc (``array-like``): (N x 3) array of the gravitational acceleration.

    Example::

        r, mu, c, s, n, m = pykep.util.load_gravity_model('gravity_models/Moon/glgm3150.txt')
        x = numpy.array([[459.5, 795.8773461, 1591.754692], [-459.5, -795.8773461, -1591.754692]])
        
        acc_zonal = pykep.util.gravity_spherical_harmonic(x, r, mu, c, s, 20, 0)
        acc_sqr = pykep.util.gravity_spherical_harmonic(x, r, mu, c, s, 20, 20)
        acc_high = pykep.util.gravity_spherical_harmonic(x, r, mu, c, s, 900, 900)

    .. note::

        This model was taken from a report by NASA:
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160011252.pdf
        This is the normalised gottlieb algorithm, as coded in MATLAB in the report and transferred to Python.
    """
    if not (len(x[0]) == 3):
        raise ValueError(f"Position must be an (N x 3) array. Shape of position is ({len(x)} x {len(x[0])}).")

    r = np.min(np.linalg.norm(x, axis=1))
    if r < r_planet:
        raise ValueError(f"Radial position is less than defined radius of central body ({np.minr}<{r_planet}).")

    acc = _gottlieb(x, r_planet, mu, c, s, n_max, m_max)

    return acc


def _gottlieb(x, r_planet, mu, c, s, n_max, m_max):
    norm1, norm2, norm11, normn10, norm1m, norm2m, normn1 = _calculate_normalisation_parameters(n_max)

    n_coor = x.shape[0]

    r = np.sqrt(x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2)

    r_inverted = 1 / r

    x_r = x[:, 0] * r_inverted
    y_r = x[:, 1] * r_inverted
    z_r = x[:, 2] * r_inverted

    rp_r = r_planet * r_inverted
    rp_rn = np.copy(rp_r)

    mu_r2 = mu * r_inverted * r_inverted

    shape = (max(2, n_max + 1), max(2, n_max + 2), n_coor)
    p = np.zeros(shape)

    p[0, 0] = 1
    p[1, 0, :] = np.sqrt(3) * z_r
    p[1, 1] = np.sqrt(3)

    # sectorial ALFs
    for n in range(2, n_max + 1):
        p[n, n] = norm11[n] * p[n - 1, n - 1] * (2 * n - 1)

    shape = (max(3, n_max + 1), n_coor)
    ctil = np.ones(shape)
    stil = np.zeros(shape)

    ctil[1] = x_r
    stil[1] = y_r

    sumh = np.zeros(n_coor)
    sumgm = np.ones(n_coor)
    sumj = np.zeros(n_coor)
    sumk = np.zeros(n_coor)

    for n in range(2, n_max + 1):
        rp_rn *= rp_r

        n2m1 = 2 * n - 1
        nm1 = n - 1
        nm2 = n - 2
        np1 = n + 1

        # tesseral ALFs
        p[n, nm1] = normn1[n, nm1] * z_r * p[n, n]

        # zonal ALFs
        p[n, 0] = (n2m1 * z_r * norm1[n] * p[nm1, 0] - nm1 * norm2[n] * p[nm2, 0]) / n
        p[n, 1] = (n2m1 * z_r * norm1m[n, 1] * p[nm1, 1] - n * norm2m[n, 1] * p[nm2, 1]) / nm1

        sumhn = normn10[n] * p[n, 1] * c[n, 0]
        sumgmn = p[n, 0] * c[n, 0] * np1

        if m_max > 0:
            for m in range(2, n - 1):
                p[n, m] = (n2m1 * z_r * norm1m[n, m] * p[nm1, m] -
                           (nm1 + m) * norm2m[n, m] * p[nm2, m]) / (n - m)

            sumjn = np.zeros(n_coor)
            sumkn = np.zeros(n_coor)

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

    _lambda = sumgm + z_r * sumh

    acc = np.zeros(x.shape)

    acc[:, 0] = - mu_r2 * (_lambda * x_r - sumj)
    acc[:, 1] = - mu_r2 * (_lambda * y_r - sumk)
    acc[:, 2] = - mu_r2 * (_lambda * z_r - sumh)

    return acc


@jit
def _calculate_normalisation_parameters(n_max):
    """
    Calculate the normalisation parameters for the normalised Gottlieb algorithm. The mapping from the parameters defined in the report
    (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20160011252.pdf) to the names used in the code looks like

        - lambda_n-1(n)       = norm1
        - lambda_n-2(n)       = norm2
        - lambda_n-1_n-1(n)   = norm11
        - lambda_n-1_m(n, m)  = norm1m
        - lambda_n-2_m(n, m)  = norm2m
        - lambda_n_m+1(n, m)  = normn1
        - lamdba_n_m+1(n, 0)  = normn10

    Args:
        - n_max (``int``): Maximum degree.

    Returns:
        Arrays of normalisation parameters.
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
