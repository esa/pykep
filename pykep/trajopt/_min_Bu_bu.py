import numpy as np

# from numba import njit
import scipy
import pygmo as pg


# @njit(cache=True)
def _stereo2cartesian(x, eij, ek):
    """Transforms the stereographic projected point back
    into the unit sphere.

    Args:
        x (array (2,)): The stereographic coordinates
        eij (array (2,3)): The basis unit vectors completing ek
        ek (array (3,)): The origin for the stereographic projection (also a unit vector)

    Returns:
        _type_: _description_
    """
    # Stereographic parametrization
    denom = 1 + x[0] * x[0] + x[1] * x[1]
    X = 2 * x[0] / denom
    Y = 2 * x[1] / denom
    Z = (1 - x[0] * x[0] - x[1] * x[1]) / denom
    return X * eij[0] + Y * eij[1] + Z * ek


# @njit(cache=True)
def _fitness_and_gradient(x, B, b, eij, ek):
    """Computes the fitness |Bu|- b u and its gradient with respect to
    the stereographic parametrization

    Args:
        x (array (2,)): The stereographic coordinates
        B (array (3,3)): The matrix in Bu
        b (array(3,1)): the vector in |BU| - b u
        eij (array (2,3)): The basis unit vectors completing ek
        ek (array (3,)): The origin for the stereographic projection (also a unit vector)

    Returns:
        float, array (2,): fitness and gradient
    """

    # Stereographic parametrization
    u = _stereo2cartesian(x, eij, ek)

    # Derivative of stereographic parametrization
    denom2 = (1 + x[0] * x[0] + x[1] * x[1]) ** 2
    dx1 = 2 * (1 - x[0] * x[0] + x[1] * x[1]) / denom2
    dy1 = -4 * x[0] * x[1] / denom2
    dz1 = -4 * x[0] / denom2
    du1 = dx1 * eij[0] + dy1 * eij[1] + dz1 * ek
    dx2 = -4 * x[0] * x[1] / denom2
    dy2 = 2 * (1 + x[0] * x[0] - x[1] * x[1]) / denom2
    dz2 = -4 * x[1] / denom2
    du2 = dx2 * eij[0] + dy2 * eij[1] + dz2 * ek

    du = np.vstack((du1, du2)).transpose()

    # Objective function and derivative
    Bu_norm = np.linalg.norm(B @ u) + 1e-18
    obj = Bu_norm - b @ u
    dobj = (B @ u) @ (B @ du) / Bu_norm - b @ du
    return obj, dobj


class _primer_vector_surrogate_udp:
    """min (|Bu| - bu)"""

    def __init__(self, B, b):
        self.B = B
        self.b = b

    def fitness(self, x):
        x = (np.array(x) / np.linalg.norm(x)).reshape(3, 1)
        f = np.linalg.norm(self.B @ x) - self.b @ x
        return [f[0]]

    def gradient(self, x):
        u = np.array(x).reshape(3, 1)
        norm_u = np.linalg.norm(u)
        uhat = u / norm_u
        grad_uhat_u = np.eye(3) / norm_u - u @ u.T / norm_u**3
        grad_f_uhat = 1.0 / (np.linalg.norm(self.B @ uhat)) * uhat.T @ self.B.T @ self.B
        grad_f_uhat -= self.b.reshape(1, 3)
        return (grad_f_uhat @ grad_uhat_u)[0]

    def get_bounds(self):
        lb = [-1, -1, -1]
        ub = [1, 1, 1]
        return (lb, ub)


# min_u (|Bu| - b u), where u is a unit vector
def minBu_bu(B, b):
    b_norm = np.linalg.norm(b)
    # Easy way out
    if b_norm < 1:
        return b_norm, np.array([0, 0, 0])

    # Compute the singular value decomposition (used to bound |Bu| as well as to provide one IG)
    svd = np.linalg.svd(B)
    # Second easy way out: smallest singular value
    sv = svd[1][-1]
    if b_norm < 1 + sv:
        return b_norm - sv, np.array([0, 0, 0])
    # Corresponding singular value vector
    u_svd = svd[2][-1, :]

    # We must construct an orthonormal basis defining the stereographic projection.
    # (we know at this point that norm b is not vanishing if we activate bound=True)
    ek = b / b_norm
    eij = scipy.linalg.null_space(ek.reshape((1, 3))).transpose()

    # The singular value vector direction is determined ...
    # ... by forcing it in the plane where b u is positive
    # This way the IG is a good one.
    if b.T @ u_svd < 0:
        u_svd = -u_svd

    # First one is with initial guess = bhat -> projected is always [0,0]
    optim1 = scipy.optimize.minimize(
        lambda x: _fitness_and_gradient(x, B, b, eij, ek),
        [0.0, 0.0],
        method="L-BFGS-B",
        jac=True,
    )

    if sv > 1e-8:
        # Second one is with svd initial guess
        xy = eij @ u_svd
        z = b @ u_svd / b_norm
        p_svd = xy / (1 + z)
        optim2 = scipy.optimize.minimize(
            lambda x: _fitness_and_gradient(x, B, b, eij, ek),
            p_svd,
            method="L-BFGS-B",
            jac=True,
        )

    else:
        optim2 = optim1
    if optim1.fun < optim2.fun:
        return - optim1.fun, _stereo2cartesian(optim1.x, eij, ek)
    else:
        return - optim2.fun, _stereo2cartesian(optim2.x, eij, ek)


# min_u (|Bu| - b u), where u is a unit vector
def minBu_bu_p(B, b):
    b_norm = np.linalg.norm(b)
    # Easy way out
    if b_norm < 1:
        return -b_norm, np.array([0, 0, 0])

    # Compute the singular value decomposition (used to bound |Bu| as well as to provide one IG)
    svd = np.linalg.svd(B)
    # Second easy way out: smallest singular value
    sv = svd[1][-1]
    if b_norm < 1 + sv:
        return -b_norm + sv, np.array([0, 0, 0])
    # Corresponding singular value vector
    u_svd = svd[2][-1, :]

    # The singular value vector direction is determined ...
    # ... by forcing it in the plane where b u is positive
    # This way the IG is a good one.
    if b.T @ u_svd < 0:
        u_svd = -u_svd

    # First one is with initial guess = bhat -> projected is always [0,0]
    prob = pg.problem(_primer_vector_surrogate_udp(B, b))
    algo = pg.algorithm(pg.nlopt(solver="slsqp"))
    pop1 = pg.population(prob)
    pop1.push_back(u_svd / np.linalg.norm(u_svd))
    pop1 = algo.evolve(pop1)

    pop2 = pg.population(prob)
    pop2.push_back(b / b_norm)
    pop2 = algo.evolve(pop2)

    if pop1.champion_f[0] < pop2.champion_f[0]:
        return - pop1.champion_f[0], pop1.champion_x / np.linalg.norm(pop1.champion_x)
    else:
        return - pop2.champion_f[0], pop2.champion_x / np.linalg.norm(pop2.champion_x)
