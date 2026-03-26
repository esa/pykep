import numpy as _np

from copy import deepcopy as _deepcopy
from math import log, exp, cos, acos, sin, pi, sqrt, atan2, asin
from scipy.special import logsumexp

def alpha2direct_py(alphas, T):
    """alpha2direct_py(x)

    Args:
        *alphas* (``array-like``): a sequence of transfer times encoded using the alpha encoding.

        *T* (:class:`float`): the total transfer time.

    Returns:
        :class:`list`:: The encoded transfer times
    """
    log_alphas = [log(item) for item in alphas]
    retval = [it / sum(log_alphas) * T for it in log_alphas]
    return retval


def direct2alpha_py(x):
    """direct2alpha(x)

    Args:
        *x* (``array-like``): a sequence of transfer times.

    Returns:
        :class:`list`:, :class:`float`: The alpha-encoded transfer times, the total transfer time (for cenvenience)
    """
    T = sum(x)
    retval = [exp(it / (-T)) for it in x]
    return retval, T


def eta2direct_py(x, max_tof):
    """eta2direct(x)

    Args:
        *x* (``array-like``): a sequence of transfer times encoded using the eta encoding.

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The encoded transfer times
    """
    n = len(x)
    # we assemble the times of flight
    T = [0] * n
    T[0] = max_tof * x[0]
    for i in range(1, len(T)):
        T[i] = (max_tof - sum(T[:i])) * x[i]
    return T


def direct2eta_py(x, max_tof):
    """direct2eta(x)

    Args:
        *x* (``array-like``):  a sequence of transfer times

        *max_tof* (:class:`float`): maximum time of flight (might actually be less according to the value of the last eta.)

    Returns:
        :class:`list`: The eta-encoded transfer times
    """
    retval = _deepcopy(x)
    retval[0] = x[0] / max_tof
    for i in range(1, len(x)):
        retval[i] = x[i] / (max_tof - sum(x[:i]))
    return retval

def uvV2cartesian(uvV):
    """This function converts the uvV encoding of a vector to cartesian coordinates

    Args:
        *uvV* (``array-like``): a sequence of 3 floats representing the vector in uvV encoding.

    Returns:
        :class:`list`:: The vector in cartesian coordinates
    """
    u, v, V = uvV
    theta = 2 * pi * u
    tmp = 2 * v - 1
    # Protecting against nans
    if (tmp*tmp) > 1:
        if tmp < -1:
            phi = pi/2
        elif tmp > 1:
            phi = -pi/2
    else:
        phi = acos(2 * v - 1) - pi / 2

    Vinfx = V * cos(phi) * cos(theta)
    Vinfy = V * cos(phi) * sin(theta)
    Vinfz = V * sin(phi)
    return [Vinfx, Vinfy, Vinfz]


def cartesian2uvV(V):
    """This function converts the cartesian coordinates of a vector to uvV encoding

    Args:
        *V* (``array-like``): a sequence of 3 floats representing the vector in cartesian coordinates.

    Returns:
        :class:`list`:: The vector in uvV encoding
    """
    v_norm = sqrt(V[0] ** 2 + V[1] ** 2 + V[2] ** 2)
    theta = atan2(V[1], V[0])
    sin_phi = V[2] / v_norm
    return [theta / 2 / pi, (-sin_phi + 1) / 2, v_norm]

def compute_softmax_and_jacobian(logits):
    """This function computes softmax and its Jacobian
    using a numerically stable formulation. The softmax transforms
    an unconstrained real-valued vector into a probability simplex (all
    weights positive and summing to one).

    Args:
        *logits* (``array-like``): a sequence of floats representing the unnormalized
            log-weights (logits) to be transformed via softmax.

    Returns:
        ``tuple``: A pair ``(weights, J)`` where:

        - *weights* (``numpy.ndarray``): the softmax-normalized weights,
          each in (0, 1) and summing to 1.
        - *J* (``numpy.ndarray``): the ``(n, n)`` Jacobian matrix of the
          softmax map, where ``J[i, j] = weights[i] * (delta_ij - weights[j])``.
    """
    w = _np.array(logits, dtype=float)

    # Subtract logsumexp for numerical stability before exponentiating
    logw = w - logsumexp(w)
    weights = _np.exp(logw)

    # Build the Jacobian of the softmax: J_ij = w_i * (δ_ij - w_j)
    # This follows from differentiating softmax_i w.r.t. input w_j
    n = len(w)
    J = _np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            J[i, j] = weights[i] * ((1.0 if i == j else 0.0) - weights[j])

    return weights, J