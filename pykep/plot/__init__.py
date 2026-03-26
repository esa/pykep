# We import symbols explicitly with an underscore to hide them from the
# imported symbols
import matplotlib.pyplot as _plt
from mpl_toolkits.mplot3d import axes3d as _axes3d
import numpy as _np

from ._planet import add_planet_orbit, add_planet, add_solar_system, add_planets
from ._lambert import add_lambert
from ._ballistic import add_ballistic_arc
from ._sf_leg import add_sf_leg
from ._mit import add_mit

def make_3Daxis(**kwargs):
    """Constructs and returns a 3D axis.  All kwargs are forwarded to the call to `figure()` in matplotlib.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis
    """
    ax = _plt.figure(**kwargs).add_subplot(projection="3d")
    return ax


def add_sun(ax, **kwargs):
    """Adds the Sun to *ax*. All kwargs are forwarded to the scatter method of *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.
    """
    kwargs.setdefault("c", "y")
    kwargs.setdefault("s", 30)

    ax.scatter(0, 0, 0, **kwargs)
    return ax


def propagate_lagrangian_theta_v(rv=[[1, 0, 0], [0, 1, 0]], thetas=[0.0, _np.pi / 2], mu=1, stm=False):
    """propagate_lagrangian_theta_v(rv = [[1,0,0], [0,1,0]], thetas = [0, pi/2], mu = 1, stm = False)

    Propagates (Keplerian) the state for an assigned difference in True anomalies. It does not compute the State Transition Matrix
    A similar API is offered as that of :func:`~pykep.propagate_lagrangian`, but not identical. The function is essentially offered only for plotting
    purposes as to avoid (or anyway alleviate) the non uniform grid effect when plotting at equal time intervals.

    Args:
          *rv* (2D array-like): Cartesian components of the initial position vector and velocity [[x0, y0, z0], [v0, vy0, vz0]]. Defaults to [[1,0,0], [0,1,0]].

          *thetas* (:class:`numpy.ndarray` (N,)): true anomaly dfferences to plot at. Defaults to [0, pi/2].

          *mu* (:class:`float`): gravitational parameter. Defaults to 1.

          *stm* (:class:`bool`): requests the computations of the State Transition Matrix

    Returns:
          :class:`numpy.ndarray` (N,6): r and v, that is the final position and velocity after the propagation.

    Examples:
        >>> import pykep as pk
        >>> import numpy as np
        >>> r0 = [1,0,0]
        >>> v0 = [0,1,0]
        >>> thetas = [0, np.pi/2]
        >>> mu = 1
        >>> [r1,v1] = pk.propagate_lagrangian_theta_v(rv=[r0,v0], thetas = thetas, mu = mu)
    """
    if stm:
        raise NotImplementedError(
            "State Transition Matrix for propagate_lagrangian_theta_v not implemented yet."
        )
    N = len(thetas)
    retval = _np.zeros((N, 6))
    # Preliminary values
    r0 = _np.array(rv[0])
    v0 = _np.array(rv[1])
    R0 = _np.linalg.norm(r0)
    V02 = _np.dot(v0, v0)
    #energy = V02 / 2 - mu / R0
    # energy will be negative for hyperbolae
    #a = -mu / 2.0 / energy
    sigma0 = _np.dot(r0, v0) / _np.sqrt(mu)
    h = _np.linalg.norm(_np.cross(r0, v0))
    p = h * h / mu
    for i, theta in enumerate(thetas):
        Rf = (
            p
            * R0
            / (R0 + (p - R0) * _np.cos(theta) - _np.sqrt(p) * sigma0 * _np.sin(theta))
        )
        F = 1 - Rf / p * (1 - _np.cos(theta))
        G = R0 * Rf / h * _np.sin(theta)
        Ft = (
            _np.sqrt(mu)
            / R0
            / p
            * (sigma0 * (1 - _np.cos(theta)) - _np.sqrt(p) * _np.sin(theta))
        )
        Gt = 1 - R0 / p * (1 - _np.cos(theta))
        retval[i][:3] = F * r0 + G * v0
        retval[i][3:] = Ft * r0 + Gt * v0
    return retval
