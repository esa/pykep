import pykep as _pk
from copy import deepcopy as _deepcopy

def add_ballistic_arc(ax, rv0, tof, mu, units=_pk.AU, N=60, **kwargs):
    """
    Add a ballistic trajectory arc to a 3D matplotlib Axes.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): The 3D Axes object to which the ballistic arc will be added.

        *rv0* (:class:`list`): Initial position and velocity vector in Cartesian coordinates.

        *tof* (:class:`float`): Time of flight for the ballistic arc.

        *mu* (:class:`float`): Gravitational parameter of the central body.

        *units* (:class:`float`): The unit conversion factor for plotting. Default is pk.AU.

        *N* (:class:`int`): The number of points to generate along the ballistic arc. Default is 60.

        *\\*\\*kwargs*: Additional keyword arguments to pass to the Axes3D.plot function.

    Notes:
        - This function visualizes a ballistic trajectory arc on the provided 3D Axes object.
        - The trajectory is computed based on the initial position, velocity, time of flight, and gravitational parameter.
        - The resulting trajectory is plotted on the given Axes object using the provided unit conversion factor.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the ballistic arc trajectory added.
    """
    # We define the integration time ...
    dt = tof / (N - 1)

    # ... and calculate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    rv = _deepcopy(rv0)
    for i in range(N):
        x[i] = rv[0][0] / units
        y[i] = rv[0][1] / units
        z[i] = rv[0][2] / units
        rv = _pk.propagate_lagrangian(rv, dt, mu, False)

    # And plot
    ax.plot(x, y, z, **kwargs)
    return ax