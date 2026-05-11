import pykep as _pk
import numpy as _np

def add_lambert(ax, lp, N: int = 60, sol: int = 0, units=_pk.AU, **kwargs):
    """Add Lambert's problem solution trajectory to a 3D matplotlib Axes.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): The 3D Axes object to which the trajectory will be added.

        *lp* (:class:`~pykep.lambert_problem`): The Lambert's problem object containing relevant information.

        *N* (:class:`int`, optional): The number of points to generate along the trajectory. Default is 60.

        *sol* (:class:`int`, optional): The solution index for Lambert's problem. Default is 0.

        *units* (:class:`float`, optional): The unit conversion factor for plotting. Default is _pk.AU.

        *\\*\\*kwargs*: Additional keyword arguments to pass to the Axes3D.plot function.

    Raises:
        ValueError: If the specified solution index (sol) is greater than twice the maximum number of revolutions (Nmax).

    Notes:
        This function visualizes the trajectory of Lambert's problem solution on the provided 3D Axes object.
        The trajectory is computed based on the initial and final positions, velocities, and gravitational parameter.
        The integration grid is defined based on the specified solution index and angle computations.
        The resulting trajectory is plotted on the given Axes object using the provided unit conversion factor.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the Lambert's problem trajectory added.
    """
    # We check that the requested arc exists
    if sol > lp.Nmax * 2:
        raise ValueError(
            "sol must be in 0 .. NMax*2 \n * Nmax is the maximum number of revolutions for which there exist a solution to the Lambert's problem \n * You can compute Nmax calling the Nmax() method of the lambert_problem object"
        )

    # We extract the relevant information from the Lambert's problem
    r0 = lp.r0
    v0 = lp.v0[sol]
    r1 = lp.r1
    mu = lp.mu

    # We define the integration grid
    if sol == 0:
        # We have the problem of computing, in the direction of motion, 
        # the true anomaly difference between the starting and the final point
        # The folloiwng works in many tested geometries
        h_unit = _np.cross(r0, v0)
        h_unit /= _np.linalg.norm(_np.cross(r0, v0))
        cross_r0_r1 = _np.cross(r0, r1)
        dot_r0_r1 = _np.dot(r0, r1)
        theta = _np.arctan2(_np.dot(h_unit, cross_r0_r1), dot_r0_r1) # in [-pi,pi]
        if theta < 0:
            theta += 2 * _np.pi # in [0 2pi]
        thetagrid = _np.linspace(0, theta, N)
    else:
        thetagrid = _np.linspace(0, 2 * _np.pi, N)
        
    # Compute the posvel at all points
    res = _pk.plot.propagate_lagrangian_theta_v(
        rv=[r0, v0], thetas=thetagrid, mu=mu, stm=False
    )
    pos = res[:, :3] / units

    # And plot
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], **kwargs)

    return ax