import pykep as pk
import numpy as _np
from copy import deepcopy as _deepcopy


def add_sf_leg(
    ax,
    sf: pk.leg.sims_flanagan,
    units=pk.AU,
    N=30,
    show_midpoints=False,
    show_gridpoints=False,
    show_throttles=False,
    use_alpha=False,
    length=0.1,
    arrow_length_ratio=0.05,
    **kwargs
):
    """
    Add a trajectory leg of Sims-Flanagan problem to a 3D matplotlib Axes.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): The 3D Axes object to which the trajectory leg will be added.

        *sf* (:class:`~pykep.leg.sims_flanagan`): The Sims-Flanagan object containing relevant information.

        *units* (:class:`float`, optional): The unit conversion factor for plotting. Default is pk.AU.

        *N* (:class:`int`, optional): The number of points to generate along each segment of the trajectory. Default is 10.

        *show_midpoints* (:class:`bool`, optional): If True, midpoints of each segment are shown. Default is False.

        *show_gridpoints* (:class:`bool`, optional): If True, gridpoints of the trajectory are shown. Default is False.

        *show_throttles* (:class:`bool`, optional): If True, thrust vectors at midpoints are shown. Default is False.

        *use_alpha* (:class:`bool`, optional): If True, Alpha encoding was used for leg segments (changes propagation times to sf.talphas). Default is False.

        *length* (:class:`float`, optional): The length of the thrust vectors when show_throttles is True. Default is 0.1.

        *arrow_length_ratio* (:class:`float`, optional): The ratio of arrow length to the total length when show_throttles is True. Default is 0.05.

        *\\*\\*kwargs*: Additional keyword arguments to pass to the Axes3D.plot function.

    Notes:
        - This function visualizes a Sims-Flanagan trajectory leg on the provided 3D Axes object.
        - Midpoints, gridpoints, and thrust vectors can be optionally shown based on the provided parameters.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the Sims-Flanagan leg added.
    """
    # We extract the number of segments from the leg.
    nseg = sf.nseg
    nseg_fwd = sf.nseg_fwd
    nseg_bck = sf.nseg_bck

    if not use_alpha:
        dt = sf.tof / nseg
        c = sf.max_thrust * dt

    # We start the forward pass of the Sims-Flanagan model------------------------------------------------------------------------
    pos_fwd = []
    pos_m_fwd = []
    throttles_fwd = []
    rv = _deepcopy(sf.rvs)
    mass_fwd = sf.ms
    # Append to plotting  data
    pos_fwd.append(rv[0])

    for i in range(nseg_fwd):
        if use_alpha:
            dt = sf.talphas[i]
            c = sf.max_thrust * dt

        # compute the dv
        throttles = sf.throttles[3 * i : 3 * i + 3]
        throttles_fwd.append(throttles)
        dv = _np.linalg.norm(throttles) * c / mass_fwd
        # plot it in a color that is proportional to the strength
        color = (
            0.25 + (0.80 - 0.25) * min(1.0, _np.linalg.norm(throttles)),
            0.41 + (0.36 - 0.41) * min(1.0, _np.linalg.norm(throttles)),
            0.88 + (0.36 - 0.88) * min(1.0, _np.linalg.norm(throttles)),
        )
        pk.plot.add_ballistic_arc(
            ax, rv, dt / 2, sf.mu, units=units, N=N, c=color, **kwargs
        )
        # propagate for dt/2
        rv = list(pk.propagate_lagrangian(rv, tof=dt / 2, mu=sf.mu, stm=False))
        # register the position as a mid-point state
        pos_m_fwd.append(rv[0])
        # add dv to the state (now dimensional)
        rv[1] = [a + b * c / mass_fwd for a, b in zip(rv[1], throttles)]
        # update the mass (increases)
        mass_fwd *= _np.exp(-dv / sf.veff)
        # 2 - propagate for dt/2
        pk.plot.add_ballistic_arc(
            ax, rv, dt / 2, sf.mu, units=units, N=N, c=color, **kwargs
        )
        rv = pk.propagate_lagrangian(rv, tof=dt / 2, mu=sf.mu, stm=False)
        pos_fwd.append(rv[0])
    pos_fwd = _np.array(pos_fwd)
    pos_m_fwd = _np.array(pos_m_fwd)
    throttles_fwd = _np.array(throttles_fwd)

    # We plot optionally gridpoints, the low-thrust or the mid points
    if show_gridpoints:
        ax.plot(
            pos_fwd[:, 0] / units, pos_fwd[:, 1] / units, pos_fwd[:, 2] / units, "k."
        )
    if show_midpoints:
        ax.plot(
            pos_m_fwd[:, 0] / units,
            pos_m_fwd[:, 1] / units,
            pos_m_fwd[:, 2] / units,
            "kx",
        )
    if show_throttles:
        ax.quiver(
            pos_m_fwd[:, 0] / units,
            pos_m_fwd[:, 1] / units,
            pos_m_fwd[:, 2] / units,
            throttles_fwd[:, 0],
            throttles_fwd[:, 1],
            throttles_fwd[:, 2],
            length=length,
            color="indianred",
            arrow_length_ratio=arrow_length_ratio,
        )

    # We start the backward pass of the Sims-Flanagan model------------------------------------------------------------------------
    pos_bck = []
    pos_bck_m = []
    throttles_bck = []
    rv = _deepcopy(sf.rvf)
    mass_bck = sf.mf
    # Append to plotting  data
    pos_bck.append(rv[0])

    for i in range(nseg_bck):
        if use_alpha:
            dt = sf.talphas[nseg - i - 1]
            c = sf.max_thrust * dt
        # compute the dv (first non dimensional)
        throttles = sf.throttles[nseg * 3 - 3 - 3 * i : nseg * 3 - 3 * i]
        throttles_bck.append(throttles)
        dv = _np.linalg.norm(throttles) * c / mass_bck
        # plot it in a color that is proportional to the strength
        color = (
            0.25 + (0.80 - 0.25) * min(1.0, _np.linalg.norm(throttles)),
            0.41 + (0.36 - 0.41) * min(1.0, _np.linalg.norm(throttles)),
            0.88 + (0.36 - 0.88) * min(1.0, _np.linalg.norm(throttles)),
        )
        pk.plot.add_ballistic_arc(
            ax, rv, -dt / 2, sf.mu, units=units, N=N, c=color, **kwargs
        )
        # propagate for dt/2
        rv = list(pk.propagate_lagrangian(rv, tof=-dt / 2, mu=sf.mu, stm=False))
        # register the position as a mid-point state
        pos_bck_m.append(rv[0])
        # add it to the state (now dimensional)
        rv[1] = [a - b * c / mass_bck for a, b in zip(rv[1], throttles)]
        # update the mass
        mass_bck *= _np.exp(dv / sf.veff)
        # 2 - propagate for dt/2
        pk.plot.add_ballistic_arc(
            ax, rv, -dt / 2, sf.mu, units=units, N=N, c=color, **kwargs
        )
        rv = pk.propagate_lagrangian(rv, tof=-dt / 2, mu=sf.mu, stm=False)
        pos_bck.append(rv[0])
    pos_bck = _np.array(pos_bck)
    pos_bck_m = _np.array(pos_bck_m)
    throttles_bck = _np.array(throttles_bck)

    # We plot optionally gridpoints, the low-thrust or the mid points
    if show_gridpoints:
        ax.plot(
            pos_bck[:, 0] / units, pos_bck[:, 1] / units, pos_bck[:, 2] / units, "k."
        )
    if show_midpoints:
        ax.plot(
            pos_bck_m[:, 0] / units,
            pos_bck_m[:, 1] / units,
            pos_bck_m[:, 2] / units,
            "kx",
        )
    if show_throttles:
        ax.quiver(
            pos_bck_m[:, 0] / units,
            pos_bck_m[:, 1] / units,
            pos_bck_m[:, 2] / units,
            throttles_bck[:, 0],
            throttles_bck[:, 1],
            throttles_bck[:, 2],
            length=length,
            color="indianred",
            arrow_length_ratio=arrow_length_ratio,
        )

    return ax
