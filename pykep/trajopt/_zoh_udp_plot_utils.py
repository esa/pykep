import numpy as _np


def set_3d_axes_equal(ax):
    """Set equal limits and aspect ratio on a three-dimensional axes."""
    limits = _np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    origin = _np.mean(limits, axis=1)
    radius = 0.5 * _np.max(_np.abs(limits[:, 1] - limits[:, 0]))
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])
    ax.set_box_aspect([1, 1, 1])


def draw_solar_sail(ax, state, alpha, beta, sail_size=0.05, state_to_cartesian=None):
    """Draw a solar-sail patch and normal at a Cartesian state.

    Args:
        ax (:class:`matplotlib.axes.Axes`): Three-dimensional axes.
        state (:class:`array-like`): Cartesian position and velocity, optionally
            followed by additional state components.
        alpha (:class:`float`): Sail cone angle from the radial direction.
        beta (:class:`float`): Sail clock angle in the transverse-normal plane.
        sail_size (:class:`float`, optional): Half-side length of the patch.
            Default is 0.05 in the plotted position units.
        state_to_cartesian (:class:`callable`, optional): Convert the dynamics
            state to Cartesian coordinates before drawing.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if state_to_cartesian is not None:
        state = state_to_cartesian(state)
    cart = _np.asarray(state, dtype=float)
    radial = cart[:3] / _np.linalg.norm(cart[:3])
    normal_orbit = _np.cross(cart[:3], cart[3:6])
    normal_orbit /= _np.linalg.norm(normal_orbit)
    transverse = _np.cross(normal_orbit, radial)
    normal = _np.column_stack((radial, transverse, normal_orbit)) @ _np.array(
        [
            _np.cos(alpha),
            _np.sin(alpha) * _np.sin(beta),
            _np.sin(alpha) * _np.cos(beta),
        ]
    )

    reference = _np.array([0.0, 0.0, 1.0])
    if abs(normal @ reference) > 0.9:
        reference = _np.array([0.0, 1.0, 0.0])
    first = _np.cross(normal, reference)
    first /= _np.linalg.norm(first)
    second = _np.cross(normal, first)
    second /= _np.linalg.norm(second)
    corners = [
        cart[:3] + sail_size * (u * first + v * second)
        for u, v in ((1, 1), (-1, 1), (-1, -1), (1, -1))
    ]
    patch = Poly3DCollection(
        [corners], alpha=0.45, facecolor="silver", edgecolor="dimgray", linewidth=0.8
    )
    ax.add_collection3d(patch)
    ax.quiver(
        *cart[:3], *(normal * sail_size * 2), color="tab:blue", arrow_length_ratio=0.3
    )
    return patch


def plot_low_thrust_segments(
    leg,
    throttles,
    ax=None,
    N=30,
    mark_segments=True,
    mark_mismatch=True,
    mismatch_at_cut=False,
    state_to_cartesian=None,
    **kwargs,
):
    """Plot forward and backward ZOH segments colored by throttle.

    Args:
        leg (:class:`~pykep.leg.zoh` or :class:`~pykep.leg.zoh_ms`): ZOH leg to plot.
        throttles (:class:`array-like`): One normalized throttle value per segment.
        ax (:class:`matplotlib.axes.Axes`, optional): Three-dimensional axes.
            A new axes is created when omitted.
        N (:class:`int`, optional): Samples per segment. Default is 30.
        mark_segments (:class:`bool`, optional): Mark segment starts. Default is True.
        mark_mismatch (:class:`bool`, optional): Mark propagated segment ends.
            Default is True.
        mismatch_at_cut (:class:`bool`, optional): Mark only the forward/backward
            join instead of every segment end. Default is False.
        state_to_cartesian (:class:`callable`, optional): Convert each dynamics
            state to Cartesian coordinates before plotting.
        **kwargs: Additional arguments for segment markers.

    Returns:
        :class:`matplotlib.axes.Axes`: Axes containing the trajectory segments.
    """
    if ax is None:
        import pykep as _pk

        ax = _pk.plot.make_3Daxis()

    forward, backward, _ = leg.get_state_info(N=N)
    throttles = _np.asarray(throttles)
    groups = (
        (forward, range(leg.nseg_fwd)),
        (backward, range(leg.nseg - 1, leg.nseg_fwd - 1, -1)),
    )
    for segments, indices in groups:
        for segment, index in zip(segments, indices):
            if not len(segment):
                continue
            if state_to_cartesian is None:
                cart = _np.asarray(segment)
            else:
                cart = _np.asarray([state_to_cartesian(state) for state in segment])
            throttle = throttles[index]
            color = (
                0.25 + 0.55 * throttle,
                0.41 - 0.05 * throttle,
                0.88 - 0.52 * throttle,
            )
            ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], c=color)
            if mark_segments:
                ax.scatter(*cart[0, :3], **kwargs)
            if mark_mismatch and not mismatch_at_cut:
                ax.scatter(*cart[-1, :3], marker="^", **kwargs)
    if mark_mismatch and mismatch_at_cut:
        for segments in (forward, backward):
            nonempty = [segment for segment in segments if len(segment)]
            if nonempty:
                endpoint = state_to_cartesian(nonempty[-1][-1]) if state_to_cartesian else nonempty[-1][-1]
                ax.scatter(*endpoint[:3], marker="^", **kwargs)
    return ax


def plot_throttle_profile(tgrid, throttles, ax=None, **kwargs):
    """Plot a piecewise-constant throttle profile over a ZOH time grid.

    Args:
        tgrid (:class:`array-like`): Segment boundary times.
        throttles (:class:`array-like`): One normalized throttle value per segment.
        ax (:class:`matplotlib.axes.Axes`, optional): Axes for the plot. A new axes
            is created when omitted.
        **kwargs: Additional arguments passed to :meth:`matplotlib.axes.Axes.step`.

    Returns:
        :class:`matplotlib.axes.Axes`: Axes containing the throttle profile.
    """
    if ax is None:
        from matplotlib import pyplot as _plt

        _, ax = _plt.subplots()
    throttles = _np.asarray(throttles)
    ax.step(tgrid, _np.r_[throttles, throttles[-1]], where="post", **kwargs)
    ax.set_xlabel("time grid")
    ax.set_ylabel("throttle value (nd)")
    return ax