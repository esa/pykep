import heyoka as _hy
import numpy as _np

def zoh_cr3bp_dyn():
    """
    The dynamics of a fixed thrust mass-varying spacecraft orbiting two main bodies
    in circular relative motion (circular restricted three-body problem, CR3BP).

    The dynamics are given by:

    .. math::
        \\left\\{
        \\begin{array}{l}
        \\dot{\\mathbf{r}} = \\mathbf{v} \\\\
        \\dot{v}_x = 2v_y + x - (1-\\mu)\\frac{x + \\mu}{r_1^3} - \\mu\\frac{x + \\mu - 1}{r_2^3} + T\\frac{i_x}{m} \\\\
        \\dot{v}_y = -2v_x + y - (1-\\mu)\\frac{y}{r_1^3} - \\mu\\frac{y}{r_2^3} + T\\frac{i_y}{m} \\\\
        \\dot{v}_z = -(1-\\mu)\\frac{z}{r_1^3} - \\mu\\frac{z}{r_2^3} + T\\frac{i_z}{m} \\\\
        \\dot{m} = -c T
        \\end{array}
        \\right.

    where:
    
    * :math:`r_1 = ||\\mathbf{r} - \\mathbf{r}_1||`, :math:`r_2 = ||\\mathbf{r} - \\mathbf{r}_2||`
    * :math:`\\mathbf{r}_1 = [-\\mu, 0, 0]`, :math:`\\mathbf{r}_2 = [1-\\mu, 0, 0]`
    * :math:`c = \\frac{1}{v_{eff}}`
    
    The state is: :math:`[x, y, z, v_x, v_y, v_z, m]`
    
    The system parameters are (in this order): :math:`[T, i_x, i_y, i_z] + [c, mu]`

    Returns:
        list[tuple[hy.expression, hy.expression]]: The dynamics in form `[(x, dx), ...]`
    """

    # The symbolic variables.
    [x, y, z, vx, vy, vz, m] = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")

    # Naming the system controls
    T_norm = _hy.par[0]
    i_x, i_y, i_z = _hy.par[1], _hy.par[2], _hy.par[3]

    # Naming the system parametes
    c = _hy.par[4]  # 1/veff
    mu = _hy.par[5]

    # Distances to the bodies.
    r_1 = _hy.sqrt(_hy.sum([pow(x + mu, 2.0), pow(y, 2.0), pow(z, 2.0)]))
    r_2 = _hy.sqrt(
        _hy.sum([pow(x - (1.0 - mu), 2.0), pow(y, 2.0), pow(z, 2.0)])
    )

    # The Equations of Motion.
    xdot = vx
    ydot = vy
    zdot = vz
    vxdot = (
        2.0 * vy
        + x
        - (1.0 - mu) * (x + mu) / (pow(r_1, 3.0))
        - mu * (x + mu - 1.0) / pow(r_2, 3.0)
        + T_norm * i_x / m
    )
    vydot = (
        -2.0 * vx
        + y
        - (1.0 - mu) * y / pow(r_1, 3.0)
        - mu * y / pow(r_2, 3.0)
        + T_norm * i_y / m
    )
    vzdot = (
        -(1.0 - mu) * z / pow(r_1, 3.0)
        - mu * z / pow(r_2, 3.0)
        + T_norm * i_z / m
    )
    mdot = (
        -c * T_norm * _hy.exp(-1.0 / m / 1e16)
    )  # the added term regularizes the dynamics keeping it differentiable

    retval = [
        (x, xdot),
        (y, ydot),
        (z, zdot),
        (vx, vxdot),
        (vy, vydot),
        (vz, vzdot),
        (m, mdot),
    ]
    return retval


# We mimick in python the C++ global caching mechanism for taylor_adaptive
# instances, so that we do not create multiple instances with
# the same tolerance.
_ta_zoh_cr3bp_cache = dict()


def get_zoh_cr3bp(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the :func:`~pykep.ta.zoh_cr3bp_dyn`
    dynamics retrieving one from a global cache if available.

    This solves the initial value problem of a constant thrust mass-varying spacecraft
    in the circular restricted three-body problem (CR3BP). Thrust direction is fixed
    in the rotating frame.

    If the requested propagator was never created, this creates it; otherwise returns cached version
    (avoiding re-jitting).

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
        Import and setup::

            import pykep as pk
            mu = pk.CR3BP_MU_EARTH_MOON
            veff = 1.
            controls = [0.001, 1., 0., 0.]
            tof = 1.

        Create propagator and propagate::

            ta = pk.ta.get_zoh_cr3bp(tol=1e-16)
            ta.time = 0.
            ta.state[:] = [1., 0., 0., 0., 1., 0., 1.]
            ta.pars[:] = controls + [1 / veff, mu]
            ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_cr3bp_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zoh_cr3bp_dyn(), init_state, tol=tol, pars=[1.0, 1.0, 0.0, 0.0, 0.0, 0.01]
        )
        _ta_zoh_cr3bp_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_cr3bp_cache[tol]


_ta_zoh_cr3bp_var_cache = dict()


def get_zoh_cr3bp_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka)
    for the :func:`~pykep.ta.zoh_cr3bp_dyn` dynamics retrieving one from
    a global cache if available.

    If the requested propagator was never created, this creates it; otherwise
    returns cached version (avoiding re-jitting).

    .. note::
       Variations are only considered with respect to initial conditions and the
       thrust parameters :math:`[T, i_x, i_y, i_z]`.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive variational propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive variational propagator.

    Examples:
        Import and setup::

            import pykep as pk
            mu = pk.CR3BP_MU_EARTH_MOON
            veff = 1.
            controls = [0.001, 1., 0., 0.]
            tof = 1.

        Create propagator and propagate::

            ta = pk.ta.get_zoh_kep_var(tol=1e-16)
            ta.time = 0.
            ta.state[:7] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
            ta.pars[:] = controls + [1 / veff, mu]
            ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_cr3bp_var_cache.get(tol) is None:
        # Cache miss, create new one.
        [x, y, z, vx, vy, vz, m] = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")
        vsys = _hy.var_ode_sys(
            zoh_cr3bp_dyn(),
            [x, y, z, vx, vy, vz, m, _hy.par[0], _hy.par[1], _hy.par[2], _hy.par[3]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zoh_cr3bp_var_cache[tol] = new_ta
        return (new_ta)
    else:
        # Cache hit, return existing.
        return _ta_zoh_cr3bp_var_cache[tol]
