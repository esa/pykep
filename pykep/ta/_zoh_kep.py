import heyoka as _hy
import numpy as _np


def zoh_kep_dyn():
    """
    The dynamics in Cartesian coordinates of a constant thrust mass-varying spacecraft
    orbiting a main central body with unitary gravitational parameter :math:`\\mu = 1`.
    The electric propulsion is characterized by its effective
    exhaust velocity :math:`v_{eff} = I_{sp} g_0`.

    The dynamics are given by:

    .. math::
        \\left\\{
        \\begin{array}{l}
        \\dot{\\mathbf{r}} = \\mathbf{v} \\\\
        \\dot{\\mathbf{v}} = -\\frac{\\mathbf{r}}{r^3} + \\frac{T}{m} \\mathbf{\\hat{i}} \\\\
        \\dot{m} = -c T
        \\end{array}
        \\right.

    where :math:`c = \\frac{1}{v_{eff}}`, :math:`\\mathbf{\\hat{i}} = [i_x, i_y, i_z]`.

    The state is: :math:`[x, y, z, v_x, v_y, v_z, m]`
    
    The system parameters are (in this order): :math:`[T, i_x, i_y, i_z] + [c]`

    Returns:
        :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(p, dp), ...]
    """
    # Estalishing the state
    x, y, z, vx, vy, vz, m = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")

    # Naming the system parametes
    c = _hy.par[4]

    # Naming the system controls
    T = _hy.par[0]
    ix, iy, iz = _hy.par[1], _hy.par[2], _hy.par[3]

    # Auxiliary expressions
    r2 = _hy.sum([x**2, y**2, z**2])

    # Building the dynamics
    xdot = vx
    ydot = vy
    zdot = vz
    vxdot = -1.0 * pow(r2, -1.5) * x + T * ix / m
    vydot = -1.0 * pow(r2, -1.5) * y + T * iy / m
    vzdot = -1.0 * pow(r2, -1.5) * z + T * iz / m
    mdot = (
        -c * T * _hy.exp(-1.0 / m / 1e16)
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
_ta_zoh_kep_cache = dict()


def get_zoh_kep(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the :func:`~pykep.ta.zoh_kep_dyn` dynamics
    retrieving one from a global cache if available.

    This solves the initial value problem of a constant thrust mass-varying spacecraft.
    Thrust direction is fixed in the inertial frame.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
        Import and setup::

            import pykep as pk
            veff = 1.32
            controls = [0.022, 0.023, -0.21, 0.1]
            tof = 1.23

        Create propagator and propagate::

            ta = pk.ta.get_zoh_kep(tol=1e-16)
            ta.time = 0.
            ta.state[:] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
            ta.pars[:] = controls + [1 / veff]
            ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_kep_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zoh_kep_dyn(), init_state, tol=tol, pars=[1.0, 1.0, 0.0, 0.0, 0.0]
        )
        _ta_zoh_kep_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_kep_cache[tol]


_ta_zoh_kep_var_cache = dict()


def get_zoh_kep_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka)
    for the :func:`~pykep.ta.zoh_kep_dyn` dynamics retrieving one from
    a global cache if available.

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
            veff = 1.32
            controls = [0.022, 0.023, -0.21, 0.1]
            tof = 1.23

        Create propagator and propagate::

            ta = pk.ta.get_zoh_kep_var(tol=1e-16)
            ta.time = 0.
            ta.state[:7] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
            ta.pars[:] = controls + [1 / veff]
            ta.propagate_until(tof)
    """

    # Lookup.
    if _ta_zoh_kep_var_cache.get(tol) is None:
        # Cache miss, create new one.
        x, y, z, vx, vy, vz, m = _hy.make_vars("x", "y", "z", "vx", "vy", "vz", "m")
        vsys = _hy.var_ode_sys(
            zoh_kep_dyn(),
            [x, y, z, vx, vy, vz, m, _hy.par[0], _hy.par[1], _hy.par[2], _hy.par[3]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zoh_kep_var_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_kep_var_cache[tol]
