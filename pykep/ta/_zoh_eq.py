import heyoka as _hy
import numpy as _np


def zoh_eq_dyn():
    """
    The dynamics in equinoctial elements of a constant thrust (RTN frame) mass-varying spacecraft
    orbiting a main central body with unitary gravitational parameter :math:`\\mu = 1`.
    The electric propulsion is characterized by its effective
    exhaust velocity :math:`v_{eff} = I_{sp} g_0`.

    The dynamics are given by:

    .. math::
        \\left\\{
        \\begin{array}{l}
        \\dot p =  \\frac 1m\\sqrt{p} \\frac{2p}{w} T_t  \\\\
        \\dot f =  \\frac 1m\\sqrt{p} \\left\\{T_r \\sin L + \\left[ (1+w)\\cos L + f \\right] \\frac{T_t}{w} - (h\\sin L-k\\cos L)\\frac{g\\cdot T_n}{w} \\right\\}  \\\\
        \\dot g =  \\frac 1m\\sqrt{p} \\left\\{ - T_r\\cos L + \\left[ (1+w)\\sin L + g \\right] \\frac{T_t}{w} + (h\\sin L-k\\cos L)\\frac{f\\cdot T_n}{w} \\right\\}  \\\\
        \\dot h =  \\frac 1m\\sqrt{p} \\frac{s^2T_n}{2w}\\cos L  \\\\
        \\dot k =  \\frac 1m\\sqrt{p} \\frac{s^2T_n}{2w}\\sin L  \\\\
        \\dot L = \\sqrt{p}\\left\\{\\left(\\frac wp\\right)^2 + \\frac 1w\\left(h\\sin L-k\\cos L\\right)\\frac{ T_n}m \\right\\} \\\\
        \\dot m = - c T \\\\
        \\end{array}
        \\right.

    where :math:`c = \\frac{1}{v_{eff}}`, :math:`w = 1 + f\\cos L + g\\sin L`, :math:`s^2 = 1 + h^2 + k^2`.

    The state is: :math:`[p,f,g,h,k,L,m]`
    
    The system parameters are (in this order): :math:`[T, i_x, i_y, i_z] + [c]`

    In compact matrix form:

    .. math::
        \\dot{\\mathbf{x}} = \\mathbf{B}(\\mathbf{x}) \\frac{\\mathbf{T}}{m} + \\mathbf{D}(\\mathbf{x})

    with:

    .. math::
        \\sqrt{\\frac{p}{1}} \\mathbf{B}(\\mathbf{x}) = 
        \\begin{bmatrix}
        0 & \\frac {2p}w & 0 \\\\
        \\sin L & [(1+w)\\cos L + f]\\frac 1w  & - \\frac gw (h\\sin L-k\\cos L)  \\\\
        -\\cos L & [(1+w)\\sin L + g]\\frac 1w  & \\frac fw (h\\sin L-k\\cos L) \\\\
        0 & 0  & \\frac 1w \\frac{s^2}{2}\\cos L \\\\
        0 & 0  & \\frac 1w \\frac{s^2}{2}\\sin L \\\\
        0 & 0  & \\frac 1w (h\\sin L - k\\cos L) 
        \\end{bmatrix}

    and:

    .. math::
        \\mathbf{D}(\\mathbf{x}) = 
        \\begin{bmatrix}
        0 & 0 & 0 & 0 & 0 & \\sqrt{\\frac{1}{p^3}}w^2
        \\end{bmatrix}^T

    Returns:
        :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(p, dp), ...]
    """

    # The state
    p, f, g, h, k, L, m = _hy.make_vars("p", "f", "g", "h", "k", "L", "m")

    # Naming the system parametes (1/veff)
    c = _hy.par[4]

    # Naming the system controls
    T = _hy.par[0]
    i_r, i_t, i_n = _hy.par[1], _hy.par[2], _hy.par[3]

    w = 1 + f * _hy.cos(L) + g * _hy.sin(L)
    s2 = 1 + h * h + k * k
    B = _np.array(
        [
            [0, 2 * p / w, 0.0],
            [
                _hy.sin(L),
                ((1 + w) * _hy.cos(L) + f) / w,
                -g / w * (h * _hy.sin(L) - k * _hy.cos(L)),
            ],
            [
                -_hy.cos(L),
                ((1 + w) * _hy.sin(L) + g) / w,
                f / w * (h * _hy.sin(L) - k * _hy.cos(L)),
            ],
            [0, 0, s2 / w / 2.0 * _hy.cos(L)],
            [0, 0, s2 / w / 2.0 * _hy.sin(L)],
            [0, 0, 1.0 / w * (h * _hy.sin(L) - k * _hy.cos(L))],
        ]
    ) * _hy.sqrt(p)
    D = _np.array([0.0, 0.0, 0.0, 0.0, 0.0, _hy.sqrt(1.0 / p / p / p) * w * w])
    T_vector = _np.array([i_r, i_t, i_n]) * T

    # Dynamics
    fx = _np.dot(B, T_vector) / m + D
    fm = (
        -c * T * _hy.exp(-1.0 / m / 1e16)
    )  # the added term regularizes the dynamics keeping it differentiable

    # Assembling retval
    retval = []
    for var, f in zip([p, f, g, h, k, L], fx):
        retval.append((var, f))
    retval.append((m, fm))
    return retval


# We mimick in python the C++ global caching mechanism for taylor_adaptive
# instances, so that we do not create multiple instances with
# the same tolerance.
_ta_zoh_eq_cache = dict()


def get_zoh_eq(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the :func:`~pykep.ta.zoh_eq_dyn` dynamics
    retrieving one from a global cache if available.

    This solves the initial value problem of a constant thrust (RTN frame) mass-varying spacecraft.
    Thrust direction is fixed in the RTN frame.

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

            ta = pk.ta.get_zoh_eq(tol=1e-16)
            ta.time = 0.
            ta.state[:] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
            ta.pars[:] = controls + [1 / veff]
            ta.propagate_until(tof)
    """

    # Lookup.
    if _ta_zoh_eq_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zoh_eq_dyn(), init_state, tol=tol, pars=[1.0, 1.0, 0.0, 0.0, 0.0]
        )
        _ta_zoh_eq_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_eq_cache[tol]


_ta_zoh_eq_var_cache = dict()


def get_zoh_eq_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the
    :func:`~pykep.ta.zoh_eq_dyn` dynamics retrieving one from a global cache if available.

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

            ta = pk.ta.get_zoh_eq_var(tol=1e-16)
            ta.time = 0.
            ta.state[:] = [1., 0., 0., 0., 1., 0., 1.]
            ta.pars[:] = controls + [1 / veff]
            ta.propagate_until(tof)
    """

    # Lookup.
    if _ta_zoh_eq_var_cache.get(tol) is None:
        # Cache miss, create new one.
        [p, f, g, h, k, L, m] = _hy.make_vars("p", "f", "g", "h", "k", "L", "m")
        vsys = _hy.var_ode_sys(
            zoh_eq_dyn(),
            [p, f, g, h, k, L, m, _hy.par[0], _hy.par[1], _hy.par[2], _hy.par[3]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zoh_eq_var_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_eq_var_cache[tol]
