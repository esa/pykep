import heyoka as _hy
import numpy as _np
from copy import deepcopy as _deepcopy

def zero_hold_eq_dyn():
    """
    The dynamics in equinoctial elements of a constant thrust (RTN frame) mass-varying spacecraft orbiting a main body.
    
    We consider the motion of a spacecraft of mass :math:`m` with a position :math:`\\mathbf{r}` and velocity :math:`\\mathbf{v}` subject only to the Sun's
    gravitational attraction in some inertial reference frame. The spacecraft also has an ion thruster with a specific impulse :math:`I_{sp}` thrusting as
    :math:`[T_r, T_t, T_n]`. We describe the spacecraft state via its mass :math:`m` and the (prograde) modified equinoctial elements
    :math:`\\mathbf{x}=\\left[p,f,g,h,k,L\\right]^T`. The dynamics are given by:
    
    .. math::
        \\begin{array}{l}
        \\dot p =  \\frac 1m\\sqrt{\\frac p\\mu} \\frac{2p}{w} T_t  \\\\
        \\dot f =  \\frac 1m\\sqrt{\\frac p\\mu} \\left\\{T_r \\sin L + \\left[ (1+w)\\cos L + f \\right] \\frac{T_t}{w} - (h\\sin L-k\\cos L)\\frac{g\\cdot T_n}{w} \\right\\}  \\\\
        \\dot g =  \\frac 1m\\sqrt{\\frac p\\mu} \\left\\{ - T_r\\cos L + \\left[ (1+w)\\sin L + g \\right] \\frac{T_t}{w} + (h\\sin L-k\\cos L)\\frac{f\\cdot T_n}{w} \\right\\}  \\\\
        \\dot h =  \\frac 1m\\sqrt{\\frac p\\mu} \\frac{s^2T_n}{2w}\\cos L  \\\\
        \\dot k =  \\frac 1m\\sqrt{\\frac p\\mu} \\frac{s^2T_n}{2w}\\sin L  \\\\
        \\dot L = \\sqrt{\\frac p\\mu}\\left\\{\\mu\\left(\\frac wp\\right)^2 + \\frac 1w\\left(h\\sin L-k\\cos L\\right)\\frac{ T_n}m \\right\\} \\\\
        \\dot m = - \\frac{|\\mathbf T|}{v_{eff}} \\\\
        \\end{array}
        
    where, :math:`w = 1 + f\\cos L + g\\sin L`, :math:`s^2 = 1 + h^2 + k^2` and :math:`[T_r, T_t, T_n]` are the radial, tangential and normal
    components of thrust direction. The gravitational parameter is denoted with :math:`\\mu` while :math:`v_{eff} = I_{sp}g_0` is 
    the effective velocity of the thruster. We rewrite the above equations in a more compact form:

    .. math::
        \\left\\{
        \\begin{array}{l}
        \\dot {\\mathbf x} =  \\mathbf B(\\mathbf x)  \\frac{\\mathbf{T}}m  + \\mathbf D(\\mathbf x) \\\\
        \\dot m = - \\frac{|\\mathbf T|}{v_{eff}} \\\\
        \\end{array}
        \\right.
        
    where:

    .. math::
        \\sqrt{\\frac \\mu p} \\mathbf B(\\mathbf x) = 
        \\left[
        \\begin{array}{ccc}
        0 & \\frac {2p}w & 0 \\
        \\sin L & [(1+w)\\cos L + f]\\frac 1w  & - \\frac gw (h\\sin L-k\\cos L)  \\
        - \\cos L & [(1+w)\\sin L + g]\\frac 1w  & \\frac fw (h\\sin L-k\\cos L) \\
        0 & 0  & \\frac 1w \\frac{s^2}{2}\\cos L \\
        0 & 0  & \\frac 1w \\frac{s^2}{2}\\sin L \\
        0 & 0  & \\frac 1w (h\\sin L - k\\cos L) \\
        \\end{array}
        \\right]

    and:

    .. math::
        \\mathbf D(\\mathbf x) = 
        \\left[
        \\begin{array}{cccccc}
        0 & 0 & 0 & 0 & 0 & \\sqrt{\\frac{\\mu}{p^3}}w^2
        \\end{array}
        \\right]^T
        
    and  :math:`[\\mu, v_{eff}, T_x, T_y, T_z]` are parameters.

    Returns:
        :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(p, dp), ...]

    """
    # The state
    p, f, g, h, k, L, m = _hy.make_vars("p", "f", "g", "h", "k", "L", "m")
    # The controls
    Tr, Tt, Tn = _hy.par[2], _hy.par[3], _hy.par[4]  # Thrust components
    mu, veff = (
        _hy.par[0],
        _hy.par[1],
    )  # Gravitational parameter and effective exhaust velocity
    # Useful expressions
    T_norm = _hy.sqrt(Tr * Tr + Tt * Tt + Tn * Tn)  # Thrust magnitude
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
    ) * _hy.sqrt(p / mu)
    D = _np.array([0.0, 0.0, 0.0, 0.0, 0.0, _hy.sqrt(mu / p / p / p) * w * w])
    T_vector = _np.array([Tr, Tt, Tn])

    # Dynamics
    fx = _np.dot(B, T_vector) / m + D
    fm = _hy.select(_hy.eq(T_norm, 0.), 0., -T_norm / veff)

    # Assembling retval
    retval = []
    for var, f in zip([p, f, g, h, k, L], fx):
        retval.append((var, f))
    retval.append((m, fm))
    return retval


# We mimick in python the C++ global caching mechanism for taylor_adaptive
# instances, so that we do not create multiple instances with
# the same tolerance.
_ta_zero_hold_eq_cache = dict()


def get_zero_hold_eq(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the zero_hold_eq dynamics retreiving
    one from a global cache and making a copy. 

    This is the initial value problem of a constant thrust (RTN frame) mass-varying spacecraft orbiting a primary. 
    Thrust direction is fixed in the RTN frame.

    If the requested propagator was never created this will create it, else it will
    return the one from the global cache, thus avoiding jitting.

    The dynamics is that returned by :func:`~pykep.ta.zero_hold_eq_dyn`.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
    >>> import pykep as pk
    >>> ta = pk.ta.get_zero_hold_eq(tol = 1e-16)
    >>> ta.time = 0.
    >>> ta.state[:] = [1.2,0.1,0.1,0.1,1.,0.123,1.]
    >>> mu = 1.122
    >>> veff = 1.32
    >>> thrust = [0.23, 0.21, 0.1]
    >>> tof = 1.23
    >>> ta.pars[:] = [mu, veff] + thrust
    >>> ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zero_hold_eq_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zero_hold_eq_dyn(), init_state, tol=tol, pars=[1.0, 1.0, 0.0, 0.0, 0.0]
        )
        _ta_zero_hold_eq_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zero_hold_eq_cache[tol]


_ta_zero_hold_eq_var_cache = dict()


def get_zero_hold_eq_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka) for the zero_hold_eq dynamics retreiving
    one from a global cache and making a copy. 

    This is the initial value problem of a constant thrust (RTN frame) mass-varying spacecraft orbiting a primary. 
    Thrust direction is fixed in the RTN frame.

    .. note:
    Variations are only considered with repsect to initial conditions and the constant thurst :math:`T_r, T_t, T_n`.
    The variations with respect to the constant thrust differ from those with respect to
    the often introduced throttles :math:`\\mathbf T = T_{max} \\mathbf u` by a factor :math:`T_{max}`.

    The dynamics is that returned by :func:`~pykep.ta.zero_hold_eq_dyn`: and also used in :func:`~pykep.ta.get_zero_hold_qe`

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator. 

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
    >>> import pykep as pk
    >>> ta = pk.ta.get_zero_hold_eq_var(tol = 1e-16)
    >>> ta.time = 0.
    >>> ta.state[:] = [1.,0.,0.,0.,1.,0.,1.]
    >>> mu = 1.
    >>> veff = 1.
    >>> thrust = [0., 0., 1e-22]
    >>> tof = 1.
    >>> ta.pars[:5] = [mu, veff] + thrust
    >>> ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zero_hold_eq_var_cache.get(tol) is None:
        # Cache miss, create new one.
        [p, f, g, h, k, L, m] = _hy.make_vars("p", "f", "g", "h", "k", "L", "m")
        vsys = _hy.var_ode_sys(
            zero_hold_eq_dyn(),
            [p, f, g, h, k, L, m, _hy.par[2], _hy.par[3], _hy.par[4]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zero_hold_eq_var_cache[tol] = new_ta
        return _deepcopy(new_ta)
    else:
        # Cache hit, return existing.
        return _deepcopy(_ta_zero_hold_eq_var_cache[tol])
