import heyoka as _hy
import numpy as _np


def zoh_ss_dyn():
    """
    The dynamics in Cartesian coordinates of a spacecraft equipped with a solar sail
    orbiting a main central body with unitary gravitational parameter :math:`\\mu = 1`.
    The body is assumed to be a star with flux :math:`C` at a reference distance `R = 1`.
    These conventions define naturally the nondimenional units of the system, which in the solar system,
    for example, are :math:`MU = \\mu_{Sun}` and :math:`L = 1 AU` (so that C is the flux at 1AU).
    
    The sail attitude is controlled through two angles (:math:`\\alpha`, :math:`\\beta`),
    defining the direction of the sail acceleration: they are commonly named cone and clock angle.
    The angles are assumed to be fixed, hence creating an acceleration having constant direction in
    the RTN reference frame :math:`[\\mathbf i_R, \\mathbf i_T, \\mathbf i_N]`.  

    The dynamics are given in non dimensional units (:matby:

    .. math::
        \\left\\{
        \\begin{array}{l}
        \\dot{\\mathbf{r}} = \\mathbf{v} \\\\
        \\dot{\\mathbf{v}} = -\\frac{\\mathbf{r}}{r^3} + \\mathbf a_{ss} \\\\
        \\end{array}
        \\right.

    where the acceleration coming from the solar sail is:
     
    .. math::
        \\mathbf a_{ss} = T \\left( \\cos\\alpha \\mathbf i_R + \\sin\\alpha\\sin\\beta \\mathbf i_T + \\sin\\alpha\\cos\\beta \\mathbf i_N\\right)
         
    and:
    
    .. math::
        T = c (1. / r^2) \\cos^2\\alpha

    The constant :math:`c` (in acceleration units) contains all physical parameters defining the sailcraft and is defined as
    :math:`c = 2 \\frac{CA}M`, where :math:`C` is central star flux at the reference unitary distance, A the sail area and M its mass.
    
    The cone angle must be constrained in :math:`\\alpha \\in [-\\frac{\\pi}2, \\frac{\\pi}2]` for the dynamics
    to be feasible (outside these bound the corresponding sail acceleration is unfeasible)

    The state is: :math:`[x, y, z, v_x, v_y, v_z]`
    
    The system parameters are (in this order): :math:`[\\alpha, \\beta] + [c]`

    .. note::
        with respect to other Zero Order Hold integrators the ``zoh_ss``has a smaller dimension (no mass equation) and only two
        controls. It CANNOT thus be used as Taylor integrator in :func:`~pykep.leg.zoh` as it does not meet its requirements.
        A dedicated leg, :func:`~pykep.leg.zoh_ss` is instead to be used.

    Returns:
        :class:`list` [ :class:`tuple` (:class:`hy::expression`, :class:`hy::expression` )]: The dynamics in the form [(p, dp), ...]
    """
    # The state
    x, y, z, vx, vy, vz = _hy.make_vars("x", "y", "z", "vx", "vy", "vz")
    # The controls
    alpha, beta = _hy.par[0], _hy.par[1]  # Sail angles
    # Auxiliary expressions
    r2 = x**2 + y**2 + z**2  # r^2
    h = _np.array(
        [y * vz - z * vy, z * vx - x * vz, x * vy - y * vx]
    )  # angular momentum vector
    ir = _np.array([x, y, z]) / _hy.sqrt(x**2 + y**2 + z**2)  # radius unit vector
    ih = h / _hy.sqrt(h[0] ** 2 + h[1] ** 2 + h[2] ** 2)  # angular momentum unit vector
    it = _np.array(
        [
            ih[1] * ir[2] - ih[2] * ir[1],
            ih[2] * ir[0] - ih[0] * ir[2],
            ih[0] * ir[1] - ih[1] * ir[0],
        ]
    )  # t unit vector
    T = (
        _hy.par[2] * (1.0 / r2) * _hy.cos(alpha) ** 2
    )  # sail acceleration magnitude 2C*A/M (C at 1 LU)
    # sail acceleration directions in the r,t,h frame
    ar = _hy.cos(alpha) * T
    at = _hy.sin(alpha) * _hy.sin(beta) * T
    ah = _hy.sin(alpha) * _hy.cos(beta) * T
    # dynamics
    dyn = (
        (x, vx),
        (y, vy),
        (z, vz),
        (vx, -1.0 * x / (r2 ** 1.5) + ar * ir[0] + at * it[0] + ah * ih[0]),
        (vy, -1.0 * y / (r2 ** 1.5) + ar * ir[1] + at * it[1] + ah * ih[1]),
        (vz, -1.0 * z / (r2 ** 1.5) + ar * ir[2] + at * it[2] + ah * ih[2]),
    )
    return dyn


# We mimick in python the C++ global caching mechanism for taylor_adaptive
# instances, so that we do not create multiple instances with
# the same tolerance.
_ta_zoh_ss_cache = dict()


def get_zoh_ss(tol: float):
    """
    Returns a Taylor adaptive propagator (Heyoka) for the :func:`~pykep.ta.zoh_ss_dyn` dynamics
    retrieving one from a global cache if available.

    This solves the initial value problem of a solar sailing spacecraft under simple hypothesis
    for the sail characteristics.

    Sail acceleration direction is fixed in the RTN frame.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive propagator.

    Examples:
        Import and setup::

            import pykep as pk
            sail_c = 0.01
            controls = [0.022, 0.023]
            tof = 1.23

        Create propagator and propagate::

            ta = pk.ta.get_zoh_ss(tol=1e-16)
            ta.time = 0.
            ta.state[:] = [1.2, 0.1, 0.1, 0.1, 1., 0.123, 1.]
            ta.pars[:] = controls + [sail_c]
            ta.propagate_until(tof)
    """
    # Lookup.
    if _ta_zoh_ss_cache.get(tol) is None:
        # Cache miss, create new one.
        init_state = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        new_ta = _hy.taylor_adaptive(
            zoh_ss_dyn(), init_state, tol=tol, pars=[0.0, 0.0, 0.0]
        )
        _ta_zoh_ss_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_ss_cache[tol]


_ta_zoh_ss_var_cache = dict()


def get_zoh_ss_var(tol: float):
    """Returns a (order 1) variational Taylor adaptive propagator (Heyoka)
    for the :func:`~pykep.ta.zoh_ss_dyn` dynamics retrieving one from
    a global cache if available.

    .. note::
       Variations are only considered with respect to initial conditions and the
       sail angles parameters :math:`[\\alpha, \\beta]`.

    Args:
        *tol* (:class:`float`): the tolerance of the Taylor adaptive variational propagator.

    Returns:
        :class:`hy::taylor_adaptive`: The Taylor adaptive variational propagator.

    Examples:
        Import and setup::

            import pykep as pk
            sail_c = 0.01
            controls = [0.022, 0.023]
            tof = 1.23

        Create propagator and propagate::

            ta = pk.ta.get_zoh_ss_var(tol=1e-16)
            ta.time = 0.
            ta.state[:7] = [1.2, 0.1, 0.1, 0.1, 1., 0.123]
            ta.pars[:] = controls + [sail_c]
            ta.propagate_until(tof)
    """

    # Lookup.
    if _ta_zoh_ss_var_cache.get(tol) is None:
        # Cache miss, create new one.
        x, y, z, vx, vy, vz = _hy.make_vars("x", "y", "z", "vx", "vy", "vz")
        vsys = _hy.var_ode_sys(
            zoh_ss_dyn(),
            [x, y, z, vx, vy, vz, _hy.par[0], _hy.par[1]],
            1,
        )
        new_ta = _hy.taylor_adaptive(
            vsys,
            tol=tol,
            compact_mode=True,
        )
        # Insert in cache and return
        _ta_zoh_ss_var_cache[tol] = new_ta
        return new_ta
    else:
        # Cache hit, return existing.
        return _ta_zoh_ss_var_cache[tol]
