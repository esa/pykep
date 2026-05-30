import pykep as _pk
import numpy as _np
import heyoka as _hy
from matplotlib import pyplot as _plt


class tops_twobody:
    """
    Two-body problems from TOPS (Trajectory Optimisation Problems in Space).

    .. note::
        These instances have fixed end-points and (mostly) fixed time of flight,
        so they are not moving-end problems.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems described in json format. In this specific
    instance the spacecraft dynamics are modeled in Cartesian coordinates under
    two-body motion and. The transcription to a Non Linear Programming priblem (NLP) uses
    a forward-backward shooting scheme, which introduces highly nonlinear
    mismatch constraints, together with throttle constraints on the segment controls.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. The formulation is built on
    :class:`~pykep.trajopt.zoh_point2point`, so the departure and arrival
    states are fixed and moving-end effects are not accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a two-body TOPS instance from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (two-body Cartesian dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_twobody_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints and
                boundary relative-velocity direction constraints are formulated as inequalities
                (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_twobody_json[prob_name]
        self.name = "Two-body Keplerian: " + prob_name
        self.extra_info = gym_problem["info"]

        # ZOH integrators for two-body Cartesian dynamics
        tol = 1e-16
        tol_var = 1e-8

        ta_twobody = _pk.ta.get_zoh_kep(tol)
        ta_twobody_var = _pk.ta.get_zoh_kep_var(tol_var)

        # Units
        self.MU = gym_problem["mu"]
        self.L = _np.linalg.norm(gym_problem["state_s"][:3])
        self.MASS = _np.linalg.norm(gym_problem["m_s"])
        # The rest is induced
        self.TIME = _np.sqrt(self.L**3 / self.MU)
        self.V = self.L / self.TIME
        self.ACC = self.V / self.TIME
        self.F = self.MASS * self.ACC

        # Scaling the problem to these units (mu will be 1. -> as expected from the pykep ta_twobody)
        states = [it / self.L for it in gym_problem["state_s"][:3]] + [
            it / self.V for it in gym_problem["state_s"][3:6]
        ]
        statef = [it / self.L for it in gym_problem["state_f"][:3]] + [
            it / self.V for it in gym_problem["state_f"][3:6]
        ]
        ms = gym_problem["m_s"] / self.MASS
        max_thrust = gym_problem["max_thrust"] / self.F
        tof_bounds = [it / self.TIME for it in gym_problem["tof_bounds"]]

        # We set the Taylor integrators parameters (veff and mu in this case)
        veff_nd = gym_problem["veff"] / self.V
        ta_twobody.pars[4] = 1.0 / veff_nd
        ta_twobody_var.pars[4] = 1.0 / veff_nd

        self.udp = _pk.trajopt.zoh_point2point(
            states=states,  # 6D state only (mass handled separately)
            statef=statef,
            ms=ms,
            max_thrust=max_thrust,
            tof_bounds=tof_bounds,
            mf_bounds=[ms / 2.0, ms],
            nseg=nseg,
            cut=cut,
            tas=[ta_twobody, ta_twobody_var],
            time_encoding=time_encoding,
            state2cart=lambda x, L=self.L, V=self.V: _np.asarray(
                [x[0] * L, x[1] * L, x[2] * L, x[3] * V, x[4] * V, x[5] * V]
            ),
            inequalities_for_tc=inequalities_for_tc,
            max_steps=1000,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        add_sun=True,
        orbit_color="lightgray",
        figsize=(15, 4),
        **kwargs,
    ):
        """Four-view plot of the TOPS two-body solution."""
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)
        for axis in axes:
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, orbit_color=orbit_color, **kwargs
            )
        for axis in axes:
            axis.set_xticks([axis.get_xticks()[0], axis.get_xticks()[-1]])
            axis.set_yticks([axis.get_yticks()[0], axis.get_yticks()[-1]])
            axis.set_zticks([axis.get_zticks()[0], axis.get_zticks()[-1]])
        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])
        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        if add_sun:
            _pk.plot.add_sun(ax3d)
        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)
        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")
        return axes

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_twobody_mb:
    """
    Two-body problems from TOPS (Trajectory Optimisation Problems in Space).

    .. note::
        These instances have moving boundaries, as they allow the starting
        epoch to move by 1/4 of the period of the departure planet.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems described in json format. In this specific
    instance the spacecraft dynamics are modeled in Cartesian coordinates under
    two-body motion and. The transcription to a Non Linear Programming priblem (NLP) uses
    a forward-backward shooting scheme, which introduces highly nonlinear
    mismatch constraints, together with throttle constraints on the segment controls.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. The formulation is built on
    :class:`~pykep.trajopt.zoh_pl2pl`, so the departure and arrival
    states are moving boundaries and moving-end effects are accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a two-body TOPS instance with moving boundaries from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (two-body Cartesian dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_twobody_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints and
                boundary relative-velocity direction constraints are formulated as inequalities
                (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_twobody_json[prob_name]
        self.name = "Two-body Keplerian (Moving Boundaries): " + prob_name
        self.extra_info = gym_problem["info"]

        # These are the units declared in the original TOPS problem definition
        L = gym_problem["L"]
        TIME = gym_problem["TIME"]
        MASS = gym_problem["MASS"]
        V = L / TIME
        MU = L**3 / TIME**2
        ACC = V**2 / L
        F = MASS * ACC

        # We may thus instantiate the moving ends as keplerian planets
        # assuming the TOPs state vectors are at epoch 0 and and at epoch tof_days respectively.
        rs = [it * L for it in gym_problem["state_s"][0:3]]
        vs = [it * V for it in gym_problem["state_s"][3:6]]
        self.pls = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(0.0), posvel=[rs, vs], mu_central_body=MU
            )
        )
        tof_days = gym_problem["tof_bounds"][0] * TIME * _pk.SEC2DAY
        rf = [it * L for it in gym_problem["state_f"][0:3]]
        vf = [it * V for it in gym_problem["state_f"][3:6]]
        self.plf = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(tof_days), posvel=[rf, vf], mu_central_body=MU
            )
        )
        # To define a consistent version of the point w point problems as moving ends, we set the bounds
        # on the departure epochs to allow for a movement of half an orbit
        period_s_day = gym_problem["period_s"] * TIME * _pk.SEC2DAY
        bound = period_s_day / 4.0  # [-bound, bound] will be half

        # We define new units for the integrator so that the initial state has unitary norm (for consistency across problems).
        # This is not strictly necessary and introduces a third unit system (TOPS, integrator, and this) and likely confusion.
        # But its hidden in the definition of the UDP, and the user should not "see" it, unless they want to
        # actually recover the results in some of the other unit system (e.g. SI or TOPS) used in the problem definition.
        L_ta = _np.linalg.norm(gym_problem["state_s"][0:3]) * L
        MU_ta = gym_problem["mu"] * MU
        MASS_ta = gym_problem["m_s"] * MASS
        TIME_ta = _np.sqrt(L_ta**3 / MU_ta)
        V_ta = L_ta / TIME_ta
        ACC_ta = V_ta**2 / L_ta
        F_ta = MASS_ta * ACC_ta

        # Lower tolerances result in higher speed (the needed tolerance depends on the orbital regime)
        tol = 1e-16
        tol_var = 1e-8
        # We instantiate ZOH Taylor integrators for Keplerian dynamics.
        ta = _pk.ta.get_zoh_kep(tol)
        ta_var = _pk.ta.get_zoh_kep_var(tol_var)

        veff_nd = gym_problem["veff"] * V / V_ta
        ta.pars[4] = 1.0 / veff_nd
        ta_var.pars[4] = 1.0 / veff_nd

        self.udp = _pk.trajopt.zoh_pl2pl(
            pls=self.pls,
            plf=self.plf,
            ms=1.0,
            nseg=nseg,
            cut=cut,
            t0_bounds=[-bound, bound],  # (MJD2000)
            tof_bounds=[it * TIME / TIME_ta for it in gym_problem["tof_bounds"]],
            mf_bounds=[2 / 3, 1.0],
            vinf_dep_bounds=[0.0 * V / V_ta, 0.0 * V / V_ta],
            vinf_arr_bounds=[0.0 * V / V_ta, 0.0 * V / V_ta],
            tas=(ta, ta_var),
            max_thrust=gym_problem["max_thrust"] * F / F_ta,
            time_encoding=time_encoding,
            inequalities_for_tc=inequalities_for_tc,
            L=L_ta,
            V=V_ta,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        return self.name

    def get_extra_info(self):
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        add_sun=True,
        orbit_color="lightgray",
        planet_color="gray",
        figsize=(15, 4),
        **kwargs,
    ):
        """Plot a TOPS moving-boundaries solution.

        This convenience method reproduces the four-view layout used in the
        documentation notebook, adds the departure and arrival reference
        orbits, and plots the transfer on top. It returns the tuple of axes.
        """
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)

        t0 = self.udp.compute_t0(x_arr)
        tof_days = x_arr[10 + 4 * self.udp.nseg] * self.udp.TIME / _pk.DAY2SEC

        for axis in axes:
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.pls, label="pls", units=self.udp.L, color=orbit_color
            )
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.plf, label="plf", units=self.udp.L, color=orbit_color
            )
            _pk.plot.add_planet(
                ax=axis, when=t0, pla=self.pls, units=self.udp.L, color=planet_color, s=20
            )
            _pk.plot.add_planet(
                ax=axis, when=t0 + tof_days, pla=self.plf, units=self.udp.L, color=planet_color, s=20
            )
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, **kwargs
            )

        for axis in axes:
            xticks = axis.get_xticks()
            yticks = axis.get_yticks()
            zticks = axis.get_zticks()
            if len(xticks) >= 2:
                axis.set_xticks([xticks[0], xticks[-1]])
            if len(yticks) >= 2:
                axis.set_yticks([yticks[0], yticks[-1]])
            if len(zticks) >= 2:
                axis.set_zticks([zticks[0], zticks[-1]])

        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])

        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        if add_sun:
            _pk.plot.add_sun(ax3d)

        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)

        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")

        return axes

    def __getattr__(self, name):
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_mee:
    """
    Mean equinoctial elements problems from TOPS.

    .. note::
        These instances have fixed end-points and (mostly) fixed time of flight,
        so they are not moving-end problems.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft dynamics are modeled in modified equinoctial
    elements under two-body motion. The transcription uses a forward-backward
    shooting scheme, which introduces highly nonlinear mismatch constraints,
    together with throttle constraints on the segment controls.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. The formulation is built on
    :class:`~pykep.trajopt.zoh_point2point`, so the departure and arrival
    states are fixed and moving-end effects are not accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a two-body MEE TOPS instance from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (two-body dynamics in modified equinoctial elements). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_mee_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints are
                formulated as inequalities (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_mee_json[prob_name]
        self.name = "Two-body MEE: " + prob_name
        self.extra_info = gym_problem["info"]

        # Instantiating the ZOH integrators
        # Tolerances used in the numerical integration
        # Low tolerances result in higher speed (the needed tolerance depends on the orbit)
        tol = 1e-16
        tol_var = 1e-8

        # We instantiate ZOH Taylor integrators for the MEE dynamics.
        ta_twobody_mee = _pk.ta.get_zoh_eq(tol)
        ta_twobody_mee_var = _pk.ta.get_zoh_eq_var(tol_var)

        # Units
        self.MU = gym_problem["mu"]
        self.L = gym_problem["state_s"][0]
        self.MASS = gym_problem["m_s"]
        # The rest is induced
        self.TIME = _np.sqrt(self.L**3 / self.MU)
        self.V = self.L / self.TIME
        self.ACC = self.V / self.TIME
        self.F = self.MASS * self.ACC

        # Scaling the problem to these units (mu will be 1. -> as expected from pykep ta_twobody_mee)
        states = [it / self.L for it in gym_problem["state_s"][:1]] + [
            it for it in gym_problem["state_s"][1:6]
        ]
        statef = [it / self.L for it in gym_problem["state_f"][:1]] + [
            it for it in gym_problem["state_f"][1:6]
        ]
        ms = gym_problem["m_s"] / self.MASS
        max_thrust = gym_problem["max_thrust"] / self.F
        tof_bounds = [it / self.TIME for it in gym_problem["tof_bounds"]]

        # We set the Taylor integrators parameters (veff and mu in this case)
        veff_nd = gym_problem["veff"] / self.V
        ta_twobody_mee.pars[4] = 1.0 / veff_nd
        ta_twobody_mee_var.pars[4] = 1.0 / veff_nd

        self.udp = _pk.trajopt.zoh_point2point(
            states=states,
            statef=statef,
            ms=ms,
            max_thrust=max_thrust,
            tof_bounds=tof_bounds,
            mf_bounds=[ms / 2.0, ms],
            nseg=nseg,
            cut=cut,
            tas=[ta_twobody_mee, ta_twobody_mee_var],
            time_encoding=time_encoding,
            state2cart=lambda x, L=self.L: _np.array(
                _pk.mee2ic([x[0] * L, x[1], x[2], x[3], x[4], x[5]], mu=1.0)
            ).flatten(),
            inequalities_for_tc=inequalities_for_tc,
            max_steps=1000,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        add_sun=True,
        orbit_color="lightgray",
        figsize=(15, 4),
        **kwargs,
    ):
        """Four-view plot of the TOPS MEE solution."""
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)
        for axis in axes:
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, orbit_color=orbit_color, **kwargs
            )
        for axis in axes:
            axis.set_xticks([axis.get_xticks()[0], axis.get_xticks()[-1]])
            axis.set_yticks([axis.get_yticks()[0], axis.get_yticks()[-1]])
            axis.set_zticks([axis.get_zticks()[0], axis.get_zticks()[-1]])
        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])
        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        if add_sun:
            _pk.plot.add_sun(ax3d)
        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)
        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")
        return axes

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_mee_mb:
    """
    Mean equinoctial elements problems from TOPS.

    .. note::
        These instances have moving boundaries, as they allow the starting
        epoch to move by 1/4 of the period of the departure planet.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft dynamics are modeled in modified equinoctial
    elements under two-body motion. The transcription uses a forward-backward
    shooting scheme, which introduces highly nonlinear mismatch constraints,
    together with throttle constraints on the segment controls.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. The formulation is built on
    :class:`~pykep.trajopt.zoh_pl2pl`, so the departure and arrival
    states are moving boundaries and moving-end effects are accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a two-body MEE TOPS instance with moving boundaries from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (two-body dynamics in modified equinoctial elements). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_mee_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints are
                formulated as inequalities (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_mee_json[prob_name]
        self.name = "Two-body MEE (Moving Boundaries): " + prob_name
        self.extra_info = gym_problem["info"]

        # Since these are mee problems we compute the number of revolutions as the difference in mean longitude divided by 2pi
        nrevs = (gym_problem["state_f"][-1] - gym_problem["state_s"][-1]) // (
            _np.pi * 2
        )

        # We define the units declared in the TOPS problem definition
        L = gym_problem["L"]
        TIME = gym_problem["TIME"]
        MASS = gym_problem["MASS"]
        V = L / TIME
        MU = L**3 / TIME**2
        ACC = V**2 / L
        F = MASS * ACC

        # And instantiate the moving ends as planets
        # We scale the non dimensional MEE state by the units to get physical units,
        # then we create the planet objects for the initial and final states.
        mees = gym_problem["state_s"].copy()
        mees[0] *= gym_problem["L"]
        self.pls = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(0.0),
                elem=mees,
                mu_central_body=gym_problem["mu"] * MU,
                el_type=_pk.el_type.MEE,
            )
        )
        # Same for the final state
        tof_days = gym_problem["tof_bounds"][0] * TIME * _pk.SEC2DAY
        meef = gym_problem["state_f"].copy()
        meef[0] *= gym_problem["L"]
        self.plf = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(tof_days),
                elem=meef,
                mu_central_body=gym_problem["mu"] * MU,
                el_type=_pk.el_type.MEE,
            )
        )

        # To let boundary move we set the bounds on the departure epoch to allow for movement of hal an orbit
        period_s_day = gym_problem["period_s"] * TIME * _pk.SEC2DAY
        bound = period_s_day / 4.0  # [-bound, bound] will be half

        # We rescale the integrator units so that the initial state has unitary norm (for consistency across problems).
        # This is not strictly necessary and introduces a third unit system (TOPS, integrator, and this) and likely confusion.
        # But its hidden in the definition of the UDP, and the user should not "see" it, unless they want to
        # actually recover the results in some of the other unit system (e.g. SI or TOPS) used in the problem definition.
        L_ta = _np.linalg.norm(gym_problem["state_s"][0]) * L
        MU_ta = gym_problem["mu"] * MU
        MASS_ta = gym_problem["m_s"] * MASS
        TIME_ta = _np.sqrt(L_ta**3 / MU_ta)
        V_ta = L_ta / TIME_ta
        ACC_ta = V_ta**2 / L_ta
        F_ta = MASS_ta * ACC_ta

        # Lower tolerances result in higher speed (the needed tolerance depends on the orbital regime)
        tol = 1e-16
        tol_var = 1e-8
        # We instantiate ZOH Taylor integrators for mee dynamics.
        ta = _pk.ta.get_zoh_eq(tol)
        ta_var = _pk.ta.get_zoh_eq_var(tol_var)

        # We set the effective exhaust velocity in the Taylor integrators
        veff_nd = gym_problem["veff"] * V / V_ta
        ta.pars[4] = 1.0 / veff_nd
        ta_var.pars[4] = 1.0 / veff_nd

        # We get the symbolic expressions
        cart2state, cart2state_J = _pk.ic2mee(
            jacobian=True
        )  # The API also requires the Jacobian
        state2cart, _ = _pk.mee2ic(
            jacobian=False
        )  # The API does not require the Jacobian
        # We jit them into cfuncs
        x, y, z, vx, vy, vz = _hy.make_vars("x", "y", "z", "vx", "vy", "vz")
        p, f, g, h, k, L_anom = _hy.make_vars("p", "f", "g", "h", "k", "L")
        cart2state_cfunc = _hy.cfunc(cart2state, vars=[x, y, z, vx, vy, vz])
        cart2state_J_cfunc = _hy.cfunc(cart2state_J, vars=[x, y, z, vx, vy, vz])
        state2cart_cfunc = _hy.cfunc(state2cart, vars=[p, f, g, h, k, L_anom])
        # We create the call signature required by the API of the zoh_pl2pl class and
        # account here for the nrevs in the mean longitude.
        cart2state_zoh = lambda x: cart2state_cfunc(x, pars=[1, 1])
        cart2state_J_zoh = lambda x: cart2state_J_cfunc(x, pars=[gym_problem["mu"], 1])
        state2cart_zoh = lambda x: state2cart_cfunc(x, pars=[1, 1])

        self.udp = _pk.trajopt.zoh_pl2pl(
            pls=self.pls,
            plf=self.plf,
            ms=gym_problem["m_s"] * MASS / MASS_ta,
            nseg=nseg,
            cut=cut,
            t0_bounds=[-0, 0],  # (MJD2000)
            tof_bounds=[it * TIME / TIME_ta for it in gym_problem["tof_bounds"]],
            mf_bounds=[0.5, 1.0],
            vinf_dep_bounds=[0.0 * V / V_ta, 0.0 * V / V_ta],
            vinf_arr_bounds=[0.0 * V / V_ta, 0.0 * V / V_ta],
            tas=(ta, ta_var),
            max_thrust=gym_problem["max_thrust"] * F / F_ta,
            time_encoding=time_encoding,
            inequalities_for_tc=inequalities_for_tc,
            cart2state=(cart2state_zoh, cart2state_J_zoh),
            state2cart=state2cart_zoh,
            L=L_ta,
            V=V_ta,
            nrevs=nrevs,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        return self.name

    def get_extra_info(self):
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        add_sun=True,
        orbit_color="lightgray",
        planet_color="gray",
        figsize=(15, 4),
        **kwargs,
    ):
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)
        t0 = self.udp.compute_t0(x_arr)
        tof_days = x_arr[10 + 4 * self.udp.nseg] * self.udp.TIME / _pk.DAY2SEC
        for axis in axes:
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.pls, label="pls", units=self.udp.L, color=orbit_color
            )
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.plf, label="plf", units=self.udp.L, color=orbit_color
            )
            _pk.plot.add_planet(
                ax=axis, when=t0, pla=self.pls, units=self.udp.L, color=planet_color, s=20
            )
            _pk.plot.add_planet(
                ax=axis, when=t0 + tof_days, pla=self.plf, units=self.udp.L, color=planet_color, s=20
            )
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, **kwargs
            )
        for axis in axes:
            xticks = axis.get_xticks()
            yticks = axis.get_yticks()
            zticks = axis.get_zticks()
            if len(xticks) >= 2:
                axis.set_xticks([xticks[0], xticks[-1]])
            if len(yticks) >= 2:
                axis.set_yticks([yticks[0], yticks[-1]])
            if len(zticks) >= 2:
                axis.set_zticks([zticks[0], zticks[-1]])
        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])
        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        if add_sun:
            _pk.plot.add_sun(ax3d)
        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)
        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")
        return axes

    def __getattr__(self, name):
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_ss:
    """
    Solar sailing problems from TOPS.

    .. note::
        These instances have fixed end-points and (mostly) fixed time of flight,
        so they are not moving-end problems.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft is controlled by a solar sail and the transcription
    uses a forward-backward shooting scheme, which introduces highly nonlinear
    mismatch constraints.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. Unlike the non solar-sail
    TOPS instances, this formulation does not include throttle constraints. The
    formulation is built on :class:`~pykep.trajopt.zoh_ss_point2point`, so the
    departure and arrival states are fixed and moving-end effects are not
    accounted for.
    """

    def __init__(self, prob_name, cut=0.5, nseg=10, time_encoding="uniform"):
        """Construct a solar-sailing TOPS instance from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (solar sail dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_ss_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        gym_problem = _pk.trajopt.gym.tops_ss_json[prob_name]
        self.name = "Solar sailing: " + prob_name
        self.extra_info = gym_problem["info"]

        # Instantiating the ZOH integrators
        # Tolerances used in the numerical integration
        # Low tolerances result in higher speed (the needed tolerance depends on the orbit)
        tol = 1e-16
        tol_var = 1e-8

        # We instantiate ZOH Taylor integrators for the SOLARSAIL dynamics.
        ta_ss = _pk.ta.get_zoh_ss(tol)
        ta_ss_var = _pk.ta.get_zoh_ss_var(tol_var)

        # Units
        self.MU = gym_problem["mu"]
        self.L = _np.linalg.norm(gym_problem["state_s"][:3])
        # The rest is induced
        self.TIME = _np.sqrt(self.L**3 / self.MU)
        self.V = self.L / self.TIME
        self.ACC = self.V / self.TIME

        # Scaling the problem to these units (mu will be 1. -> as expected from pykep ta_ss)
        states = [it / self.L for it in gym_problem["state_s"][:3]] + [
            it / self.V for it in gym_problem["state_s"][3:6]
        ]
        statef = [it / self.L for it in gym_problem["state_f"][:3]] + [
            it / self.V for it in gym_problem["state_f"][3:6]
        ]
        c_sail = gym_problem["c_sail"] / self.ACC
        tof_bounds = [it / self.TIME for it in gym_problem["tof_bounds"]]

        # We set the Taylor integrators parameter (sail acceleration constant in this case)
        ta_ss.pars[2] = c_sail
        ta_ss_var.pars[2] = c_sail

        self.udp = _pk.trajopt.zoh_ss_point2point(
            states=states,
            statef=statef,
            tof_bounds=tof_bounds,
            nseg=nseg,
            cut=cut,
            tas=[ta_ss, ta_ss_var],
            time_encoding=time_encoding,
            max_steps=1000,
        )

    def get_nic(self):
        """Override: solar_sailing problems have no inequality constraints."""
        return 0

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding))

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info

    def plot(self, x, N=100, mark_segments=True, figsize=(6, 6), **kwargs):
        """Single-view plot of the TOPS solar sail solution."""
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        self.udp.plot(
            x=x_arr, ax=ax, N=N, mark_segments=mark_segments, **kwargs
        )
        ax.set_xticks([ax.get_xticks()[0], ax.get_xticks()[-1]])
        ax.set_yticks([ax.get_yticks()[0], ax.get_yticks()[-1]])
        ax.set_zticks([ax.get_zticks()[0], ax.get_zticks()[-1]])
        ax.view_init(elev=45, azim=270)
        ax.set_title("3D")
        return ax

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_ss_mb:
    """
    Solar sailing problems from TOPS.

    .. note::
        These instances have moving boundaries, as they allow the starting
        epoch to move by 1/4 of the period of the departure planet.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft is controlled by a solar sail and the transcription
    uses a forward-backward shooting scheme, which introduces highly nonlinear
    mismatch constraints.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. Unlike the non solar-sail
    TOPS instances, this formulation does not include throttle constraints. The
    formulation is built on :class:`~pykep.trajopt.zoh_ss_pl2pl`, so the
    departure and arrival states are moving boundaries and moving-end effects are
    accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a solar-sailing TOPS instance with moving boundaries from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (solar sail dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_ss_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, boundary relative-velocity
                direction constraints are formulated as inequalities, otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_ss_json[prob_name]
        self.name = "Solar Sailing (Moving Boundaries): " + prob_name
        self.extra_info = gym_problem["info"]

        # These are the units declared in the original TOPS problem definition
        L = gym_problem["L"]
        TIME = gym_problem["TIME"]
        V = L / TIME
        MU = L**3 / TIME**2
        ACC = V**2 / L

        # We may thus instantiate the moving ends as keplerian planets
        # assuming the TOPs state vectors are at epoch 0 and and at epoch tof_days respectively.
        rs = [it * L for it in gym_problem["state_s"][0:3]]
        vs = [it * V for it in gym_problem["state_s"][3:6]]
        self.pls = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(0.0), posvel=[rs, vs], mu_central_body=MU
            )
        )
        tof_days = gym_problem["tof_bounds"][0] * TIME * _pk.SEC2DAY
        rf = [it * L for it in gym_problem["state_f"][0:3]]
        vf = [it * V for it in gym_problem["state_f"][3:6]]
        self.plf = _pk.planet(
            _pk.udpla.keplerian(
                when=_pk.epoch(tof_days), posvel=[rf, vf], mu_central_body=MU
            )
        )
        # To define a consistent version of the point w point problems as moving ends, we set the bounds
        # on the departure epochs to allow for a movement of half an orbit
        period_s_day = gym_problem["period_s"] * TIME * _pk.SEC2DAY
        bound = period_s_day / 4.0  # [-bound, bound] will be half

        # Lower tolerances result in higher speed (the needed tolerance depends on the orbital regime)
        tol = 1e-16
        tol_var = 1e-8
        # We instantiate ZOH Taylor integrators for Keplerian dynamics.
        ta = _pk.ta.get_zoh_ss(tol)
        ta_var = _pk.ta.get_zoh_ss_var(tol_var)

        # We set the Taylor integrators parameter (sail acceleration constant in this case)
        ta.pars[2] = gym_problem["c_sail"]
        ta_var.pars[2] = gym_problem["c_sail"]

        self.udp = _pk.trajopt.zoh_ss_pl2pl(
            pls=self.pls,
            plf=self.plf,
            nseg=nseg,
            cut=cut,
            t0_bounds=[-bound, bound],
            tof_bounds=gym_problem["tof_bounds"],
            vinf_dep_bounds=[0.0, 0.0],
            vinf_arr_bounds=[0.0, 0.0],
            tas=(ta, ta_var),
            time_encoding=time_encoding,
            inequalities_for_tc=inequalities_for_tc,
            L=L,
            V=V,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        return self.name

    def get_extra_info(self):
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        add_sun=True,
        orbit_color="lightgray",
        planet_color="gray",
        figsize=(6, 6),
        **kwargs,
    ):
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        t0 = self.udp.compute_t0(x_arr)
        tof_days = x_arr[9 + 2 * self.udp.nseg] * self.udp.TIME / _pk.DAY2SEC
        _pk.plot.add_planet_orbit(
            ax=ax, pla=self.pls, label="pls", units=self.udp.L, color=orbit_color
        )
        _pk.plot.add_planet_orbit(
            ax=ax, pla=self.plf, label="plf", units=self.udp.L, color=orbit_color
        )
        _pk.plot.add_planet(
            ax=ax, when=t0, pla=self.pls, units=self.udp.L, color=planet_color, s=20
        )
        _pk.plot.add_planet(
            ax=ax, when=t0 + tof_days, pla=self.plf, units=self.udp.L, color=planet_color, s=20
        )
        self.udp.plot(
            x=x_arr, ax=ax, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, **kwargs
        )
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()
        zticks = ax.get_zticks()
        if len(xticks) >= 2:
            ax.set_xticks([xticks[0], xticks[-1]])
        if len(yticks) >= 2:
            ax.set_yticks([yticks[0], yticks[-1]])
        if len(zticks) >= 2:
            ax.set_zticks([zticks[0], zticks[-1]])
        ax.view_init(elev=45, azim=270)
        ax.set_aspect("equal")
        if add_sun:
            _pk.plot.add_sun(ax)
        ax.set_title("3D")
        return ax

    def __getattr__(self, name):
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_cr3bp:
    """
    CR3BP problems from TOPS.

    .. note::
        These instances have fixed end-points and (mostly) fixed time of flight,
        so they are not moving-end problems.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft dynamics are modeled in the circular restricted
    three-body problem. The transcription uses a forward-backward shooting
    scheme, which introduces highly nonlinear mismatch constraints.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. In this CR3BP formulation
    throttle constraints are not exposed as separate inequalities. The
    formulation is built on :class:`~pykep.trajopt.zoh_point2point`, so the
    departure and arrival states are fixed and moving-end effects are not
    accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=60,
        time_encoding="uniform",
        inequalities_for_tc=False,
    ):
        """Construct a CR3BP TOPS instance from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (circular restricted three-body problem dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_cr3bp_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. CR3BP problems typically require
                more segments than two-body problems. Defaults to 60.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints are
                formulated as inequalities (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to False.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_cr3bp_json[prob_name]
        self.name = "CR3BP: " + prob_name
        self.extra_info = gym_problem["info"]

        # Instantiating the ZOH integrators
        # Tolerances used in the numerical integration
        # Low tolerances result in higher speed (the needed tolerance depends on the orbit)
        tol = 1e-16
        tol_var = 1e-8

        # For consistency we provide units (these are not the final ones rather the ones
        # transforming data to the usual cr3bp units)
        self.L = 1.0
        self.V = 1.0
        self.TIME = self.L / self.V
        self.ACC = self.V / self.TIME

        # We instantiate ZOH Taylor integrators for the CR3BP dynamics.
        ta_cr3bp = _pk.ta.get_zoh_cr3bp(tol)
        ta_cr3bp_var = _pk.ta.get_zoh_cr3bp_var(tol_var)

        # We set the Taylor integrators parameters (veff and mu in this case)
        veff_nd = gym_problem["veff"]
        mu_cr3bp = gym_problem["mu_cr3bp"]
        ta_cr3bp.pars[4] = 1.0 / veff_nd
        ta_cr3bp_var.pars[4] = 1.0 / veff_nd
        ta_cr3bp.pars[5] = mu_cr3bp
        ta_cr3bp_var.pars[5] = mu_cr3bp

        self.udp = _pk.trajopt.zoh_point2point(
            states=gym_problem["state_s"],
            statef=gym_problem["state_f"],
            ms=gym_problem["m_s"],
            max_thrust=gym_problem["max_thrust"],
            tof_bounds=gym_problem["tof_bounds"],
            mf_bounds=[gym_problem["m_s"] / 2.0, gym_problem["m_s"]],
            nseg=nseg,
            cut=cut,
            tas=[ta_cr3bp, ta_cr3bp_var],
            time_encoding=time_encoding,
            inequalities_for_tc=inequalities_for_tc,
            max_steps=1000,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        orbit_color="lightgray",
        figsize=(15, 4),
        **kwargs,
    ):
        """Four-view plot of the TOPS CR3BP solution."""
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)
        for axis in axes:
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, orbit_color=orbit_color, **kwargs
            )
        for axis in axes:
            axis.set_xticks([axis.get_xticks()[0], axis.get_xticks()[-1]])
            axis.set_yticks([axis.get_yticks()[0], axis.get_yticks()[-1]])
            axis.set_zticks([axis.get_zticks()[0], axis.get_zticks()[-1]])
        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])
        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)
        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")
        return axes

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)


class tops_cr3bp_mb:
    """
    CR3BP problems from TOPS.

    .. note::
        These instances have moving boundaries, as they allow the starting
        epoch to move by 1/4 of the period of the departure planet.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft dynamics are modeled in the circular restricted
    three-body problem. The transcription uses a forward-backward shooting
    scheme, which introduces highly nonlinear mismatch constraints.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. In this CR3BP formulation
    throttle constraints are not exposed as separate inequalities. The
    formulation is built on :class:`~pykep.trajopt.zoh_pl2pl`, so the
    departure and arrival states are moving boundaries and moving-end effects are
    accounted for.
    """

    def __init__(
        self,
        prob_name,
        cut=0.5,
        nseg=10,
        time_encoding="uniform",
        inequalities_for_tc=True,
    ):
        """Construct a CR3BP TOPS instance with moving boundaries from a predefined gym case.

        Args:
            *prob_name* (:class:`str`): Name of the predefined gym problem case from the TOPS database
                (circular restricted three-body problem dynamics). Available problem names can be found in the keys of
                :data:`~pykep.trajopt.gym.tops_cr3bp_json`. For example: ``'P0'``, ``'P1'``, ``'P2'``, etc.

            *cut* (:class:`float`, optional): Cut parameter for the forward/backward split in the ZOH transcription.
                Defaults to 0.5.

            *nseg* (:class:`int`, optional): Number of constant-control segments for the zero-order hold transcription.
                Each segment contributes to the NLP dimensionality and complexity. CR3BP problems typically require
                more segments than two-body problems. Defaults to 10.

            *time_encoding* (:class:`str`, optional): Time-grid encoding scheme. Options are:

                - ``'uniform'``: equally spaced segment boundaries over the total time of flight.
                - ``'softmax'``: variable segment lengths controlled by softmax weights added to the decision vector.

                Defaults to ``'uniform'``.

            *inequalities_for_tc* (:class:`bool`, optional): If True, throttle constraints and
                boundary relative-velocity direction constraints are formulated as inequalities
                (``|i_u|**2 - 1 <= 0``), otherwise as equalities.
                Defaults to True.
        """
        self.prob_name = prob_name
        self.cut = cut
        self.nseg = nseg
        self.time_encoding = time_encoding
        self.inequalities_for_tc = inequalities_for_tc
        gym_problem = _pk.trajopt.gym.tops_cr3bp_json[prob_name]
        self.name = "CR3BP (Moving Boundaries): " + prob_name
        self.extra_info = gym_problem["info"]

        # Units
        L = gym_problem["L"]
        TIME = gym_problem["TIME"]
        MASS = gym_problem["MASS"]
        V = L / TIME
        MU = L**3 / TIME**2
        ACC = V**2 / L
        F = MASS * ACC

        # We now need to define planets in the CR3BP
        self.pls = _pk.planet(
            _pk.udpla.cr3bp(
                when=_pk.epoch(0.0),
                state_nd=gym_problem["state_s"],
                mu_cr3bp=gym_problem["mu_cr3bp"],
                TIME=TIME,
                L=L,
                name="pls",
                tol=1e-16,
            )
        )
        tof_days = gym_problem["tof_bounds"][0] * TIME * _pk.SEC2DAY
        self.plf = _pk.planet(
            _pk.udpla.cr3bp(
                when=_pk.epoch(tof_days),
                state_nd=gym_problem["state_f"],
                mu_cr3bp=gym_problem["mu_cr3bp"],
                TIME=TIME,
                L=L,
                name="plf",
                tol=1e-16,
            )
        )

        # To let boundary move we set the bounds on the departure epoch to allow for movement of hal an orbit
        period_s_day = gym_problem["period_s"] * TIME * _pk.SEC2DAY
        bound = period_s_day / 4.0  # [-bound, bound] will be half

        # Lower tolerances result in higher speed (the needed tolerance depends on the orbital regime)
        tol = 1e-16
        tol_var = 1e-8
        # We instantiate ZOH Taylor integrators for mee dynamics.
        ta = _pk.ta.get_zoh_cr3bp(tol)
        ta_var = _pk.ta.get_zoh_cr3bp_var(tol_var)

        # We set the effective exhaust velocity in the Taylor integrators
        veff_nd = gym_problem["veff"]
        mu_cr3bp = gym_problem["mu_cr3bp"]
        ta.pars[4] = 1.0 / veff_nd
        ta_var.pars[4] = 1.0 / veff_nd
        ta.pars[5] = mu_cr3bp
        ta_var.pars[5] = mu_cr3bp

        self.udp = _pk.trajopt.zoh_pl2pl(
            pls=self.pls,
            plf=self.plf,
            ms=gym_problem["m_s"],
            nseg=nseg,
            cut=cut,
            t0_bounds=[-bound, bound],
            tof_bounds=gym_problem["tof_bounds"],
            mf_bounds=[gym_problem["m_s"] / 2.0, gym_problem["m_s"]],
            vinf_dep_bounds=[0.0, 0.0],
            vinf_arr_bounds=[0.0, 0.0],
            tas=(ta, ta_var),
            max_thrust=gym_problem["max_thrust"],
            time_encoding="uniform",
            inequalities_for_tc=inequalities_for_tc,
            L=L,
            V=V,
        )

    # __reduce__ is needed because lambdas passed to self.udp are not picklable.
    def __reduce__(self):
        return (self.__class__, (self.prob_name, self.cut, self.nseg, self.time_encoding, self.inequalities_for_tc))

    def get_name(self):
        return self.name

    def get_extra_info(self):
        return self.extra_info

    def plot(
        self,
        x,
        N=100,
        mark_segments=True,
        mark_mismatch=True,
        orbit_color="lightgray",
        planet_color="gray",
        figsize=(15, 4),
        **kwargs,
    ):
        x_arr = _np.asarray(x)
        fig = _plt.figure(figsize=figsize, layout="constrained")
        ax3d = fig.add_subplot(1, 4, 1, projection="3d")
        ax_xy = fig.add_subplot(1, 4, 2, projection="3d")
        ax_xz = fig.add_subplot(1, 4, 3, projection="3d")
        ax_yz = fig.add_subplot(1, 4, 4, projection="3d")
        axes = (ax3d, ax_xy, ax_xz, ax_yz)
        gym_problem = _pk.trajopt.gym.tops_cr3bp_json[self.prob_name]
        TIME = gym_problem["TIME"]
        period_s_days = gym_problem["period_s"] * TIME * _pk.SEC2DAY
        period_f_days = gym_problem["period_f"] * TIME * _pk.SEC2DAY
        t0 = self.udp.compute_t0(x_arr)
        tof_days = x_arr[10 + 4 * self.udp.nseg] * self.udp.TIME / _pk.DAY2SEC
        for axis in axes:
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.pls, label="pls", units=self.udp.L, color=orbit_color, plot_range=(0, period_s_days), N=N * 5
            )
            _pk.plot.add_planet_orbit(
                ax=axis, pla=self.plf, label="plf", units=self.udp.L, color=orbit_color, plot_range=(0, period_f_days), N=N * 5
            )
            _pk.plot.add_planet(
                ax=axis, when=t0, pla=self.pls, units=self.udp.L, color=planet_color, s=20
            )
            _pk.plot.add_planet(
                ax=axis, when=t0 + tof_days, pla=self.plf, units=self.udp.L, color=planet_color, s=20
            )
            self.udp.plot(
                x=x_arr, ax=axis, N=N, mark_segments=mark_segments, mark_mismatch=mark_mismatch, **kwargs
            )
        for axis in axes:
            xticks = axis.get_xticks()
            yticks = axis.get_yticks()
            zticks = axis.get_zticks()
            if len(xticks) >= 2:
                axis.set_xticks([xticks[0], xticks[-1]])
            if len(yticks) >= 2:
                axis.set_yticks([yticks[0], yticks[-1]])
            if len(zticks) >= 2:
                axis.set_zticks([zticks[0], zticks[-1]])
        ax_xy.set_zticks([])
        ax_xz.set_yticks([])
        ax_yz.set_xticks([])
        ax3d.view_init(elev=20, azim=270)
        ax3d.set_aspect("equal")
        ax_xy.view_init(elev=90, azim=-90)
        ax_xz.view_init(elev=0, azim=-90)
        ax_yz.view_init(elev=0, azim=180)
        ax3d.set_title("3D")
        ax_xy.set_title("xy - view")
        ax_xz.set_title("xz - view")
        ax_yz.set_title("yz - view")
        return axes

    def __getattr__(self, name):
        udp = object.__getattribute__(self, "udp")
        return getattr(udp, name)
