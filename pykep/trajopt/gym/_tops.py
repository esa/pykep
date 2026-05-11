import pykep as _pk
import numpy as _np
from scipy import sparse

class tops_twobody:
    """
    Two-body low-thrust benchmark from the TOPS
    (Trajectory Optimisation Problems in Space) database.

    The TOPS benchmark problems, from Trajectory Optimisation Problems in
    Space, are low-thrust trajectory optimization
    problems transcribed into nonlinear programs by a direct method. In this
    instance the spacecraft dynamics are modeled in Cartesian coordinates under
    two-body motion. The transcription uses a forward-backward shooting scheme,
    which introduces highly nonlinear mismatch constraints, together with
    throttle constraints on the segment controls.

    The control is represented with a zero-order hold over each segment, so the
    resulting NLP size and difficulty depend on the number of segments and on
    the selected time encoding. Different choices of ``nseg`` and
    ``time_encoding`` therefore provide a convenient way to regulate the
    dimensionality and complexity of the benchmark. The formulation is built on
    :class:`~pykep.trajopt.zoh_point2point`, so the departure and arrival
    states are fixed and moving-end effects are not accounted for.
    """

    def __init__(self, prob_name, cut=0.5, nseg=10, time_encoding="uniform"):
        """Construct a two-body TOPS instance from a predefined gym case."""
        self.prob_name = prob_name
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
        states = [it / self.L for it in gym_problem["state_s"][:3]] + [it / self.V for it in gym_problem["state_s"][3:6]]
        statef = [it / self.L for it in gym_problem["state_f"][:3]] + [it / self.V for it in gym_problem["state_f"][3:6]]
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
            inequalities_for_tc=True,
            max_steps=1000,
        )

    def get_name(self):
        """Return the problem name."""
        return self.name
    
    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info
    
    def known_hessians(self, x):
        """
        Return known constraint Hessians as ``(k, hess_func)`` pairs.

        Index contract:
        - ``k`` is the index in the full fitness ordering
          ``[obj, eq_0, ..., eq_{nec-1}, ineq_0, ..., ineq_{nic-1}]``.
        - The API is agnostic to equality/inequality type at this level.

        Sparse contract:
                - ``hess_func(x)`` can return any sparse structure agreed with the
                    consumer (solver-side adapter).
                - For this UDP we return a SciPy ``csc_matrix``.

        In this two-body formulation the known Hessians are those of throttle
        constraints, each depending only on ``(ix, iy, iz)`` of one segment and
        equal to ``2 * I_3`` on that block.
        """
        first_throttle_idx = 1 + self.udp.get_nec()
        known_idxs = [first_throttle_idx + i for i in range(self.udp.nseg - 1)]
        throttle_diag = _np.array([2.0, 2.0, 2.0], dtype=float)
        retval = []
        for i in known_idxs:
            def hess_func(x, idx=i):
                # Segment control is (T, ix, iy, iz); Hessian is nonzero only on ix,iy,iz.
                control_start_idx = 1 + (idx - first_throttle_idx) * 4
                rows = _np.array(
                    [control_start_idx + 1, control_start_idx + 2, control_start_idx + 3],
                    dtype=int,
                )
                cols = rows.copy()
                return sparse.csc_matrix(
                    (throttle_diag.copy(), (rows, cols)), shape=(len(x), len(x))
                )
            retval.append((i, hess_func))
        return retval


    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, 'udp')
        return getattr(udp, name)


class tops_mee:
    """
    Two-body low-thrust benchmark in modified equinoctial elements from the
    TOPS (Trajectory Optimisation Problems in Space) database.

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

    def __init__(self, prob_name, cut=0.5, nseg=10, time_encoding="uniform"):
        """Construct a two-body MEE TOPS instance from a predefined gym case."""
        self.prob_name = prob_name
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
        states = [it / self.L for it in gym_problem["state_s"][:1]] + [it for it in gym_problem["state_s"][1:6]]
        statef = [it / self.L for it in gym_problem["state_f"][:1]] + [it for it in gym_problem["state_f"][1:6]]
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
            inequalities_for_tc=True,
            max_steps=1000,
        )

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info
    
    def known_hessians(self, x):
        """
        Return known constraint Hessians as ``(k, hess_func)`` pairs.

        Index contract:
        - ``k`` is the index in the full fitness ordering
          ``[obj, eq_0, ..., eq_{nec-1}, ineq_0, ..., ineq_{nic-1}]``.
        - The API is agnostic to equality/inequality type at this level.

        Sparse contract:
                - ``hess_func(x)`` can return any sparse structure agreed with the
                    consumer (solver-side adapter).
                - For this UDP we return a SciPy ``csc_matrix``.

        In this two-body formulation the known Hessians are those of throttle
        constraints, each depending only on ``(ix, iy, iz)`` of one segment and
        equal to ``2 * I_3`` on that block.
        """
        first_throttle_idx = 1 + self.udp.get_nec()
        known_idxs = [first_throttle_idx + i for i in range(self.udp.nseg - 1)]
        throttle_diag = _np.array([2.0, 2.0, 2.0], dtype=float)
        retval = []
        for i in known_idxs:
            def hess_func(x, idx=i):
                # Segment control is (T, ix, iy, iz); Hessian is nonzero only on ix,iy,iz.
                control_start_idx = 1 + (idx - first_throttle_idx) * 4
                rows = _np.array(
                    [control_start_idx + 1, control_start_idx + 2, control_start_idx + 3],
                    dtype=int,
                )
                cols = rows.copy()
                return sparse.csc_matrix(
                    (throttle_diag.copy(), (rows, cols)), shape=(len(x), len(x))
                )
            retval.append((i, hess_func))
        return retval

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, 'udp')
        return getattr(udp, name)


class tops_ss:
    """
    Solar-sailing benchmark from the TOPS
    (Trajectory Optimisation Problems in Space) database.

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
        """Construct a solar-sailing TOPS instance from a predefined gym case."""
        self.prob_name = prob_name
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
        self.L = gym_problem["L"]
        # The rest is induced
        self.TIME = _np.sqrt(self.L**3 / self.MU)
        self.V = self.L / self.TIME
        self.ACC = self.V / self.TIME

        # Scaling the problem to these units (mu will be 1. -> as expected from pykep ta_ss)
        states = [it / self.L for it in gym_problem["state_s"][:3]] + [it / self.V for it in gym_problem["state_s"][3:6]]
        statef = [it / self.L for it in gym_problem["state_f"][:3]] + [it / self.V for it in gym_problem["state_f"][3:6]]
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

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, 'udp')
        return getattr(udp, name)


class tops_cr3bp:
    """
    CR3BP low-thrust benchmark from the TOPS
    (Trajectory Optimisation Problems in Space) database.

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

    def __init__(self, prob_name, cut=0.5, nseg=60, time_encoding="uniform"):
        """Construct a CR3BP TOPS instance from a predefined gym case."""
        self.prob_name = prob_name
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
        self.L = 1.
        self.V = 1.
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
            inequalities_for_tc=False,
            max_steps=1000,
        )

    def get_name(self):
        """Return the problem name."""
        return self.name

    def get_extra_info(self):
        """Return additional information for the selected gym case."""
        return self.extra_info
    
    def known_hessians(self, x):
        """
        Return known constraint Hessians as ``(k, hess_func)`` pairs.

        Index contract:
        - ``k`` is the index in the full fitness ordering
          ``[obj, eq_0, ..., eq_{nec-1}, ineq_0, ..., ineq_{nic-1}]``.
        - The API is agnostic to equality/inequality type at this level.

        Sparse contract:
                - ``hess_func(x)`` can return any sparse structure agreed with the
                    consumer (solver-side adapter).
                - For this UDP we return a SciPy ``csc_matrix``.

        In this two-body formulation the known Hessians are those of throttle
        constraints, each depending only on ``(ix, iy, iz)`` of one segment and
        equal to ``2 * I_3`` on that block.
        """
        first_throttle_idx = 1 + self.udp.get_nec()
        known_idxs = [first_throttle_idx + i for i in range(self.udp.nseg - 1)]
        throttle_diag = _np.array([2.0, 2.0, 2.0], dtype=float)
        retval = []
        for i in known_idxs:
            def hess_func(x, idx=i):
                # Segment control is (T, ix, iy, iz); Hessian is nonzero only on ix,iy,iz.
                control_start_idx = 1 + (idx - first_throttle_idx) * 4
                rows = _np.array(
                    [control_start_idx + 1, control_start_idx + 2, control_start_idx + 3],
                    dtype=int,
                )
                cols = rows.copy()
                return sparse.csc_matrix(
                    (throttle_diag.copy(), (rows, cols)), shape=(len(x), len(x))
                )
            retval.append((i, hess_func))
        return retval

    def __getattr__(self, name):
        """Forward any undefined attribute/method calls to self.udp."""
        udp = object.__getattribute__(self, 'udp')
        return getattr(udp, name)
