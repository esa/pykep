import numpy as _np
import pykep as _pk
from . import _zoh_udp_plot_utils as _zoh_plot


class zoh_ms_ss_pl2pl:
    """Represents a solar-sail planet-to-planet transfer using multiple shooting.

    The decision vector preserves the :class:`~pykep.trajopt.zoh_ss_pl2pl`
    prefix and appends the six-component interior shooting states::

        x = [t0_frac, vinf_dep, idep_x, idep_y, idep_z,
             vinf_arr, iarr_x, iarr_y, iarr_z] + controls + [tof] (+ [weights]) + nodes

    Each segment control is ``[alpha, beta]``: sail cone and clock angles in
    the RTN frame. The objective is time of flight; segment defects and the
    boundary excess-velocity direction constraints complete the NLP.
    """

    def __init__(
        self,
        pls=None,
        plf=None,
        t0_bounds=None,
        tof_bounds=None,
        vinf_dep_bounds=None,
        vinf_arr_bounds=None,
        nseg=10,
        cut=0.6,
        tas=None,
        inequalities_for_tc=False,
        time_encoding="uniform",
        w_bounds_softmax=None,
        max_steps=None,
        L=_pk.AU,
        V=_pk.EARTH_VELOCITY,
        state_bounds=None,
    ):
        """Initialize a planet-to-planet solar-sail multiple-shooting problem.

        Args:
            pls (:class:`~pykep.planet`, optional): Departure ephemeris provider.
            plf (:class:`~pykep.planet`, optional): Arrival ephemeris provider.
            t0_bounds (:class:`list`, optional): Departure epoch bounds in MJD2000 days.
                Default is ``[6700, 6800]``.
            tof_bounds (:class:`list`, optional): Time-of-flight bounds in integrator units.
                Default is ``[3.4, 8.6]``.
            vinf_dep_bounds (:class:`list`, optional): Departure excess-speed bounds in
                integrator velocity units. Default is ``[0, 0.2]``.
            vinf_arr_bounds (:class:`list`, optional): Arrival excess-speed bounds in
                integrator velocity units. Default is ``[0, 0.2]``.
            nseg (:class:`int`, optional): Number of shooting segments. Default is 10.
            cut (:class:`float`, optional): Forward/backward split passed to
                :class:`~pykep.leg.zoh_ms`. Default is 0.6.
            tas (:class:`tuple`, optional): Nominal and variational solar-sail integrators.
                Default uses solar-sail dynamics without gradients.
            inequalities_for_tc (:class:`bool`, optional): Form direction-norm constraints
                as inequalities instead of equalities. Default is False.
            time_encoding (:class:`str`, optional): ``'uniform'`` or ``'softmax'`` segment
                durations. Default is ``'uniform'``.
            w_bounds_softmax (:class:`list`, optional): Bounds for softmax weights.
                Default is ``[-1, 1]``.
            max_steps (:class:`int`, optional): Maximum integrator steps per segment.
            L (:class:`float`, optional): Ephemeris length scale. Default is :data:`~pykep.AU`.
            V (:class:`float`, optional): Ephemeris velocity scale. Default is
                :data:`~pykep.EARTH_VELOCITY`.
            state_bounds (:class:`tuple`, optional): Lower and upper bounds for each
                six-component interior state. Default is ``([-100] * 6, [100] * 6)``.
        """
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")
        if time_encoding not in ("uniform", "softmax"):
            raise NotImplementedError("Only uniform and softmax time encodings are supported")

        self.pls = _pk.planet(_pk.udpla.jpl_lp(body="EARTH")) if pls is None else pls
        self.plf = _pk.planet(_pk.udpla.jpl_lp(body="MARS")) if plf is None else plf
        self.t0_bounds = [6700.0, 6800.0] if t0_bounds is None else t0_bounds
        self.tof_bounds = [3.4, 8.6] if tof_bounds is None else tof_bounds
        self.vinf_dep_bounds = [0.0, 0.2] if vinf_dep_bounds is None else vinf_dep_bounds
        self.vinf_arr_bounds = [0.0, 0.2] if vinf_arr_bounds is None else vinf_arr_bounds
        self.nseg = nseg
        self.inequalities_for_tc = inequalities_for_tc
        self.time_encoding = time_encoding
        self.w_bounds_softmax = [-1.0, 1.0] if w_bounds_softmax is None else w_bounds_softmax
        self.max_steps = max_steps
        self.L = L
        self.V = V
        self.TIME = self.L / self.V
        self.ACC = self.V / self.TIME
        bounds = _np.asarray(
            ([-100.0] * 6, [100.0] * 6) if state_bounds is None else state_bounds,
            dtype=float,
        )
        if (
            bounds.shape != (2, 6)
            or not _np.all(_np.isfinite(bounds))
            or _np.any(bounds[0] > bounds[1])
        ):
            raise ValueError("state_bounds must contain finite, ordered six-component bounds")
        self.state_bounds = (bounds[0].tolist(), bounds[1].tolist())
        self._node_lb, self._node_ub = self.state_bounds
        self._single_shooting_size = 10 + 2 * nseg + (nseg if time_encoding == "softmax" else 0)
        self._dtgrid_dtof_uniform = _np.linspace(0.0, 1.0, nseg + 1)
        self._softmax_cumsum_matrix = _np.vstack(
            [_np.zeros(nseg), _np.tril(_np.ones((nseg, nseg)))]
        )
        if tas is None:
            tas = (_pk.ta.get_zoh_ss(1e-10), None)
        self.with_gradient = tas[1] is not None

        controls = _np.random.uniform(-_np.pi / 2, _np.pi / 2, 2 * nseg)
        t0_mid = sum(self.t0_bounds) / 2.0
        tof_mid = sum(self.tof_bounds) / 2.0
        rs, vs = self.pls.eph(t0_mid)
        rf, vf = self.plf.eph(t0_mid + self.TIME * tof_mid / _pk.DAY2SEC)
        state0 = _np.r_[_np.asarray(rs) / self.L, _np.asarray(vs) / self.V]
        statef = _np.r_[_np.asarray(rf) / self.L, _np.asarray(vf) / self.V]
        node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
        self.leg = _pk.leg.zoh_ms(
            state0.tolist() + node * (nseg - 1) + statef.tolist(),
            controls.tolist(),
            _np.linspace(0.0, tof_mid, nseg + 1).tolist(),
            cut,
            tas,
            max_steps=max_steps,
            dim_dynamics=6,
            dim_controls=2,
        )
        self._gradient_pattern = self._make_gradient_sparsity()
        self._gradient_index = {pair: index for index, pair in enumerate(self._gradient_pattern)}

    def _expected_nx(self):
        return self._single_shooting_size + 6 * (self.nseg - 1)

    def compute_t0(self, x):
        """Convert the fractional departure-time variable to MJD2000 days."""
        return self.t0_bounds[0] + x[0] * self.TIME * _pk.SEC2DAY

    def _set_leg_from_x(self, x):
        if len(x) != self._expected_nx():
            raise ValueError(
                f"Invalid decision vector length: got {len(x)}, expected {self._expected_nx()}"
            )
        t0 = self.compute_t0(x)
        vinf_dep, idep = x[1], _np.asarray(x[2:5])
        vinf_arr, iarr = x[5], _np.asarray(x[6:9])
        controls = list(x[9 : 9 + 2 * self.nseg])
        self.leg.controls = controls
        tof_col = 9 + 2 * self.nseg
        tof = x[tof_col]
        tf = t0 + self.TIME * tof / _pk.DAY2SEC
        rs, vs = self.pls.eph(t0)
        rf, vf = self.plf.eph(tf)
        state0 = _np.r_[
            _np.asarray(rs) / self.L,
            _np.asarray(vs) / self.V + vinf_dep * idep,
        ]
        statef = _np.r_[
            _np.asarray(rf) / self.L,
            _np.asarray(vf) / self.V + vinf_arr * iarr,
        ]
        interior = list(x[self._single_shooting_size :])
        self.leg.states = state0.tolist() + interior + statef.tolist()
        if self.time_encoding == "uniform":
            tgrid = tof * self._dtgrid_dtof_uniform
        else:
            weights = x[tof_col + 1 : self._single_shooting_size]
            fractions, _ = _pk.compute_softmax_and_jacobian(weights)
            tgrid = _np.r_[0.0, _np.cumsum(tof * fractions)]
        self.leg.tgrid = tgrid.tolist()
        return (
            rs,
            vs,
            rf,
            vf,
            _np.asarray(rs) / self.L,
            _np.asarray(vs) / self.V,
            _np.asarray(rf) / self.L,
            _np.asarray(vf) / self.V,
        )

    def get_bounds(self):
        """Return bounds in decision-vector order, including interior states."""
        lb = (
            [0.0, self.vinf_dep_bounds[0], -1.0, -1.0, -1.0,
             self.vinf_arr_bounds[0], -1.0, -1.0, -1.0]
            + ([-_np.pi / 2 + 1e-3, 0.0] * self.nseg)
            + [self.tof_bounds[0]]
        )
        ub = (
            [(self.t0_bounds[1] - self.t0_bounds[0]) * _pk.DAY2SEC / self.TIME,
             self.vinf_dep_bounds[1], 1.0, 1.0, 1.0,
             self.vinf_arr_bounds[1], 1.0, 1.0, 1.0]
            + ([_np.pi / 2 - 1e-3, 2.0 * _np.pi] * self.nseg)
            + [self.tof_bounds[1]]
        )
        if self.time_encoding == "softmax":
            lb += [self.w_bounds_softmax[0]] * self.nseg
            ub += [self.w_bounds_softmax[1]] * self.nseg
        return lb + self._node_lb * (self.nseg - 1), ub + self._node_ub * (self.nseg - 1)

    def guess_internal_states(self, x, ballistic=False):
        """Return a complete decision vector with propagated interior states.

        Args:
            x (:class:`list`): A single- or multiple-shooting decision vector.
            ballistic (:class:`bool`, optional): Propagate with zero sail controls.
                Default is False.

        Raises:
            ValueError: If the vector is invalid or a propagated node is outside its bounds.

        Returns:
            :class:`list`: The input prefix followed by flattened interior states.
        """
        x = list(x)
        if len(x) == self._single_shooting_size:
            node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
            x += node * (self.nseg - 1)
        elif len(x) != self._expected_nx():
            raise ValueError("Invalid decision vector length")
        self._set_leg_from_x(x)
        self.leg.set_initial_guess(ballistic=ballistic)
        nodes = _np.asarray(self.leg.states[6:-6]).reshape(self.nseg - 1, 6)
        if (
            not _np.all(_np.isfinite(nodes))
            or _np.any(nodes < self._node_lb)
            or _np.any(nodes > self._node_ub)
        ):
            raise ValueError("Propagated interior nodes are outside their bounds or non-finite")
        return x[: self._single_shooting_size] + nodes.ravel().tolist()

    def fitness(self, x):
        """Evaluate time of flight, segment defects, and boundary directions."""
        self._set_leg_from_x(x)
        idep, iarr = _np.asarray(x[2:5]), _np.asarray(x[6:9])
        directions = [idep @ idep - 1.0, iarr @ iarr - 1.0]
        return [x[9 + 2 * self.nseg]] + self.leg.compute_defects() + directions

    def _make_gradient_sparsity(self):
        pattern = {(0, 9 + 2 * self.nseg)}
        tof_col = 9 + 2 * self.nseg
        for segment in range(self.nseg):
            for component in range(6):
                row = 1 + 6 * segment + component
                for node in (segment, segment + 1):
                    if 0 < node < self.nseg:
                        for state_component in range(6):
                            pattern.add(
                                (row, self._single_shooting_size + 6 * (node - 1) + state_component)
                            )
                    elif node == 0:
                        for column in range(5):
                            pattern.add((row, column))
                    else:
                        for column in (0, 5, 6, 7, 8, tof_col):
                            pattern.add((row, column))
                for control in range(2):
                    pattern.add((row, 9 + 2 * segment + control))
                pattern.add((row, tof_col))
                if self.time_encoding == "softmax":
                    for weight in range(self.nseg):
                        pattern.add((row, tof_col + 1 + weight))
        for component in range(3):
            pattern.add((1 + 6 * self.nseg, 2 + component))
            pattern.add((2 + 6 * self.nseg, 6 + component))
        return sorted(pattern)

    def gradient_sparsity(self):
        """Return the sparse fitness-Jacobian pattern in pygmo ordering."""
        return list(self._gradient_pattern)

    def gradient(self, x):
        """Return sparse derivatives of the objective and constraints.

        Raises:
            RuntimeError: If no variational integrator was supplied.
            NotImplementedError: If either ephemeris lacks an acceleration method.
        """
        if not self.with_gradient:
            raise RuntimeError("Gradient computation requires a variational integrator")
        if not hasattr(self.pls, "acc") or not hasattr(self.plf, "acc"):
            raise NotImplementedError("Gradient computation requires ephemeris acceleration methods")
        rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd = self._set_leg_from_x(x)
        nx = len(x)
        tof_col = 9 + 2 * self.nseg
        tof = x[tof_col]
        state_pattern, control_pattern, time_pattern = self.leg.defects_grad_sparsity()
        state_gradient, control_gradient, time_gradient = self.leg.compute_defects_grad()
        state_jacobian = _np.zeros((6 * self.nseg, 6 * (self.nseg + 1)))
        control_jacobian = _np.zeros((6 * self.nseg, 2 * self.nseg))
        time_jacobian = _np.zeros((6 * self.nseg, self.nseg + 1))
        for (row, column), value in zip(state_pattern, state_gradient):
            state_jacobian[int(row), int(column)] += value
        for (row, column), value in zip(control_pattern, control_gradient):
            control_jacobian[int(row), int(column)] += value
        for (row, column), value in zip(time_pattern, time_gradient):
            time_jacobian[int(row), int(column)] += value

        as_nd = _np.asarray(self.pls.acc(self.compute_t0(x))) / self.ACC
        af_nd = _np.asarray(self.plf.acc(self.compute_t0(x) + self.TIME * tof / _pk.DAY2SEC)) / self.ACC
        boundary0 = _np.zeros((6, nx))
        boundary1 = _np.zeros((6, nx))
        boundary0[:3, 0] = vs_nd * _pk.DAY2SEC / self.TIME
        boundary0[3:6, 0] = as_nd * _pk.DAY2SEC / self.TIME
        boundary1[:3, 0] = vf_nd * _pk.DAY2SEC / self.TIME
        boundary1[3:6, 0] = af_nd * _pk.DAY2SEC / self.TIME
        boundary1[:3, tof_col] = vf_nd
        boundary1[3:6, tof_col] = af_nd
        boundary0[3:6, 1] = x[2:5]
        boundary0[3:6, 2:5] = x[1] * _np.eye(3)
        boundary1[3:6, 5] = x[6:9]
        boundary1[3:6, 6:9] = x[5] * _np.eye(3)

        gradient = _np.zeros((len(self.fitness(x)), nx))
        gradient[0, tof_col] = 1.0
        for node in range(self.nseg + 1):
            node_jacobian = state_jacobian[:, 6 * node : 6 * (node + 1)]
            if node == 0:
                gradient[1 : 1 + 6 * self.nseg] += node_jacobian @ boundary0
            elif node == self.nseg:
                gradient[1 : 1 + 6 * self.nseg] += node_jacobian @ boundary1
            else:
                start = self._single_shooting_size + 6 * (node - 1)
                gradient[1 : 1 + 6 * self.nseg, start : start + 6] += node_jacobian

        gradient[1 : 1 + 6 * self.nseg, 9 : 9 + 2 * self.nseg] += control_jacobian
        if self.time_encoding == "uniform":
            dtgrid_dtof = self._dtgrid_dtof_uniform
        else:
            fractions, jacobian = _pk.compute_softmax_and_jacobian(
                x[tof_col + 1 : self._single_shooting_size]
            )
            dtgrid_dtof = _np.r_[0.0, _np.cumsum(fractions)]
            dtgrid_dw = tof * (self._softmax_cumsum_matrix @ jacobian)
        gradient[1 : 1 + 6 * self.nseg, tof_col] += time_jacobian @ dtgrid_dtof
        if self.time_encoding == "softmax":
            gradient[1 : 1 + 6 * self.nseg, tof_col + 1 : self._single_shooting_size] += (
                time_jacobian @ dtgrid_dw
            )

        gradient[1 + 6 * self.nseg, 2:5] = 2.0 * _np.asarray(x[2:5])
        gradient[2 + 6 * self.nseg, 6:9] = 2.0 * _np.asarray(x[6:9])
        gradient[:, 0] *= self.TIME * _pk.SEC2DAY
        return _np.asarray(
            [gradient[row, column] for row, column in self._gradient_pattern]
        )

    def has_gradient(self):
        """Return whether variational dynamics are available."""
        return self.with_gradient

    def get_nec(self):
        """Return the number of equality constraints."""
        return 6 * self.nseg + (0 if self.inequalities_for_tc else 2)

    def get_nic(self):
        """Return the number of inequality constraints."""
        return 2 if self.inequalities_for_tc else 0

    def pretty(self, x):
        """Print transfer epochs and the multiple-shooting leg."""
        self._set_leg_from_x(x)
        tof_days = x[9 + 2 * self.nseg] * self.TIME / _pk.DAY2SEC
        t0 = self.compute_t0(x)
        print("\nSolar Sail ZOH multiple-shooting transfer")
        print(f"Departure: {self.pls.name}\nArrival: {self.plf.name}")
        print(f"Launch epoch: {t0:.5f} MJD2000")
        print(f"Arrival epoch: {t0 + tof_days:.5f} MJD2000")
        print(self.leg)

    def plot(
        self,
        x,
        ax=None,
        N=30,
        mark_segments=True,
        mark_mismatch=True,
        plot_sail=True,
        sail_size=0.05,
        **kwargs,
    ):
        """Plot the propagated segments and their sail orientations.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.
            ax (:class:`matplotlib.axes.Axes`, optional): Three-dimensional axes.
            N (:class:`int`, optional): Samples per segment. Default is 30.
            mark_segments (:class:`bool`, optional): Mark segment starts. Default is True.
            mark_mismatch (:class:`bool`, optional): Mark propagated segment ends. Default is True.
            plot_sail (:class:`bool`, optional): Draw the sail and its normal. Default is True.
            sail_size (:class:`float`, optional): Half-side length of the sail patch.
            **kwargs: Additional arguments for segment markers.

        Returns:
            :class:`matplotlib.axes.Axes`: The trajectory axes.
        """
        self._set_leg_from_x(x)
        if ax is None:
            ax = _pk.plot.make_3Daxis()
        forward, backward, _ = self.leg.get_state_info(N=N)
        angles = _np.asarray(x[9 : 9 + 2 * self.nseg]).reshape(self.nseg, 2)

        groups = (
            (forward, range(self.leg.nseg_fwd)),
            (backward, range(self.nseg - 1, self.leg.nseg_fwd - 1, -1)),
        )
        for segments, indices in groups:
            for segment, index in zip(segments, indices):
                if not len(segment):
                    continue
                cart = _np.asarray(segment)
                ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], c="tab:orange")
                if mark_segments:
                    ax.scatter(*cart[0, :3], **kwargs)
                if mark_mismatch:
                    ax.scatter(*cart[-1, :3], marker="^", **kwargs)
                alpha, beta = angles[index]
                if plot_sail and abs(abs(alpha) - _np.pi / 2) > 1e-3:
                    _zoh_plot.draw_solar_sail(
                        ax, segment[len(segment) // 2], alpha, beta, sail_size
                    )
        _zoh_plot.set_3d_axes_equal(ax)
        return ax
