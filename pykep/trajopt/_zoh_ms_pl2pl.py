import numpy as _np
import pykep as _pk
from . import _zoh_udp_plot_utils as _zoh_plot


class zoh_ms_pl2pl:
    """Represents the optimal low-thrust transfer between two planets using multiple shooting.

    This problem uses :class:`~pykep.leg.zoh_ms` to connect planet ephemerides with a low-thrust
    trajectory. It optimizes departure epoch, final mass, boundary excess velocities, controls,
    time of flight, and the interior states at each shooting node.

    The decision vector is::

        x = [t0_frac, mf, vinf_dep, idep_x, idep_y, idep_z,
             vinf_arr, iarr_x, iarr_y, iarr_z] + controls + [tof] (+ [weights]) + nodes

    where:

    - ``t0_frac`` is the fractional departure epoch in integrator time units.
    - ``mf`` is the final spacecraft mass.
    - ``vinf_dep`` and ``vinf_arr`` are the magnitudes of the boundary excess velocities;
      ``idep`` and ``iarr`` are their direction vectors.
    - Each control is ``[throttle, ix, iy, iz]``; throttle is scaled by ``max_thrust``.
    - ``weights`` are present only for softmax time encoding.
    - ``nodes`` are a flat tail of consecutive seven-component interior states, including mass,
        in chronological order.

    The objective maximizes final mass. Segment defects are equalities; throttle and boundary
    direction constraints may be formulated as equalities or inequalities.
    """
        
    def __init__(
        self,
        pls=None,
        plf=None,
        ms=1.0,
        max_thrust=0.22,
        t0_bounds=None,
        tof_bounds=None,
        mf_bounds=None,
        vinf_dep_bounds=None,
        vinf_arr_bounds=None,
        nseg=10,
        cut=0.6,
        tas=None,
        cart2state=None,
        state2cart=None,
        inequalities_for_tc=False,
        time_encoding="uniform",
        w_bounds_softmax=None,
        max_steps=None,
        nrevs=None,
        L=_pk.AU,
        V=_pk.EARTH_VELOCITY,
        state_bounds=None,
        interior_mass_bounds=None,
    ):
        """Initialize a planet-to-planet multiple-shooting problem.

        Args:
            pls (:class:`~pykep.planet`, optional): Departure ephemeris provider.
            plf (:class:`~pykep.planet`, optional): Arrival ephemeris provider.
            ms (:class:`float`, optional): Initial spacecraft mass in dynamics units. Default is 1.
            max_thrust (:class:`float`, optional): Maximum thrust in dynamics units. Default is 0.22.
            t0_bounds (:class:`list`, optional): Departure epoch bounds in MJD2000 days. Default is ``[6700, 6800]``.
            tof_bounds (:class:`list`, optional): Time-of-flight bounds in integrator units. Default is ``[3.4, 8.6]``.
            mf_bounds (:class:`list`, optional): Final mass bounds. Default is ``[0.2, 1]``.
            vinf_dep_bounds (:class:`list`, optional): Departure excess-speed bounds in integrator units. Default is ``[0, 0.2]``.
            vinf_arr_bounds (:class:`list`, optional): Arrival excess-speed bounds in integrator units. Default is ``[0, 0.2]``.
            nseg (:class:`int`, optional): Number of shooting segments. Default is 10.
            cut (:class:`float`, optional): Forward/backward segment split passed to :class:`~pykep.leg.zoh_ms`. Default is 0.6.
            tas (:class:`tuple`, optional): Nominal and variational Taylor integrators. Default uses Kepler dynamics without gradients.
            cart2state (:class:`list`, optional): State transform and its Jacobian, as in :class:`~pykep.trajopt.zoh_pl2pl`.
            state2cart (:class:`callable`, optional): Convert a dynamics state to Cartesian coordinates for plotting.
            inequalities_for_tc (:class:`bool`, optional): Use inequalities for throttle and boundary direction norms. Default is False.
            time_encoding (:class:`str`, optional): ``'uniform'`` or ``'softmax'`` segment durations. Default is ``'uniform'``.
            w_bounds_softmax (:class:`list`, optional): Softmax weight bounds. Default is ``[-1, 1]``.
            max_steps (:class:`int`, optional): Maximum integrator steps per segment. Default is None.
            nrevs (:class:`int`, optional): Additional revolutions for mean-longitude state representations. Default is None.
            L (:class:`float`, optional): Ephemeris length scale. Default is :data:`~pykep.AU`.
            V (:class:`float`, optional): Ephemeris velocity scale. Default is :data:`~pykep.EARTH_VELOCITY`.
            state_bounds (:class:`tuple`, optional): Interior-node state lower and upper bounds. Default is ``([-100] * 6, [100] * 6)``.
            interior_mass_bounds (:class:`tuple`, optional): Interior-node mass bounds. Default is ``(0.1 * ms, 4 * ms)``.
        """
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")
        if time_encoding not in ("uniform", "softmax"):
            raise NotImplementedError("Only uniform and softmax time encodings are supported")
        if cart2state is not None:
            if len(cart2state) != 2 or not callable(cart2state[0]):
                raise ValueError("cart2state must contain a state transform and its Jacobian")
            if tas is not None and tas[1] is not None and not callable(cart2state[1]):
                raise ValueError("cart2state Jacobian is required when using gradients")
        if state2cart is not None and not callable(state2cart):
            raise ValueError("state2cart must be a callable function")

        self.pls = _pk.planet(_pk.udpla.jpl_lp(body="EARTH")) if pls is None else pls
        self.plf = _pk.planet(_pk.udpla.jpl_lp(body="MARS")) if plf is None else plf
        self.initial_mass = ms
        self.max_thrust = max_thrust
        self.t0_bounds = [6700.0, 6800.0] if t0_bounds is None else t0_bounds
        self.tof_bounds = [3.4, 8.6] if tof_bounds is None else tof_bounds
        self.mf_bounds = [0.2, 1.0] if mf_bounds is None else mf_bounds
        self.vinf_dep_bounds = [0.0, 0.2] if vinf_dep_bounds is None else vinf_dep_bounds
        self.vinf_arr_bounds = [0.0, 0.2] if vinf_arr_bounds is None else vinf_arr_bounds
        self.nseg = nseg
        self.with_gradient = tas is not None and tas[1] is not None
        self.cart2state = cart2state
        self.state2cart = state2cart
        self.inequalities_for_tc = inequalities_for_tc
        self.time_encoding = time_encoding
        self.w_bounds_softmax = [-1.0, 1.0] if w_bounds_softmax is None else w_bounds_softmax
        self.max_steps = max_steps
        self.nrevs = nrevs
        self.L = L
        self.V = V
        self.TIME = self.L / self.V
        self.ACC = self.V / self.TIME
        self.state_bounds = ([-100.0] * 6, [100.0] * 6) if state_bounds is None else state_bounds
        state_bounds_array = _np.asarray(self.state_bounds, dtype=float)
        if (
            state_bounds_array.shape != (2, 6)
            or not _np.all(_np.isfinite(state_bounds_array))
            or _np.any(state_bounds_array[0] > state_bounds_array[1])
        ):
            raise ValueError("state_bounds must contain finite, ordered six-component bounds")
        self.state_bounds = (state_bounds_array[0].tolist(), state_bounds_array[1].tolist())
        if interior_mass_bounds is None:
            interior_mass_bounds = (0.1 * ms, 4.0 * ms)
        mass_bounds = _np.asarray(interior_mass_bounds, dtype=float)
        if (
            mass_bounds.shape != (2,)
            or not _np.all(_np.isfinite(mass_bounds))
            or mass_bounds[0] <= 0.0
            or mass_bounds[0] > mass_bounds[1]
        ):
            raise ValueError("interior_mass_bounds must contain positive, finite, ordered bounds")
        self.interior_mass_bounds = mass_bounds.tolist()
        self._node_lb = self.state_bounds[0] + [mass_bounds[0]]
        self._node_ub = self.state_bounds[1] + [mass_bounds[1]]
        # Preserve the single-shooting prefix; multiple shooting appends seven-value interior nodes.
        self._single_shooting_size = 11 + 4 * nseg + (nseg if time_encoding == "softmax" else 0)
        self._dtgrid_dtof_uniform = _np.linspace(0.0, 1.0, nseg + 1)
        self._softmax_cumsum_matrix = _np.vstack(
            [_np.zeros(nseg), _np.tril(_np.ones((nseg, nseg)))]
        )
        if tas is None:
            tas = (_pk.ta.get_zoh_kep(1e-10), None)

        controls = _np.random.uniform(-1.0, 1.0, 4 * nseg)
        controls[0::4] = _np.abs(controls[0::4]) * max_thrust
        t0_mid = sum(self.t0_bounds) / 2.0
        tof_mid = sum(self.tof_bounds) / 2.0
        rs, vs = self.pls.eph(t0_mid)
        rf, vf = self.plf.eph(t0_mid + self.TIME * tof_mid / _pk.DAY2SEC)
        state0 = self._planet_state(rs, vs, _np.zeros(3), self.initial_mass)
        statef = self._planet_state(rf, vf, _np.zeros(3), self.initial_mass)
        node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
        # zoh_ms stores the departure node, chronological interior nodes, and arrival node in one flat vector.
        states = state0 + node * (nseg - 1) + statef
        tgrid = _np.linspace(0.0, tof_mid, nseg + 1).tolist()
        self.leg = _pk.leg.zoh_ms(
            states, controls.tolist(), tgrid, cut, tas,
            max_steps=max_steps, dim_dynamics=7, dim_controls=4,
        )
        self._gradient_pattern = self._make_gradient_sparsity()
        self._gradient_index = {pair: index for index, pair in enumerate(self._gradient_pattern)}

    def _planet_state(self, r, v, dv, mass):
        state = _np.r_[_np.asarray(r) / self.L, _np.asarray(v) / self.V + dv]
        if self.cart2state is not None:
            state = _np.asarray(self.cart2state[0](state), dtype=float)
        return state.tolist() + [mass]

    def _expected_nx(self):
        return self._single_shooting_size + 7 * (self.nseg - 1)

    def compute_t0(self, x):
        """Convert the fractional departure-time variable to MJD2000 days.

        Args:
            x (:class:`list`): Decision vector whose first element is the fractional epoch.

        Returns:
            :class:`float`: Departure epoch in MJD2000 days.
        """
        return self.t0_bounds[0] + x[0] * self.TIME * _pk.SEC2DAY

    def _set_leg_from_x(self, x):
        if len(x) != self._expected_nx():
            raise ValueError(
                f"Invalid decision vector length: got {len(x)}, expected {self._expected_nx()}"
            )
        # Decode the shared planet-to-planet prefix and convert normalized throttle to dynamics units.
        t0 = self.compute_t0(x)
        vinf_dep, idep = x[2], _np.asarray(x[3:6])
        vinf_arr, iarr = x[6], _np.asarray(x[7:10])
        controls = _np.array(x[10 : 10 + 4 * self.nseg], copy=True)
        controls[0::4] *= self.max_thrust
        self.leg.controls = controls.tolist()
        tof_col = 10 + 4 * self.nseg
        tof = x[tof_col]
        tf = t0 + self.TIME * tof / _pk.DAY2SEC
        # Boundary states move with the trial epochs and include the free excess-velocity vectors.
        rs, vs = self.pls.eph(t0)
        rf, vf = self.plf.eph(tf)
        state0 = self._planet_state(rs, vs, vinf_dep * idep, self.initial_mass)
        statef = self._planet_state(rf, vf, vinf_arr * iarr, x[1])
        if self.cart2state is not None and self.nrevs is not None:
            if statef[-2] < state0[-2]:
                statef[-2] += 2.0 * _np.pi
            statef[-2] += self.nrevs * 2.0 * _np.pi
        interior = list(x[self._single_shooting_size :])
        self.leg.states = state0 + interior + statef
        # Build segment times either uniformly or from positive softmax fractions.
        if self.time_encoding == "uniform":
            tgrid = tof * self._dtgrid_dtof_uniform
        else:
            weights = x[tof_col + 1 : self._single_shooting_size]
            fractions, _ = _pk.compute_softmax_and_jacobian(weights)
            tgrid = _np.r_[0.0, _np.cumsum(tof * fractions)]
        self.leg.tgrid = tgrid.tolist()
        return rs, vs, rf, vf, _np.asarray(rs) / self.L, _np.asarray(vs) / self.V, _np.asarray(rf) / self.L, _np.asarray(vf) / self.V

    def get_bounds(self):
        """Return lower and upper bounds in decision-vector order.

        Returns:
            :class:`tuple`: Lower- and upper-bound lists, including the planet-to-planet
            variables, optional softmax weights, and flattened interior nodes.
        """
        lb = (
            [0.0, self.mf_bounds[0], self.vinf_dep_bounds[0], -1.0, -1.0, -1.0,
             self.vinf_arr_bounds[0], -1.0, -1.0, -1.0]
            + [0.0, -1.0, -1.0, -1.0] * self.nseg
            + [self.tof_bounds[0]]
        )
        ub = (
            [(self.t0_bounds[1] - self.t0_bounds[0]) * _pk.DAY2SEC / self.TIME,
             self.mf_bounds[1], self.vinf_dep_bounds[1], 1.0, 1.0, 1.0,
             self.vinf_arr_bounds[1], 1.0, 1.0, 1.0]
            + [1.0, 1.0, 1.0, 1.0] * self.nseg
            + [self.tof_bounds[1]]
        )
        if self.time_encoding == "softmax":
            lb += [self.w_bounds_softmax[0]] * self.nseg
            ub += [self.w_bounds_softmax[1]] * self.nseg
        return lb + self._node_lb * (self.nseg - 1), ub + self._node_ub * (self.nseg - 1)

    def guess_internal_states(self, x, ballistic=False):
        """Return a complete decision vector with propagated interior nodes.

        Args:
            x (:class:`list`): Planet-to-planet single- or multiple-shooting vector.
            ballistic (:class:`bool`, optional): Use zero controls for node propagation. Default is False.

        Raises:
            ValueError: If the vector length is invalid or a guessed node is outside its bounds.

        Returns:
            :class:`list`: Complete decision vector, preserving the input prefix and
            appending the flattened seven-component interior nodes.
        """
        x = list(x)
        if len(x) == self._single_shooting_size:
            node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
            x += node * (self.nseg - 1)
        elif len(x) != self._expected_nx():
            raise ValueError("Invalid decision vector length")
        self._set_leg_from_x(x)
        self.leg.set_initial_guess(ballistic=ballistic)
        nodes = _np.asarray(self.leg.states[7:-7]).reshape(self.nseg - 1, 7)
        if (
            not _np.all(_np.isfinite(nodes))
            or _np.any(nodes < self._node_lb)
            or _np.any(nodes > self._node_ub)
        ):
            raise ValueError("Propagated interior nodes are outside their bounds or non-finite")
        return x[:self._single_shooting_size] + nodes.ravel().tolist()

    def _compute_throttle_constraints(self):
        controls = _np.asarray(self.leg.controls).reshape(self.nseg, 4)
        return (_np.sum(controls[:, 1:] ** 2, axis=1) - 1.0).tolist()

    def fitness(self, x):
        """Evaluate the objective, segment defects, and boundary/control constraints.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.

        Returns:
            :class:`list`: Objective ``-mf``, followed by ``7 * nseg`` segment defects,
            ``nseg`` throttle constraints, and two boundary-direction norm constraints.
        """
        self._set_leg_from_x(x)
        idep, iarr = _np.asarray(x[3:6]), _np.asarray(x[7:10])
        directions = [idep @ idep - 1.0, iarr @ iarr - 1.0]
        return [-x[1]] + self.leg.compute_defects() + self._compute_throttle_constraints() + directions

    def _make_gradient_sparsity(self):
        pattern = {(0, 1)}
        tof_col = 10 + 4 * self.nseg
        for segment in range(self.nseg):
            for component in range(7):
                row = 1 + 7 * segment + component
                for node in (segment, segment + 1):
                    if 0 < node < self.nseg:
                        for state_component in range(7):
                            pattern.add((row, self._single_shooting_size + 7 * (node - 1) + state_component))
                    elif node == 0:
                        for column in (0, 2, 3, 4, 5):
                            pattern.add((row, column))
                    else:
                        for column in (0, 1, 6, 7, 8, 9, tof_col):
                            pattern.add((row, column))
                for control in range(4):
                    pattern.add((row, 10 + 4 * segment + control))
                pattern.add((row, tof_col))
                if self.time_encoding == "softmax":
                    for weight in range(self.nseg):
                        pattern.add((row, tof_col + 1 + weight))
            for direction in range(1, 4):
                pattern.add((1 + 7 * self.nseg + segment, 10 + 4 * segment + direction))
        direction_row = 1 + 8 * self.nseg
        for column in (3, 4, 5):
            pattern.add((direction_row, column))
        for column in (7, 8, 9):
            pattern.add((direction_row + 1, column))
        return sorted(pattern)

    def gradient_sparsity(self):
        """Return the sparse fitness-Jacobian pattern in pygmo ordering.

        Returns:
            :class:`list`: Sorted ``(row, column)`` indices corresponding to
            :meth:`gradient` values.
        """
        return list(self._gradient_pattern)

    def gradient(self, x):
        """Return sparse derivatives of the objective and constraints.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.

        Raises:
            RuntimeError: If no variational integrator was supplied.

        Returns:
            :class:`numpy.ndarray`: Nonzero derivatives ordered as in :meth:`gradient_sparsity`.
        """
        if not self.with_gradient:
            raise RuntimeError("Gradient computation requires a variational integrator")
        rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd = self._set_leg_from_x(x)
        nrows, ncols = len(self.fitness(x)), len(x)
        gradient = _np.zeros((nrows, ncols))
        gradient[0, 1] = -1.0
        # Collect the leg's derivatives by state, control, and grid variable before chaining them.
        state_sp, control_sp, tgrid_sp = self.leg.defects_grad_sparsity()
        state_grad, control_grad, tgrid_grad = self.leg.compute_defects_grad()
        J_states = _np.zeros((7 * self.nseg, 7 * (self.nseg + 1)))
        J_controls = _np.zeros((7 * self.nseg, 4 * self.nseg))
        J_tgrid = _np.zeros((7 * self.nseg, self.nseg + 1))
        for (row, column), value in zip(state_sp, state_grad):
            J_states[int(row), int(column)] += value
        for (row, column), value in zip(control_sp, control_grad):
            J_controls[int(row), int(column)] += value
        for (row, column), value in zip(tgrid_sp, tgrid_grad):
            J_tgrid[int(row), int(column)] += value

        t0 = self.compute_t0(x)
        tof_col = 10 + 4 * self.nseg
        tof = x[tof_col]
        as_nd = _np.asarray(self.pls.acc(t0)) / self.ACC
        af_nd = _np.asarray(self.plf.acc(t0 + self.TIME * tof / _pk.DAY2SEC)) / self.ACC
        # Boundary-node sensitivities include ephemeris motion, final mass, and time of flight.
        B0 = _np.zeros((7, ncols))
        B1 = _np.zeros((7, ncols))
        B0[:3, 0] = vs_nd * _pk.DAY2SEC / self.TIME
        B0[3:6, 0] = as_nd * _pk.DAY2SEC / self.TIME
        B1[:3, 0] = vf_nd * _pk.DAY2SEC / self.TIME
        B1[3:6, 0] = af_nd * _pk.DAY2SEC / self.TIME
        B1[:3, tof_col] = vf_nd
        B1[3:6, tof_col] = af_nd
        B1[6, 1] = 1.0

        dep_velocity = vs_nd + x[2] * _np.asarray(x[3:6])
        arr_velocity = vf_nd + x[6] * _np.asarray(x[7:10])
        if self.cart2state is None:
            J0 = J1 = _np.eye(6)
        else:
            # Transform Cartesian boundary sensitivities into the dynamics state coordinates.
            J0 = _np.asarray(self.cart2state[1](_np.r_[rs_nd, dep_velocity])).reshape(6, 6)
            J1 = _np.asarray(self.cart2state[1](_np.r_[rf_nd, arr_velocity])).reshape(6, 6)
        dep_velocity_columns = J0[:, 3:6]
        arr_velocity_columns = J1[:, 3:6]
        B0[:6, 0] = J0 @ B0[:6, 0]
        B1[:6, 0] = J1 @ B1[:6, 0]
        B1[:6, tof_col] = J1 @ B1[:6, tof_col]
        B0[:6, 2] = dep_velocity_columns @ x[3:6]
        B0[:6, 3:6] = dep_velocity_columns * x[2]
        B1[:6, 6] = arr_velocity_columns @ x[7:10]
        B1[:6, 7:10] = arr_velocity_columns * x[6]

        # Interior nodes map directly to decision variables; endpoint nodes use their boundary chain rules.
        for node in range(self.nseg + 1):
            node_jac = J_states[:, 7 * node : 7 * (node + 1)]
            if node == 0:
                gradient[1 : 1 + 7 * self.nseg] += node_jac @ B0
            elif node == self.nseg:
                gradient[1 : 1 + 7 * self.nseg] += node_jac @ B1
            else:
                columns = slice(self._single_shooting_size + 7 * (node - 1), self._single_shooting_size + 7 * node)
                gradient[1 : 1 + 7 * self.nseg, columns] += node_jac

        for segment in range(self.nseg):
            rows = slice(1 + 7 * segment, 1 + 7 * (segment + 1))
            columns = slice(10 + 4 * segment, 10 + 4 * (segment + 1))
            gradient[rows, columns] += J_controls[7 * segment : 7 * (segment + 1), 4 * segment : 4 * (segment + 1)]
            gradient[rows, 10 + 4 * segment] *= self.max_thrust

        # Chain grid-time derivatives to tof and, when present, the softmax segment weights.
        if self.time_encoding == "uniform":
            dtgrid_dtof = self._dtgrid_dtof_uniform
        else:
            fractions, jacobian = _pk.compute_softmax_and_jacobian(x[tof_col + 1 : self._single_shooting_size])
            dtgrid_dtof = _np.r_[0.0, _np.cumsum(fractions)]
            dtgrid_dw = tof * (self._softmax_cumsum_matrix @ jacobian)
        gradient[1 : 1 + 7 * self.nseg, tof_col] += J_tgrid @ dtgrid_dtof
        if self.time_encoding == "softmax":
            gradient[1 : 1 + 7 * self.nseg, tof_col + 1 : self._single_shooting_size] += J_tgrid @ dtgrid_dw

        # Add throttle and boundary-direction norm derivatives after the segment defects.
        controls = _np.asarray(self.leg.controls).reshape(self.nseg, 4)
        for segment in range(self.nseg):
            row = 1 + 7 * self.nseg + segment
            gradient[row, 10 + 4 * segment + 1 : 10 + 4 * segment + 4] = 2.0 * controls[segment, 1:]
        gradient[1 + 8 * self.nseg, 3:6] = 2.0 * _np.asarray(x[3:6])
        gradient[2 + 8 * self.nseg, 7:10] = 2.0 * _np.asarray(x[7:10])
        gradient[:, 0] *= self.TIME * _pk.SEC2DAY
        return _np.asarray([gradient[row, column] for row, column in self._gradient_pattern])

    def has_gradient(self):
        """Return whether a variational integrator is available for gradients.

        Returns:
            :class:`bool`: True when the problem was constructed with variational dynamics.
        """
        return self.with_gradient

    def get_nec(self):
        """Return the number of equality constraints.

        Returns:
            :class:`int`: Segment defects plus throttle and boundary-direction constraints
            when those constraints are not configured as inequalities.
        """
        return 7 * self.nseg + (0 if self.inequalities_for_tc else self.nseg + 2)

    def get_nic(self):
        """Return the number of inequality constraints.

        Returns:
            :class:`int`: ``nseg + 2`` when throttle and boundary-direction constraints
            are inequalities, otherwise zero.
        """
        return self.nseg + 2 if self.inequalities_for_tc else 0

    def pretty(self, x):
        """Print the transfer epochs, free velocities, and multiple-shooting leg.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.
        """
        self._set_leg_from_x(x)
        t0 = self.compute_t0(x)
        tof_days = x[10 + 4 * self.nseg] * self.TIME / _pk.DAY2SEC
        print("\nLow-thrust ZOH multiple-shooting transfer")
        print(f"Departure: {self.pls.name}\nArrival: {self.plf.name}")
        print(f"Launch epoch: {t0:.5f} MJD2000")
        print(f"Arrival epoch: {t0 + tof_days:.5f} MJD2000")
        print(self.leg)

    def plot(self, x, ax=None, N=30, mark_segments=True, mark_mismatch=True,
             orbit_color=None, tof=None, **kwargs):
        """Plot the propagated multiple-shooting segments and optional boundary orbits.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.
            ax (:class:`matplotlib.axes.Axes`, optional): Three-dimensional axes. A new
                axes is created when omitted.
            N (:class:`int`, optional): Samples per segment. Default is 30.
            mark_segments (:class:`bool`, optional): Mark segment starts. Default is True.
            mark_mismatch (:class:`bool`, optional): Mark propagated segment endpoints.
                Default is True.
            orbit_color (:class:`str`, optional): Plot the boundary ballistic orbits in this color.
            tof (:class:`float`, optional): Duration for boundary-orbit plots; defaults to
                the decision-vector time of flight.
            **kwargs: Additional keyword arguments for segment markers.

        Returns:
            :class:`matplotlib.axes.Axes`: Axes containing the trajectory plot.
        """
        self._set_leg_from_x(x)
        if ax is None:
            ax = _pk.plot.make_3Daxis()

        def to_cartesian(state):
            if self.state2cart is None:
                return state
            return list(self.state2cart(state[:6])) + [state[6]]

        if orbit_color is not None:
            if tof is None:
                tof = x[10 + 4 * self.nseg]
            for state, duration in ((self.leg.states[:7], tof), (self.leg.states[-7:], -tof)):
                self.leg.ta.state[:] = state
                self.leg.ta.pars[0] = 0.0
                self.leg.ta.time = 0.0
                samples = self.leg.ta.propagate_grid(_np.linspace(0.0, duration, N * self.nseg))[-1]
                cart = _np.asarray([to_cartesian(sample) for sample in samples])
                ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], color=orbit_color, alpha=0.5)

        throttles = _np.asarray(x[10 : 10 + 4 * self.nseg : 4])
        return _zoh_plot.plot_low_thrust_segments(
            self.leg,
            throttles,
            ax=ax,
            N=N,
            mark_segments=mark_segments,
            mark_mismatch=mark_mismatch,
            state_to_cartesian=to_cartesian,
            **kwargs,
        )

    def plot_throttle(self, x, ax=None, **kwargs):
        """Plot the normalized piecewise-constant throttle profile.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.
            ax (:class:`matplotlib.axes.Axes`, optional): Axes for the plot. A new axes
                is created when omitted.
            **kwargs: Additional keyword arguments passed to :meth:`matplotlib.axes.Axes.step`.

        Returns:
            :class:`matplotlib.axes.Axes`: Axes containing the throttle profile.
        """
        self._set_leg_from_x(x)
        throttles = _np.asarray(x[10 : 10 + 4 * self.nseg : 4])
        return _zoh_plot.plot_throttle_profile(self.leg.tgrid, throttles, ax=ax, **kwargs)