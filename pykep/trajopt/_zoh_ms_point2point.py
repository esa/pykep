import numpy as _np
import pykep as _pk
from matplotlib import pyplot as _plt


class zoh_ms_point2point:
    """Optimal low-thrust transfer between fixed states using multiple shooting.

    The decision vector is ``[mf] + controls + [tof] + weights + nodes``,
    where weights are present only for the softmax time encoding and nodes are
    the seven-component interior states in chronological order. Each control
    is ``[throttle, ix, iy, iz]``; the throttle is scaled by ``max_thrust``.
    The objective is to maximize final mass. Segment defects are equalities,
    and direction norms are equalities or inequalities as in
    :class:`~pykep.trajopt.zoh_point2point`.

    Args:
        states (:class:`list`, optional): Initial six-component dynamics state.
            Default is ``[1.2, 0, -0.01, 0.01, 1, -0.01]``.
        statef (:class:`list`, optional): Final six-component dynamics state.
            Default is ``[1, 0, 0, 0.01, 1.1, 0]``.
        ms (:class:`float`, optional): Initial mass. Default is 1.
        max_thrust (:class:`float`, optional): Maximum thrust. Default is 0.22.
        tof_bounds (:class:`list`, optional): Time-of-flight bounds. Default is
            ``[3.4, 8.6]``.
        mf_bounds (:class:`list`, optional): Final-mass bounds. Default is
            ``[0.2, 1]``.
        nseg (:class:`int`, optional): Number of segments. Default is 10.
        cut (:class:`float`, optional): Fraction of segments propagated forward.
            Default is 0.6.
        tas (:class:`tuple`, optional): Nominal and variational integrators.
            A variational integrator is needed for gradients.
        time_encoding (:class:`str`, optional): ``'uniform'`` or ``'softmax'``.
            Default is ``'uniform'``.
        w_bounds_softmax (:class:`list`, optional): Softmax weight bounds.
            Default is ``[-1, 1]``.
        inequalities_for_tc (:class:`bool`, optional): Use direction-norm
            inequalities instead of equalities. Default is False.
        max_steps (:class:`int`, optional): Maximum integration steps per segment.
            Default is None.
        state2cart (:class:`callable`, optional): Convert a six-component
            dynamics state to Cartesian position and velocity for plotting.
        state_bounds (:class:`tuple`, optional): Lower and upper bounds for the
            six non-mass components of every interior node, in integrator units.
            Default is ``([-10] * 6, [10] * 6)``; choose bounds appropriate for
            the dynamics.
        interior_mass_bounds (:class:`tuple`, optional): Positive lower and upper
            mass bounds for every interior node, in integrator units. Default is
            ``(0.1 * ms, 4.0 * ms)``.
    """

    def __init__(
        self,
        states=None,
        statef=None,
        ms=1.0,
        max_thrust=0.22,
        tof_bounds=None,
        mf_bounds=None,
        nseg=10,
        cut=0.6,
        tas=(_pk.ta.get_zoh_kep(1e-10), None),
        time_encoding="uniform",
        w_bounds_softmax=None,
        inequalities_for_tc=False,
        max_steps=None,
        state2cart=None,
        state_bounds=None,
        interior_mass_bounds=None,
    ):
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")
        if time_encoding not in ("uniform", "softmax"):
            raise NotImplementedError("Only uniform and softmax time encodings are supported")
        if state2cart is not None and not callable(state2cart):
            raise ValueError("state2cart must be a callable function")

        states = [1.2, 0.0, -0.01, 0.01, 1.0, -0.01] if states is None else list(states)
        statef = [1.0, 0.0, 0.0, 0.01, 1.1, 0.0] if statef is None else list(statef)
        if len(states) != 6 or len(statef) != 6:
            raise ValueError("states and statef must each contain six components")
        if state_bounds is None:
            state_bounds = ([-100.0] * 6, [100.0] * 6)
        bounds = _np.asarray(state_bounds, dtype=float)
        if bounds.shape != (2, 6) or not _np.all(_np.isfinite(bounds)) or _np.any(bounds[0] > bounds[1]):
            raise ValueError("state_bounds must contain finite, ordered six-component lower and upper bounds")
        if interior_mass_bounds is None:
            interior_mass_bounds = (0.1 * ms, 4.0 * ms)
        mass_bounds = _np.asarray(interior_mass_bounds, dtype=float)
        if (mass_bounds.shape != (2,) or not _np.all(_np.isfinite(mass_bounds))
                or mass_bounds[0] <= 0.0 or mass_bounds[0] > mass_bounds[1]):
            raise ValueError("interior_mass_bounds must contain positive, finite, ordered bounds")

        self.nseg = nseg
        self.tof_bounds = [3.4, 8.6] if tof_bounds is None else tof_bounds
        self.mf_bounds = [0.2, 1.0] if mf_bounds is None else mf_bounds
        self.max_thrust = max_thrust
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = [-1.0, 1.0] if w_bounds_softmax is None else w_bounds_softmax
        self.inequalities_for_tc = inequalities_for_tc
        self.max_steps = max_steps
        self.state2cart = state2cart
        self.state_bounds = (bounds[0].tolist(), bounds[1].tolist())
        self.interior_mass_bounds = mass_bounds.tolist()
        self._node_lb = self.state_bounds[0] + [mass_bounds[0]]
        self._node_ub = self.state_bounds[1] + [mass_bounds[1]]
        self._state0 = states + [ms]
        self._statef = statef
        self._single_shooting_size = 2 + 4 * nseg + (nseg if time_encoding == "softmax" else 0)
        self._dtgrid_dtof_uniform = _np.linspace(0.0, 1.0, nseg + 1)
        self._softmax_cumsum_matrix = _np.vstack(
            [_np.zeros(nseg), _np.tril(_np.ones((nseg, nseg)))]
        )

        controls = _np.random.uniform(-1.0, 1.0, 4 * nseg)
        controls[0::4] = _np.abs(controls[0::4]) * max_thrust
        tgrid = _np.linspace(0.0, sum(self.tof_bounds) / 2.0, nseg + 1)
        node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
        self.leg = _pk.leg.zoh_ms(
            self._state0 + node * (nseg - 1) + self._statef + [ms],
            controls.tolist(), tgrid, cut, tas, max_steps=max_steps,
        )
        self._gradient_pattern = self._make_gradient_sparsity()
        self._gradient_index = {pair: index for index, pair in enumerate(self._gradient_pattern)}

    def _make_gradient_sparsity(self):
        pattern = {(0, 0)}
        tof_col = 1 + 4 * self.nseg
        for segment in range(self.nseg):
            for component in range(7):
                row = 1 + 7 * segment + component
                for node in (segment, segment + 1):
                    if 0 < node < self.nseg:
                        for state_component in range(7):
                            pattern.add((row, self._single_shooting_size + 7 * (node - 1) + state_component))
                    elif node == self.nseg:
                        pattern.add((row, 0))
                for control in range(4):
                    pattern.add((row, 1 + 4 * segment + control))
                pattern.add((row, tof_col))
                if self.time_encoding == "softmax":
                    for weight in range(self.nseg):
                        pattern.add((row, tof_col + 1 + weight))
            for direction in range(1, 4):
                pattern.add((1 + 7 * self.nseg + segment, 1 + 4 * segment + direction))
        return sorted(pattern)

    def _set_leg_from_x(self, x):
        if len(x) != self._single_shooting_size + 7 * (self.nseg - 1):
            raise ValueError("Invalid decision vector length")
        self.leg.states = self._state0 + list(x[self._single_shooting_size:]) + self._statef + [x[0]]
        controls = list(x[1 : 1 + 4 * self.nseg])
        controls[0::4] = [value * self.max_thrust for value in controls[0::4]]
        self.leg.controls = controls
        tof = x[1 + 4 * self.nseg]
        if self.time_encoding == "uniform":
            tgrid = tof * self._dtgrid_dtof_uniform
        else:
            weights = x[2 + 4 * self.nseg : self._single_shooting_size]
            fractions, _ = _pk.compute_softmax_and_jacobian(weights)
            tgrid = _np.r_[0.0, _np.cumsum(tof * fractions)]
        self.leg.tgrid = tgrid.tolist()

    def get_bounds(self):
        lb = [self.mf_bounds[0]] + [0.0, -1.0, -1.0, -1.0] * self.nseg + [self.tof_bounds[0]]
        ub = [self.mf_bounds[1]] + [1.0, 1.0, 1.0, 1.0] * self.nseg + [self.tof_bounds[1]]
        if self.time_encoding == "softmax":
            lb += [self.w_bounds_softmax[0]] * self.nseg
            ub += [self.w_bounds_softmax[1]] * self.nseg
        return lb + self._node_lb * (self.nseg - 1), ub + self._node_ub * (self.nseg - 1)

    def guess_internal_states(self, x, ballistic=False):
        """Return a complete decision vector with guessed interior states.

        Args:
            x (:class:`list`): A counterpart ZOH decision vector or a complete
                multiple-shooting vector. The prefix is preserved.
            ballistic (:class:`bool`, optional): Propagate with zero controls when
                constructing the nodes. Default is False.

        Raises:
            ValueError: If the vector length is invalid or a propagated interior
                node is outside its bounds. Widen `state_bounds` or
                `interior_mass_bounds` if needed.

        Returns:
            :class:`list`: Complete multiple-shooting decision vector.
        """
        x = list(x)
        full_size = self._single_shooting_size + 7 * (self.nseg - 1)
        if len(x) == self._single_shooting_size:
            node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
            x += node * (self.nseg - 1)
        elif len(x) != full_size:
            raise ValueError("Invalid decision vector length")
        self._set_leg_from_x(x)
        self.leg.set_initial_guess(ballistic=ballistic)
        nodes = self.leg.states[7:-7]
        if _np.any(_np.asarray(nodes) < self._node_lb * (self.nseg - 1)) or _np.any(
            _np.asarray(nodes) > self._node_ub * (self.nseg - 1)
        ) or not _np.all(_np.isfinite(nodes)):
            raise ValueError("Propagated interior nodes are outside their bounds or non-finite")
        return x[:self._single_shooting_size] + nodes

    def _compute_throttle_constraints(self):
        controls = _np.asarray(self.leg.controls).reshape(self.nseg, 4)
        return (_np.sum(controls[:, 1:] ** 2, axis=1) - 1.0).tolist()

    def compute_tc_grad(self):
        """Return the throttle-constraint Jacobian with respect to leg controls."""
        controls = _np.asarray(self.leg.controls).reshape(self.nseg, 4)
        result = _np.zeros((self.nseg, 4 * self.nseg))
        for segment in range(self.nseg):
            result[segment, 4 * segment + 1 : 4 * segment + 4] = 2.0 * controls[segment, 1:]
        return result

    def fitness(self, x):
        self._set_leg_from_x(x)
        return [-x[0]] + self.leg.compute_defects() + self._compute_throttle_constraints()

    def gradient_sparsity(self):
        """Return the sorted pygmo fitness-Jacobian sparsity pattern."""
        return list(self._gradient_pattern)

    def gradient(self, x):
        """Return fitness derivatives in the order of `gradient_sparsity`."""
        if not self.with_gradient:
            raise RuntimeError("Gradient computation requires a variational integrator")
        self._set_leg_from_x(x)
        values = _np.zeros(len(self._gradient_pattern))
        values[self._gradient_index[(0, 0)]] = -1.0
        patterns = self.leg.defects_grad_sparsity()
        derivatives = self.leg.compute_defects_grad()
        for (row, column), value in zip(patterns[0], derivatives[0]):
            node, component = divmod(int(column), 7)
            if 0 < node < self.nseg:
                dest = self._single_shooting_size + 7 * (node - 1) + component
            elif node == self.nseg and component == 6:
                dest = 0
            else:
                continue
            values[self._gradient_index[(int(row) + 1, dest)]] += value
        for (row, column), value in zip(patterns[1], derivatives[1]):
            column = int(column)
            scale = self.max_thrust if column % 4 == 0 else 1.0
            values[self._gradient_index[(int(row) + 1, column + 1)]] += value * scale

        tof_col = 1 + 4 * self.nseg
        if self.time_encoding == "uniform":
            dtgrid_dtof = self._dtgrid_dtof_uniform
        else:
            fractions, jacobian = _pk.compute_softmax_and_jacobian(x[tof_col + 1 : self._single_shooting_size])
            dtgrid_dtof = _np.r_[0.0, _np.cumsum(fractions)]
            dtgrid_dw = x[tof_col] * (self._softmax_cumsum_matrix @ jacobian)
        for (row, column), value in zip(patterns[2], derivatives[2]):
            row, column = int(row) + 1, int(column)
            values[self._gradient_index[(row, tof_col)]] += value * dtgrid_dtof[column]
            if self.time_encoding == "softmax":
                for weight in range(self.nseg):
                    values[self._gradient_index[(row, tof_col + 1 + weight)]] += value * dtgrid_dw[column, weight]
        for segment in range(self.nseg):
            row = 1 + 7 * self.nseg + segment
            for direction in range(1, 4):
                column = 1 + 4 * segment + direction
                values[self._gradient_index[(row, column)]] = 2.0 * x[column]
        return values

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        return 7 * self.nseg + (0 if self.inequalities_for_tc else self.nseg)

    def get_nic(self):
        return self.nseg if self.inequalities_for_tc else 0

    def pretty(self, x):
        """Print the leg corresponding to the decision vector."""
        self._set_leg_from_x(x)
        print(self.leg)

    def plot(self, x, ax=None, N=30, mark_segments=True, mark_mismatch=True,
             orbit_color=None, tof=None, **kwargs):
        """Plot the propagated segments and optionally their boundary orbits.

        Args:
            x (:class:`list`): Multiple-shooting decision vector.
            ax (:class:`matplotlib.axes.Axes`, optional): Three-dimensional axes.
            N (:class:`int`, optional): Number of samples per segment. Default is 30.
            mark_segments (:class:`bool`, optional): Mark segment starts. Default is True.
            mark_mismatch (:class:`bool`, optional): Mark propagated defect endpoints.
                Default is True.
            orbit_color (:class:`str`, optional): Color of the ballistic boundary orbits.
            tof (:class:`float`, optional): Boundary-orbit duration; defaults to the
                decision vector's time of flight.
            **kwargs: Additional arguments for the segment markers.

        Returns:
            :class:`matplotlib.axes.Axes`: The trajectory axes.
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
                tof = x[1 + 4 * self.nseg]
            ta = self.leg.ta
            for state, duration in ((self.leg.states[:7], tof), (self.leg.states[-7:], -tof)):
                ta.state[:] = state
                ta.pars[0] = 0.0
                ta.time = 0.0
                samples = ta.propagate_grid(_np.linspace(0.0, duration, N * self.nseg))[-1]
                if len(samples):
                    cart = _np.asarray([to_cartesian(sample) for sample in samples])
                    ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], color=orbit_color, alpha=0.5)

        forward, backward, _ = self.leg.get_state_info(N=N)
        throttles = _np.asarray(x[1 : 1 + 4 * self.nseg : 4])
        for segments, indices in ((forward, range(self.leg.nseg_fwd)),
                                  (backward, range(self.nseg - 1, self.leg.nseg_fwd - 1, -1))):
            for segment, index in zip(segments, indices):
                if not len(segment):
                    continue
                cart = _np.asarray([to_cartesian(state) for state in segment])
                throttle = throttles[index]
                color = (0.25 + 0.55 * throttle, 0.41 - 0.05 * throttle, 0.88 - 0.52 * throttle)
                ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], c=color)
                if mark_segments:
                    ax.scatter(*cart[0, :3], **kwargs)
                if mark_mismatch:
                    ax.scatter(*cart[-1, :3], marker="^", **kwargs)
        return ax

    def plot_throttle(self, x, ax=None, **kwargs):
        """Plot the normalized piecewise-constant throttle profile."""
        if ax is None:
            _, ax = _plt.subplots()
        self._set_leg_from_x(x)
        throttles = _np.asarray(x[1 : 1 + 4 * self.nseg : 4])
        ax.step(self.leg.tgrid, _np.r_[throttles, throttles[-1]], where="post", **kwargs)
        ax.set_xlabel("time grid")
        ax.set_ylabel("throttle value (nd)")
        return ax