import numpy as _np
import pykep as _pk
from . import _zoh_udp_plot_utils as _zoh_plot


class zoh_ms_ss_point2point:
    """Represents a fixed-endpoint solar-sail transfer using multiple shooting.

    The decision vector preserves the :class:`~pykep.trajopt.zoh_ss_point2point`
    prefix and appends the six-component interior shooting states::

        x = controls + [tof] (+ [weights]) + interior_states

    Each segment control is ``[alpha, beta]``: sail cone and clock angles in
    the RTN frame. The objective is time of flight and the constraints are the
    six-component defect at each segment boundary.

    Args:
        states (:class:`list`, optional): Initial six-component state. Defaults to
            ``[1.2, 0, -0.01, 0.01, 1, -0.01]``.
        statef (:class:`list`, optional): Final six-component state. Defaults to
            ``[1, 0, 0, 0.01, 1.1, 0]``.
        tof_bounds (:class:`list`, optional): Time-of-flight bounds in integrator units.
            Defaults to ``[3.4, 8.6]``.
        nseg (:class:`int`, optional): Number of shooting segments. Default is 10.
        cut (:class:`float`, optional): Forward/backward split passed to
            :class:`~pykep.leg.zoh_ms`. Default is 0.6.
        tas (:class:`tuple`, optional): Nominal and variational solar-sail integrators.
            A variational integrator enables gradients.
        time_encoding (:class:`str`, optional): ``'uniform'`` or ``'softmax'`` segment
            durations. Default is ``'uniform'``.
        w_bounds_softmax (:class:`list`, optional): Softmax-weight bounds. Default is
            ``[-1, 1]``.
        max_steps (:class:`int`, optional): Maximum integrator steps per segment.
        state_bounds (:class:`tuple`, optional): Lower and upper bounds for each
            six-component interior state. Default is ``([-100] * 6, [100] * 6)``.
    """

    def __init__(
        self,
        states=None,
        statef=None,
        tof_bounds=None,
        nseg=10,
        cut=0.6,
        tas=(_pk.ta.get_zoh_ss(1e-10), None),
        time_encoding="uniform",
        w_bounds_softmax=None,
        max_steps=None,
        state_bounds=None,
    ):
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")
        if time_encoding not in ("uniform", "softmax"):
            raise NotImplementedError("Only uniform and softmax time encodings are supported")

        states = [1.2, 0.0, -0.01, 0.01, 1.0, -0.01] if states is None else list(states)
        statef = [1.0, 0.0, 0.0, 0.01, 1.1, 0.0] if statef is None else list(statef)
        if len(states) != 6 or len(statef) != 6:
            raise ValueError("states and statef must each contain six components")
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

        self.nseg = nseg
        self.tof_bounds = [3.4, 8.6] if tof_bounds is None else tof_bounds
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = [-1.0, 1.0] if w_bounds_softmax is None else w_bounds_softmax
        self.max_steps = max_steps
        self.state_bounds = (bounds[0].tolist(), bounds[1].tolist())
        self._single_shooting_size = 2 * nseg + 1 + (nseg if time_encoding == "softmax" else 0)
        self._dtgrid_dtof_uniform = _np.linspace(0.0, 1.0, nseg + 1)
        self._softmax_cumsum_matrix = _np.vstack(
            [_np.zeros(nseg), _np.tril(_np.ones((nseg, nseg)))]
        )
        self._node_lb, self._node_ub = self.state_bounds
        controls = _np.random.uniform(-_np.pi / 2, _np.pi / 2, 2 * nseg)
        tof_mid = sum(self.tof_bounds) / 2.0
        node = ((_np.asarray(self._node_lb) + self._node_ub) / 2.0).tolist()
        self.leg = _pk.leg.zoh_ms(
            states + node * (nseg - 1) + statef,
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

    def _set_leg_from_x(self, x):
        if len(x) != self._expected_nx():
            raise ValueError("Invalid decision vector length")
        self.leg.states = list(self.leg.states[:6]) + list(
            x[self._single_shooting_size :]
        ) + list(self.leg.states[-6:])
        self.leg.controls = list(x[: 2 * self.nseg])
        tof_col = 2 * self.nseg
        tof = x[tof_col]
        if self.time_encoding == "uniform":
            tgrid = tof * self._dtgrid_dtof_uniform
        else:
            weights = x[tof_col + 1 : self._single_shooting_size]
            fractions, _ = _pk.compute_softmax_and_jacobian(weights)
            tgrid = _np.r_[0.0, _np.cumsum(tof * fractions)]
        self.leg.tgrid = tgrid.tolist()

    def _make_gradient_sparsity(self):
        pattern = {(0, 2 * self.nseg)}
        tof_col = 2 * self.nseg
        for segment in range(self.nseg):
            for component in range(6):
                row = 1 + 6 * segment + component
                for node in (segment, segment + 1):
                    if 0 < node < self.nseg:
                        for state_component in range(6):
                            pattern.add(
                                (row, self._single_shooting_size + 6 * (node - 1) + state_component)
                            )
                for control in range(2):
                    pattern.add((row, 2 * segment + control))
                pattern.add((row, tof_col))
                if self.time_encoding == "softmax":
                    for weight in range(self.nseg):
                        pattern.add((row, tof_col + 1 + weight))
        return sorted(pattern)

    def get_bounds(self):
        """Return bounds in decision-vector order, including interior states."""
        lb = ([-_np.pi / 2 + 1e-3, 0.0] * self.nseg) + [self.tof_bounds[0]]
        ub = ([_np.pi / 2 - 1e-3, 2.0 * _np.pi] * self.nseg) + [self.tof_bounds[1]]
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
        """Evaluate time of flight and the multiple-shooting defects."""
        self._set_leg_from_x(x)
        return [x[2 * self.nseg]] + self.leg.compute_defects()

    def gradient_sparsity(self):
        """Return the sparse fitness-Jacobian pattern in pygmo ordering."""
        return list(self._gradient_pattern)

    def gradient(self, x):
        """Return sparse derivatives of the objective and segment defects.

        Raises:
            RuntimeError: If no variational integrator was supplied.
        """
        if not self.with_gradient:
            raise RuntimeError("Gradient computation requires a variational integrator")
        self._set_leg_from_x(x)
        values = _np.zeros(len(self._gradient_pattern))
        values[self._gradient_index[(0, 2 * self.nseg)]] = 1.0
        state_pattern, control_pattern, time_pattern = self.leg.defects_grad_sparsity()
        state_gradient, control_gradient, time_gradient = self.leg.compute_defects_grad()
        for (row, column), value in zip(state_pattern, state_gradient):
            node, component = divmod(int(column), 6)
            if 0 < node < self.nseg:
                dest = self._single_shooting_size + 6 * (node - 1) + component
                values[self._gradient_index[(int(row) + 1, dest)]] += value
        for (row, column), value in zip(control_pattern, control_gradient):
            values[self._gradient_index[(int(row) + 1, int(column))]] += value

        tof_col = 2 * self.nseg
        if self.time_encoding == "uniform":
            dtgrid_dtof = self._dtgrid_dtof_uniform
        else:
            fractions, jacobian = _pk.compute_softmax_and_jacobian(
                x[tof_col + 1 : self._single_shooting_size]
            )
            dtgrid_dtof = _np.r_[0.0, _np.cumsum(fractions)]
            dtgrid_dw = x[tof_col] * (self._softmax_cumsum_matrix @ jacobian)
        for (row, column), value in zip(time_pattern, time_gradient):
            row, column = int(row) + 1, int(column)
            values[self._gradient_index[(row, tof_col)]] += value * dtgrid_dtof[column]
            if self.time_encoding == "softmax":
                for weight in range(self.nseg):
                    dest = tof_col + 1 + weight
                    values[self._gradient_index[(row, dest)]] += value * dtgrid_dw[column, weight]
        return values

    def has_gradient(self):
        """Return whether variational dynamics are available."""
        return self.with_gradient

    def get_nec(self):
        """Return the number of equality constraints."""
        return 6 * self.nseg

    def get_nic(self):
        """Return the number of inequality constraints."""
        return 0

    def pretty(self, x):
        """Print the multiple-shooting leg corresponding to a decision vector."""
        self._set_leg_from_x(x)
        print(self.leg)

    def plot(
        self,
        x,
        ax=None,
        N=30,
        to_cartesian=lambda state: state,
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
            to_cartesian (:class:`callable`, optional): Convert a dynamics state to Cartesian coordinates.
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
        angles = _np.asarray(x[: 2 * self.nseg]).reshape(self.nseg, 2)

        groups = (
            (forward, range(self.leg.nseg_fwd)),
            (backward, range(self.nseg - 1, self.leg.nseg_fwd - 1, -1)),
        )
        for segments, indices in groups:
            for segment, index in zip(segments, indices):
                if not len(segment):
                    continue
                cart = _np.asarray([to_cartesian(state) for state in segment])
                ax.plot(cart[:, 0], cart[:, 1], cart[:, 2], c="tab:orange")
                if mark_segments:
                    ax.scatter(*cart[0, :3], **kwargs)
                if mark_mismatch:
                    ax.scatter(*cart[-1, :3], marker="^", **kwargs)
                alpha, beta = angles[index]
                if plot_sail and abs(abs(alpha) - _np.pi / 2) > 1e-3:
                    _zoh_plot.draw_solar_sail(
                        ax,
                        segment[len(segment) // 2],
                        alpha,
                        beta,
                        sail_size,
                        state_to_cartesian=to_cartesian,
                    )
        _zoh_plot.set_3d_axes_equal(ax)
        return ax