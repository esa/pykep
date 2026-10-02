import copy as _copy

import numpy as _np
import heyoka as _hy


class zoh_ms:
    """Generic zero-order-hold trajectory leg. (multiple shooting)

    This class propagates a state of dimension ``dim_dynamics`` using piecewise-constant
    controls of size ``dim_controls`` over a user-supplied time grid and forms
    defects at each segment boundary, but preserving the fwd-bck logic.

    A transfer is feasible when all defects are zero.
    Any additional constraints on controls are intentionally left to the calling code.
    """

    def __init__(
        self,
        states,
        controls,
        tgrid,
        cut,
        tas, 
        max_steps=None,
        dim_dynamics=7,
        dim_controls=4,
    ):
        """zoh_ms(states, controls, tgrid, cut, tas, max_steps=None, dim_dynamics=7, dim_controls=4)

        Args:
            states (:class:`list`): Flat vector of the ``nseg + 1`` nodes, of length ``(nseg + 1) * dim_dynamics``.

            controls (:class:`list`): Flat vector of the piecewise-constant controls, of length ``nseg * dim_controls``.

            tgrid (:class:`list`): Time grid of ``nseg + 1`` points.

            cut (:class:`float`): Fraction of segments, in :math:`[0, 1]`, propagated forward.

            tas (:class:`tuple`): Pair ``(ta, ta_var)`` of :class:`heyoka.taylor_adaptive` integrators. The first
            has ``dim_dynamics`` states and at least ``dim_controls`` parameters (controls first); the second is its
            variational counterpart w.r.t. states and controls, or None.

            max_steps (:class:`int`, optional): Maximum number of integration steps per segment. Default is None (no limit).

            dim_dynamics (:class:`int`, optional): Dimension of the state. Default is 7.

            dim_controls (:class:`int`, optional): Dimension of the control. Default is 4.

        Raises:
            ValueError: If the integrators, ``states``, ``controls`` and ``tgrid`` have inconsistent dimensions.

        Notes:
            ``states``, ``controls`` and ``tgrid`` are stored as passed and are meant to be lists, mirroring
            ``std::vector``. Other sequences, such as NumPy arrays, also work, but methods that write
            these attributes (e.g. :meth:`set_initial_guess`) store lists. The integrators in ``tas`` are
            deep-copied; propagation mutates the leg's copies, not the supplied integrators.
        """
        # We store the constructor args
        self.states = states
        self.controls = controls
        self.tgrid = tgrid
        self.cut = cut
        self.max_steps = max_steps
        self.dim_dynamics = dim_dynamics
        self.dim_controls = dim_controls

        # Store the integrators
        self.ta = _copy.deepcopy(tas[0])
        self.ta_var = _copy.deepcopy(tas[1])

        # Save non-control parameter values for cfunc calls
        self.pars_no_control = self.ta.pars[self.dim_controls :].tolist()

        # Convenient quantities
        self.nseg = len(self.controls) // self.dim_controls
        self.nseg_fwd = int(self.nseg * cut)
        self.nseg_bck = self.nseg - self.nseg_fwd

        self._validate_leg()

        if self.ta_var is not None:
            self.ic_var = _np.hstack(
                (
                    _np.eye(self.dim_dynamics, self.dim_dynamics),
                    _np.zeros((self.dim_dynamics, self.dim_controls)),
                )
            ).flatten()

        # Compile dynamics cfunc used in gradient computations
        sys = self.ta.sys
        vars = [it[0] for it in sys]
        dyn = [it[1] for it in sys]
        self.dyn_cfunc = _hy.cfunc_dbl(dyn, vars, compact_mode=True)

    def _validate_leg(self):
        """Validate the dimensions of the mutable leg data."""
        if len(self.ta.state) != self.dim_dynamics:
            raise ValueError(
                f"Attempting to use a zoh_leg with a Taylor Adaptive integrator state dimension of {len(self.ta.state)}, while {self.dim_dynamics} is required"
            )
        if len(self.ta.pars) < self.dim_controls:
            raise ValueError(
                f"Attempting to use a zoh_leg with a Taylor Adaptive integrator parameters dimension of {len(self.ta.pars)}, while >={self.dim_controls} is required"
            )
        if self.ta_var is not None:
            expected_var_dim = (
                self.dim_dynamics
                + self.dim_dynamics * self.dim_dynamics
                + self.dim_dynamics * self.dim_controls
            )
            if len(self.ta_var.state) != expected_var_dim:
                raise ValueError(
                    f"Attempting to use a zoh_leg with a variational Taylor Adaptive integrator state dimension of {len(self.ta_var.state)}, while {expected_var_dim} is required"
                )
            if len(self.ta_var.pars) != len(self.ta.pars):
                raise ValueError(
                    "While using a zoh_leg, the number of parameters in the variational and non-variational integrators must be equal"
                )
        if len(self.controls) % self.dim_controls > 0:
            raise ValueError(
                f"Attempting to use a zoh_leg with a number of controls ({len(self.controls)}) that is not a multiple of the control dimension ({self.dim_controls})"
            )
        if len(self.controls) != self.nseg * self.dim_controls:
            raise ValueError(
                f"Attempting to use a zoh_leg with a controls array of length {len(self.controls)}, while {self.nseg * self.dim_controls} is required"
            )
        if len(self.tgrid) != self.nseg + 1:
            raise ValueError(
                f"Attempting to use a zoh_leg with a time grid of length {len(self.tgrid)}, while {self.nseg + 1} is required"
            )
        if len(self.states) != (self.nseg + 1) * self.dim_dynamics:
            raise ValueError(
                f"Attempting to use a zoh_leg with a states array of length {len(self.states)}, while {(self.nseg + 1) * self.dim_dynamics} is required"
            )

    def _propagate_until(self, ta, t_end):
        """Propagate safely and restore state if the integrator fails.
           NOTE: this behaviour should actually be implemented in heyoka directly
           for consistency with the propagate_grid one.
        """
        previous_time = ta.time
        previous_state = ta.state.copy()
        try:
            if self.max_steps is not None:
                ta.propagate_until(t_end, max_steps=self.max_steps)
            else:
                ta.propagate_until(t_end)
        except Exception:
            ta.time = previous_time
            ta.state[:] = previous_state
            return False
        return True

    def compute_defects(self):
        """Propagates each segment from its starting node and returns the defects.

        Segments with index :math:`i < n_{fwd}` are propagated forward from node :math:`i` to node
        :math:`i+1`, the remaining ones backward from node :math:`i+1` to node :math:`i`.
        Following the fwd-bck convention of :class:`~pykep.leg.zoh`, the defects are:

        .. math::
            \\mathbf d_i = \\begin{cases}
            \\boldsymbol\\phi_i(\\mathbf x_i) - \\mathbf x_{i+1} & i < n_{fwd} \\\\
            \\mathbf x_i - \\boldsymbol\\phi_i^{-1}(\\mathbf x_{i+1}) & i \\ge n_{fwd}
            \\end{cases}

        so that the initial node is never the target of a forward segment and the final node
        is never the target of a backward one.

        Notes:
            If the integration of a segment fails, the integrator state is restored to the
            segment starting node, which then enters the defect in place of the propagated state.

            Since the flow is invertible, :math:`\\boldsymbol\\phi_i(\\mathbf x_i) = \\mathbf x_{i+1}
            \\iff \\mathbf x_i = \\boldsymbol\\phi_i^{-1}(\\mathbf x_{i+1})`: the set of feasible
            trajectories does not depend on ``cut``. Away from feasibility, though, forward and
            backward defects differ (:math:`\\mathbf d_{bck} \\approx -\\Phi_i^{-1}\\mathbf d_{fwd}`),
            so ``cut`` affects constraint tolerances, the conditioning of the Jacobian seen by the
            solver and, in case of integration failures, whether a segment appears feasible.

        Returns:
            :class:`list`: Flattened defects vector of length ``nseg * dim_dynamics``, ordered by segment.
        """
        self._validate_leg()
        c = self.dim_controls
        d = self.dim_dynamics
        defects = []

        # Forward segments: node i -> node i+1
        for i in range(self.nseg_fwd):
            self.ta.time = self.tgrid[i]
            self.ta.state[:] = self.states[d * i : d * i + d]
            self.ta.pars[:c] = self.controls[c * i : c * i + c]
            self._propagate_until(self.ta, self.tgrid[i + 1])
            defects.append(self.ta.state - self.states[d * (i + 1) : d * (i + 1) + d])

        # Backward segments: node i+1 -> node i
        for i in range(self.nseg_fwd, self.nseg):
            self.ta.time = self.tgrid[i + 1]
            self.ta.state[:] = self.states[d * (i + 1) : d * (i + 1) + d]
            self.ta.pars[:c] = self.controls[c * i : c * i + c]
            self._propagate_until(self.ta, self.tgrid[i])
            defects.append(self.states[d * i : d * i + d] - self.ta.state)

        return _np.concatenate(defects).tolist()

    def defects_grad_sparsity(self):
        """Returns the sparsity patterns of the defects gradients.

        Each pattern lists the ``(row, col)`` indices of the structurally nonzero entries of one
        Jacobian block, in the order used by :meth:`compute_defects_grad` (pagmo convention).
        Rows index the defects, columns index ``states``, ``controls`` and ``tgrid`` respectively.

        For each segment, the node the propagation starts from enters through the full STM
        :math:`\\Phi_i`, the other one through :math:`\\pm\\mathbf I`, of which only the diagonal is stored.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: ``(sp_states, sp_controls, sp_tgrid)``,
            integer arrays of shapes ``(nseg * (d * d + d), 2)``, ``(nseg * d * c, 2)`` and ``(nseg * d * 2, 2)``.
        """
        d = self.dim_dynamics
        c = self.dim_controls
        sp_states, sp_controls, sp_tgrid = [], [], []

        # States: full STM row on the starting node, diagonal of +-I on the other
        for i in range(self.nseg):
            for k in range(d):
                r = d * i + k
                if i < self.nseg_fwd:
                    sp_states += [(r, d * i + j) for j in range(d)] + [(r, d * (i + 1) + k)]
                else:
                    sp_states += [(r, d * i + k)] + [(r, d * (i + 1) + j) for j in range(d)]

        # Controls: dense d x c block on the segment control
        for i in range(self.nseg):
            for k in range(d):
                sp_controls += [(d * i + k, c * i + j) for j in range(c)]

        # Tgrid: dense d x 2 block on the segment start and end times
        for i in range(self.nseg):
            for k in range(d):
                sp_tgrid += [(d * i + k, i), (d * i + k, i + 1)]

        return tuple(
            _np.array(sp, dtype=_np.int64).reshape(-1, 2)
            for sp in (sp_states, sp_controls, sp_tgrid)
        )

    def compute_defects_grad(self):
        """Computes the gradients of the defects w.r.t. ``states``, ``controls`` and ``tgrid``.

        Calling :math:`\\Phi_i` and :math:`\\mathbf C_i` the sensitivities of the segment flow
        w.r.t. its starting node and controls, and :math:`\\mathbf f` the dynamics, the nonzero blocks are:

        .. math::
            \\begin{array}{c|ccccc}
             & \\partial\\mathbf x_i & \\partial\\mathbf x_{i+1} & \\partial\\mathbf u_i & \\partial t_i & \\partial t_{i+1} \\\\
            \\hline
            i < n_{fwd} & \\Phi_i & -\\mathbf I & \\mathbf C_i & -\\Phi_i\\mathbf f(\\mathbf x_i) & \\mathbf f(\\boldsymbol\\phi_i) \\\\
            i \\ge n_{fwd} & \\mathbf I & -\\Phi_i & -\\mathbf C_i & -\\mathbf f(\\boldsymbol\\phi_i^{-1}) & \\Phi_i\\mathbf f(\\mathbf x_{i+1})
            \\end{array}

        Raises:
            RuntimeError: If no variational integrator was provided.

        Notes:
            If the integration of a segment fails, the segment defect does not depend on the flow
            and the corresponding blocks reduce to :math:`\\Phi_i = \\mathbf I`, :math:`\\mathbf C_i = 0`
            and zero time derivatives.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: ``(grad_states, grad_controls, grad_tgrid)``,
            flat arrays of nonzero values ordered as in :meth:`defects_grad_sparsity`.
        """
        if self.ta_var is None:
            raise RuntimeError(
                "compute_defects_grad requires a variational integrator (tas[1] must not be None)"
            )
        self._validate_leg()
        d = self.dim_dynamics
        c = self.dim_controls
        grad_states, grad_controls, grad_tgrid = [], [], []

        for i in range(self.nseg):
            fwd = i < self.nseg_fwd
            start, end = (i, i + 1) if fwd else (i + 1, i)
            x0 = _np.array(self.states[d * start : d * start + d], dtype=float)
            u = self.controls[c * i : c * i + c]
            pars = [*u, *self.pars_no_control]

            self.ta_var.time = self.tgrid[start]
            self.ta_var.state[:d] = x0
            self.ta_var.state[d:] = self.ic_var
            self.ta_var.pars[:c] = u
            success = self._propagate_until(self.ta_var, self.tgrid[end])

            stm = self.ta_var.state[d:].reshape(d, d + c).copy()
            M, C = stm[:, :d], stm[:, d:]

            # Flow derivatives w.r.t. its start and end times
            if success:
                dphi_dt0 = -M @ self.dyn_cfunc(x0, pars=pars)
                dphi_dt1 = self.dyn_cfunc(self.ta_var.state[:d], pars=pars)
            else:
                dphi_dt0 = dphi_dt1 = _np.zeros(d)

            if fwd:
                grad_states.append(_np.hstack((M, -_np.ones((d, 1)))))
                grad_controls.append(C)
                grad_tgrid.append(_np.column_stack((dphi_dt0, dphi_dt1)))
            else:
                grad_states.append(_np.hstack((_np.ones((d, 1)), -M)))
                grad_controls.append(-C)
                grad_tgrid.append(-_np.column_stack((dphi_dt1, dphi_dt0)))

        return tuple(
            _np.concatenate([g.ravel() for g in grads])
            for grads in (grad_states, grad_controls, grad_tgrid)
        )

    def set_initial_guess(self, ballistic=False):
        """Sets interior nodes by propagating inward from the current endpoint states.

        The endpoint states and time grid are read from ``self.states`` and ``self.tgrid``.
        Existing interior nodes are replaced, while the endpoints and time grid are preserved.

        Args:
            ballistic (:class:`bool`, optional): If True, propagate the initial-guess mesh with zero controls,
                without changing the leg's stored controls. Default is False.

        Raises:
            ValueError: If the controls, states, or time grid have inconsistent dimensions.

        Notes:
            With ``ballistic=True``, mismatches may occur at each segment. With ``ballistic=False``,
            only the segment where the forward and backward propagations meet may have a mismatch.
        """
        self._validate_leg()
        c = self.dim_controls
        d = self.dim_dynamics

        controls = [0.0] * (c * self.nseg) if ballistic else self.controls

        # Forward propagation fills nodes 1 through nseg_fwd.
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.states[:d]

        for i in range(min(self.nseg_fwd, self.nseg - 1)):
            self.ta.pars[:c] = controls[c * i : c * i + c]
            self._propagate_until(self.ta, self.tgrid[i + 1])
            self.states[d*(i + 1) : d*(i + 1) + d] = self.ta.state

        # Backward propagation fills nodes nseg-1 through nseg_fwd+1.
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.states[-d:]

        for i in range(self.nseg - 1, self.nseg_fwd, -1):
            self.ta.pars[:c] = controls[c * i : c * i + c]
            self._propagate_until(self.ta, self.tgrid[i])
            self.states[d*i : d*i + d] = self.ta.state

    def get_state_info(self, N=5):
        """Returns sampled state histories for each forward and backward segment.

        Each entry is an ``N`` by ``dim_dynamics`` array. Forward segments are
        ordered from node 0 toward the cut and sampled in increasing-time
        order. Backward segments are ordered from the final node toward the cut
        and sampled in decreasing-time order.

        Args:
            N (:class:`int`, optional): Number of sampling points per segment,
                including both endpoints. Default is 5.

        Returns:
            tuple[list, list, bool]: ``(state_fwd, state_bck, success)``.
            ``success`` is ``True`` only if every requested propagation returns
            all ``N`` samples. A failed propagation may return fewer samples.
        """
        c = self.dim_controls
        state_fwd = []
        state_bck = []
        success = True

        for i in range(self.nseg_fwd):
            start = self.dim_dynamics * i
            state = self.states[start : start + self.dim_dynamics]
            self.ta.time = self.tgrid[i]
            self.ta.state[:] = state
            self.ta.pars[:c] = self.controls[c * i : c * i + c]
            plot_grid = _np.linspace(self.tgrid[i], self.tgrid[i + 1], N)
            if self.max_steps is not None:
                sol = self.ta.propagate_grid(plot_grid, max_steps=self.max_steps)[-1]
            else:
                sol = self.ta.propagate_grid(plot_grid)[-1]
            state_fwd.append(sol)
            if len(sol) < N:
                success = False
                break

        for i in range(self.nseg - 1, self.nseg_fwd - 1, -1):
            start = self.dim_dynamics * (i + 1)
            state = self.states[start : start + self.dim_dynamics]
            self.ta.time = self.tgrid[i + 1]
            self.ta.state[:] = state
            self.ta.pars[:c] = self.controls[c * i : c * i + c]
            plot_grid = _np.linspace(self.tgrid[i + 1], self.tgrid[i], N)
            if self.max_steps is not None:
                sol = self.ta.propagate_grid(plot_grid, max_steps=self.max_steps)[-1]
            else:
                sol = self.ta.propagate_grid(plot_grid)[-1]
            state_bck.append(sol)
            if len(sol) < N:
                success = False
                break

        return state_fwd, state_bck, success

