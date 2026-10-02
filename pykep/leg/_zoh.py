import copy as _copy

import numpy as _np
import heyoka as _hy

class zoh:
    """Generic zero-order-hold trajectory leg. (fwd-bck shooting)

    This class propagates a state of dimension ``dim_dynamics`` using piecewise-constant
    controls of size ``dim_controls`` over a user-supplied time grid and forms
    one single defect (mismatch constraint) in the middle.

    A transfer is feasible when all compnents of the mismatch constraints are zero.
    Any additional constraints on controls are intentionally left to the calling UDP.
    """

    def __init__(
        self,
        state0,
        controls,
        state1,
        tgrid,
        cut,
        tas, 
        max_steps=None,
        dim_dynamics=7,
        dim_controls=4,
    ):
        """zoh(state0, controls, state1, tgrid, cut, tas, max_steps=None, dim_dynamics=7, dim_controls=4)

        Args:
            state0 (:class:`list`): Initial state, of length ``dim_dynamics``.

            controls (:class:`list`): Flat vector of the piecewise-constant controls, of length ``nseg * dim_controls``.

            state1 (:class:`list`): Final state, of length ``dim_dynamics``.

            tgrid (:class:`list`): Time grid of ``nseg + 1`` points.

            cut (:class:`float`): Fraction of segments, in :math:`[0, 1]`, propagated forward.

            tas (:class:`tuple`): Pair ``(ta, ta_var)`` of :class:`heyoka.taylor_adaptive` integrators. The first
            has ``dim_dynamics`` states and at least ``dim_controls`` parameters (controls first); the second is its
            variational counterpart w.r.t. states and controls, or None.

            max_steps (:class:`int`, optional): Maximum number of integration steps per segment. Default is None (no limit).

            dim_dynamics (:class:`int`, optional): Dimension of the state. Default is 7.

            dim_controls (:class:`int`, optional): Dimension of the control. Default is 4.

        Raises:
            ValueError: If the integrators, endpoint states, ``controls`` or ``tgrid`` have inconsistent dimensions.

        Notes:
            ``state0``, ``controls``, ``state1`` and ``tgrid`` are stored as passed and are meant to be lists,
            mirroring ``std::vector``. Other sequences, such as NumPy arrays, also work. The integrators in
            ``tas`` are deep-copied; propagation mutates the leg's copies, not the supplied integrators.
        """
        # We store the constructor args
        self.state0 = state0
        self.controls = controls
        self.state1 = state1
        self.tgrid = tgrid
        self.cut = cut
        self.max_steps = max_steps
        self.dim_dynamics = dim_dynamics
        self.dim_controls = dim_controls

        # Store the integrators
        self.ta = _copy.deepcopy(tas[0])
        self.ta_var = _copy.deepcopy(tas[1])

        # Convenient quantities
        if self.dim_controls <= 0:
            raise ValueError("The control dimension must be positive")
        self.nseg = len(self.controls) // self.dim_controls
        self.nseg_fwd = int(self.nseg * cut)
        self.nseg_bck = self.nseg - self.nseg_fwd

        # Guard against mutations to the leg data before using it.
        self._validate_leg()

        # Save non-control parameter values for cfunc calls
        self.pars_no_control = self.ta.pars[self.dim_controls :].tolist()

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
        """Validate integrator dimensions and the mutable leg data."""
        if self.dim_dynamics <= 0:
            raise ValueError("The dynamics dimension must be positive")
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
        if len(self.state0) != self.dim_dynamics:
            raise ValueError(
                f"Attempting to use a zoh_leg with an initial state of length {len(self.state0)}, while {self.dim_dynamics} is required"
            )
        if len(self.state1) != self.dim_dynamics:
            raise ValueError(
                f"Attempting to use a zoh_leg with a final state of length {len(self.state1)}, while {self.dim_dynamics} is required"
            )
        if len(self.tgrid) != self.nseg + 1:
            raise ValueError(
                f"Attempting to use a zoh_leg with a time grid of length {len(self.tgrid)}, while {self.nseg + 1} is required"
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

    def compute_mismatch_constraints(self):
        """Propagates forward/backward and returns the state mismatch at the midpoint.

        Returns:
            :class:`list`: Mismatch vector of length ``dim_dynamics``.
        """
        # Guard against mutations to the leg data before using it.
        self._validate_leg()
        c = self.dim_controls

        # Forward segments (up to failure or cut)
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            start = c * i
            self.ta.pars[:c] = self.controls[start : start + c]
            if not self._propagate_until(self.ta, self.tgrid[i + 1]):
                break
        state_fwd = self.ta.state.copy()

        # Backward segments (up to failure or cut)
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            start = c * (self.nseg - 1 - i)
            self.ta.pars[:c] = self.controls[start : start + c]
            if not self._propagate_until(self.ta, self.tgrid[-2 - i]):
                break

        state_bck = self.ta.state
        return (state_fwd - state_bck).tolist()

    def compute_mc_grad(self):
        """Computes gradients of mismatch constraints.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
            ``(dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid)`` with shapes
            ``(d,d)``, ``(d,d)``, ``(d,c*nseg)``, ``(d,nseg+1)``.
        """
        # Guard against mutations to the leg data before using it.
        self._validate_leg()
        if self.ta_var is None:
            raise RuntimeError(
                "compute_mc_grad requires a variational integrator (tas[1] must not be None)"
            )

        d = self.dim_dynamics
        c = self.dim_controls
        stm_cols = d + c

        # STMs -> forward
        M_seg_fwd = []
        M_fwd = []
        C_seg_fwd = []
        C_fwd = _np.zeros((d, c * self.nseg_fwd))
        dyn_fwd = []

        self.ta_var.time = self.tgrid[0]
        self.ta_var.state[:d] = self.state0
        successful_fwd = 0
        for i in range(self.nseg_fwd):
            self.ta_var.state[d:] = self.ic_var
            start = c * i
            self.ta_var.pars[:c] = self.controls[start : start + c]
            if not self._propagate_until(self.ta_var, self.tgrid[i + 1]):
                break

            seg_stm = self.ta_var.state[d:].reshape(d, stm_cols)
            M_seg_fwd.append(seg_stm[:, :d].copy())
            C_seg_fwd.append(seg_stm[:, d:].copy())
            dyn_fwd.append(
                self.dyn_cfunc(
                    self.ta_var.state[:d],
                    pars=[*self.controls[start : start + c], *self.pars_no_control],
                )
            )
            successful_fwd += 1

        cur = _np.eye(d)
        for M in reversed(M_seg_fwd):
            cur = cur @ M
            M_fwd.append(cur)
        M_fwd = list(reversed(M_fwd)) + [_np.eye(d)]

        dmc_dx0 = M_fwd[0]

        i = 0
        for M, C in zip(M_fwd[1:], C_seg_fwd):
            C_fwd[:, c * i : c * i + c] = M @ C
            i += 1

        dmcdtgrid = _np.zeros((d, self.nseg + 1))
        if successful_fwd > 0:
            dmcdtgrid[:, 0] = -M_fwd[1] @ dyn_fwd[0]
            dmcdtgrid[:, successful_fwd] = M_fwd[-1] @ dyn_fwd[-1]
            for i in range(1, successful_fwd):
                dmcdtgrid[:, i] = M_fwd[i + 1] @ (
                    M_seg_fwd[i] @ dyn_fwd[i - 1] - dyn_fwd[i]
                )

        # STMs -> backward
        M_seg_bck = []
        M_bck = []
        C_seg_bck = []
        C_bck = _np.zeros((d, c * self.nseg_bck))
        dyn_bck = []

        self.ta_var.time = self.tgrid[-1]
        self.ta_var.state[:d] = self.state1
        successful_bck = 0
        for i in range(self.nseg_bck):
            self.ta_var.state[d:] = self.ic_var
            start = c * (self.nseg - 1 - i)
            self.ta_var.pars[:c] = self.controls[start : start + c]
            if not self._propagate_until(self.ta_var, self.tgrid[-2 - i]):
                break

            seg_stm = self.ta_var.state[d:].reshape(d, stm_cols)
            M_seg_bck.append(seg_stm[:, :d].copy())
            C_seg_bck.append(seg_stm[:, d:].copy())
            dyn_bck.append(
                self.dyn_cfunc(
                    self.ta_var.state[:d],
                    pars=[*self.controls[start : start + c], *self.pars_no_control],
                )
            )
            successful_bck += 1

        cur = _np.eye(d)
        for M in reversed(M_seg_bck):
            cur = cur @ M
            M_bck.append(cur)
        M_bck = list(reversed(M_bck)) + [_np.eye(d)]

        dmc_dx1 = -M_bck[0]

        i = 0
        for M, C in zip(M_bck[1:], C_seg_bck):
            start = c * self.nseg_bck - c * (i + 1)
            C_bck[:, start : start + c] = M @ C
            i += 1

        if successful_bck > 0:
            dmcdtgrid[:, -1] = M_bck[1] @ dyn_bck[0]
            dmcdtgrid[:, self.nseg - successful_bck] -= M_bck[-1] @ dyn_bck[-1]
            for i in range(1, successful_bck):
                dmcdtgrid[:, -1 - i] = -M_bck[i + 1] @ (
                    M_seg_bck[i] @ dyn_bck[i - 1] - dyn_bck[i]
                )

        dmc_dcontrols = _np.hstack((C_fwd, -C_bck))
        return dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid

    def get_state_info(self, N=50):
        """Returns sampled state histories on each forward and backward segment.

        Each segment contains ``N`` states, including its endpoints. Forward
        segments are ordered from the initial time toward the cut, and their
        samples are in increasing-time order. Backward segments are ordered
        from the final time toward the cut, and their samples are in
        decreasing-time order. The state components retain the ordering used
        by the integrator.

        Args:
            N (:class:`int`, optional): Number of sampling points per segment,
                including both endpoints. Default is 50.

        Returns:
            tuple[list, list, bool]: ``(state_fwd, state_bck, success)``.
            Each history is a list of per-segment state arrays; either list may
            be shorter than its expected length if propagation fails.
            ``success`` is ``True`` only if every requested segment returns all
            ``N`` samples.
        """
        # Guard against mutations to the leg data before using it.
        self._validate_leg()
        c = self.dim_controls
        success = True

        state_fwd = []
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            start = c * i
            self.ta.pars[:c] = self.controls[start : start + c]
            plot_grid_fwd = _np.linspace(self.tgrid[i], self.tgrid[i + 1], N)
            if self.max_steps is not None:
                sol_fwd = self.ta.propagate_grid(plot_grid_fwd, max_steps=self.max_steps)[-1]
            else:
                sol_fwd = self.ta.propagate_grid(plot_grid_fwd)[-1]
            state_fwd.append(sol_fwd)
            if len(sol_fwd) < N:
                success = False
                break

        state_bck = []
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            start = c * (self.nseg - 1 - i)
            self.ta.pars[:c] = self.controls[start : start + c]
            plot_grid_bck = _np.linspace(self.tgrid[-1 - i], self.tgrid[-2 - i], N)
            if self.max_steps is not None:
                sol_bck = self.ta.propagate_grid(plot_grid_bck, max_steps=self.max_steps)[-1]
            else:
                sol_bck = self.ta.propagate_grid(plot_grid_bck)[-1]
            state_bck.append(sol_bck)
            if len(sol_bck) < N:
                success = False
                break

        return state_fwd, state_bck, success
