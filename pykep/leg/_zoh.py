import numpy as _np
import heyoka as _hy


class zoh:
    """Generic zero-order-hold trajectory leg.

    This class propagates a state of dimension ``dim_dynamics`` using piecewise-constant
    controls of size ``dim_controls`` over a user-supplied non-uniform time grid.

    A transfer is feasible when the mismatch constraints are zero. Any additional
    constraints on controls are intentionally left to the calling UDP.
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
        self.ta = tas[0]
        self.ta_var = tas[1]

        # Save non-control parameter values for cfunc calls
        self.pars_no_control = self.ta.pars[self.dim_controls :].tolist()

        # Convenient quantities
        self.nseg = len(self.controls) // self.dim_controls
        self.nseg_fwd = int(self.nseg * cut)
        self.nseg_bck = self.nseg - self.nseg_fwd

        # Sanity checks on integrators
        if len(self.ta.state) != self.dim_dynamics:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator state dimension of {len(self.ta.state)}, while {self.dim_dynamics} is required"
            )
        if len(self.ta.pars) < self.dim_controls:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator parameters dimension of {len(self.ta.pars)}, while >={self.dim_controls} is required"
            )

        if self.ta_var is not None:
            expected_var_dim = (
                self.dim_dynamics
                + self.dim_dynamics * self.dim_dynamics
                + self.dim_dynamics * self.dim_controls
            )
            if len(self.ta_var.state) != expected_var_dim:
                raise ValueError(
                    f"Attempting to construct a zoh_leg with a variational Taylor Adaptive integrator state dimension of {len(self.ta_var.state)}, while {expected_var_dim} is required"
                )
            if len(self.ta_var.pars) != len(self.ta.pars):
                raise ValueError(
                    "While constructing a zoh_leg, the number of parameters in the variational and non-variational integrators must be equal"
                )
            self.ic_var = _np.hstack(
                (
                    _np.eye(self.dim_dynamics, self.dim_dynamics),
                    _np.zeros((self.dim_dynamics, self.dim_controls)),
                )
            ).flatten()

        # Sanity checks on grids/controls
        if len(self.controls) % self.dim_controls > 0:
            raise ValueError(
                f"In a zoh_leg controls must be multiple of {self.dim_controls}: control * nseg"
            )
        if len(tgrid) != self.nseg + 1:
            raise ValueError(
                f"The t_grid and the controls have incompatible lengths. It must be nseg*{self.dim_controls} and nseg+1"
            )

        # Compile dynamics cfunc used in gradient computations
        sys = self.ta.sys
        vars = [it[0] for it in sys]
        dyn = [it[1] for it in sys]
        self.dyn_cfunc = _hy.cfunc_dbl(dyn, vars, compact_mode=True)

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

    def get_state_info(self, N=2):
        """Returns sampled state histories on each forward/backward segment.

        Returns:
            tuple[list, list, bool]: ``(state_fwd, state_bck, success)``.
        """
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
