import numpy as _np
import heyoka as _hy

class zoh:
    """
    This class implements an interplanetary low-thrust transfer between a starting and final state
    in the augmented state-space :math:`[\\mathbf{r}, \\mathbf{v}, m]`. The transfer is modelled
    as a sequence of non-uniform segments along which a continuous and constant (zero-order hold)
    control acts. The time intervals defining these segments are also provided in `tgrid`.

    The formulation generalises :class:`pykep.leg.sims_flanagan` to arbitrary dynamics and non-uniform
    time grids. The dynamics are assumed to be zero-order hold and must be provided as compatible
    Taylor-adaptive integrators (`tas`). 
    
    .. note::
       The requirements on the `tas` passed are: a) the first four *heyoka* parameters
       must be :math:`T, i_x, i_y, i_z`, b) the system dimension must be 7 c) for the variational 
       integrator, variations on the state and the four parameters only are considered. These 
       requirements are all fulfilled by :class:`pykep.ta.zoh_kep`, :class:`pykep.ta.zoh_eq`,
       :class:`pykep.ta.zoh_cr3bp` and their variational versions.

    A transfer is feasible when the state mismatch equality constraints are satisfied. In the
    intended usage, throttle equality constraints are also enforced to ensure a proper thrust
    representation as :math:`T \\hat{\\mathbf{i}}` with :math:`|\\hat{\\mathbf{i}}| = 1`.

    .. math::
       i_x^2 + i_y^2 + i_z^2 = 1, \\quad \\forall \\text{segments}
    """

    def __init__(
        self,
        state0,
        controls,
        state1,
        tgrid,
        cut,
        tas,
        max_steps = None
    ):
        """
        .. note::
           The variational integrator `ta_var` can be ``None`` or must have state dimension 84 (7 + 7x7 STM + 7x4 control
           sensitivity) and with the same dynamics as the nominal integrator `ta`. It is the user
           that must ensure the suitability of the integrators.

        Args:
            *state0* (:class:`array-like`): Initial state :math:`[\\mathbf{r}_0, \\mathbf{v}_0, m_0]` (length 7)
            
            *controls* (:class:`array-like`): Control parameters :math:`[T, i_x, i_y, i_z] \\times n_\\text{seg}`
            
            *state1* (:class:`array-like`): Final state :math:`[\\mathbf{r}_1, \\mathbf{v}_1, m_1]` (length 7)
            
            *tgrid* (:class:`array-like`): Non-uniform time grid (length nseg+1)
            
            *cut* (:class:`float`): Forward/backward segment split ratio (0 ≤ cut ≤ 1)
            
            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators

                - `ta`: Nominal dynamics (state dim 7, pars ≥ 4)
                
                - `ta_var`: Variational dynamics (state dim 84, same pars)
                
            *max_steps* (:class:`int`): Maximum number of steps for the integrator. If not ``None``, it is passed to the integrator before each propagation call.

        Raises:
            :class:`ValueError`: If state/parameter dimensions mismatch or input lengths are incompatible.

        Examples:
        
        .. code-block:: python 
        
           import numpy as np
           import heyoka as hy
           # Define states, controls, time grid
           state0 = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0])
           state1 = np.array([1.2, 0.1, 0.0, 0.0, 0.9, 0.1, 0.95])
           controls = np.array([0.022, 0.7, 0.7, 0.1, 0.025, -0.3, 0.8, 0.4, 0.015, -0.2, 0.8, 0.4])
           tgrid = np.array([0.0, 0.5, 1.0, 1.23])
           # Get integrators
           ta = pk.ta.get_zoh_eq(tol=1e-16)
           ta_var = pk.ta.get_zoh_eq_var(tol=1e-16)
           # Construct leg (50/50 split)
           leg = zoh(state0, controls, state1, tgrid, cut=0.5, tas=(ta, ta_var))
        """

        # We store the constructor args
        self.state0 = state0
        self.controls = controls
        self.state1 = state1
        self.tgrid = tgrid
        self.cut = cut
        self.max_steps = max_steps
        
        # Store the tas
        self.ta = tas[0]
        self.ta_var = tas[1]

        # And save the non control parameter values for cfunc calls
        self.pars_no_control = self.ta.pars[4:].tolist()

        # We compute convenient quantities
        self.nseg = len(self.controls) // 4
        self.nseg_fwd = int(self.nseg * cut)
        self.nseg_bck = self.nseg - self.nseg_fwd

        # Sanity checks
        # On the tas
        if len(self.ta.state) != 7:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator state dimension of {len(self.ta.state)}, while 7 is required"
            )
        if len(self.ta.pars) <= 4:
            raise ValueError(
                f"Attempting to construct a zoh_leg with a Taylor Adaptive integrator parameters dimension of {len(self.ta.pars)}, while >=4 is required"
            )
        if self.ta_var: # we skip these lines if no variational integrator is provided
            if len(self.ta_var.state) != 7 + 7 * 7 + 7 * 4:
                raise ValueError(
                    f"Attempting to construct a zoh_leg with a variational Taylor Adaptive integrator state dimension of {len(self.ta_var.state)}, while 84 is required"
                )
            if len(self.ta_var.pars) != len(self.ta.pars):
                raise ValueError(
                    f"While constructing a zoh_leg, the number of parameters in the variational version of the Taylor integrator and the non variational version were detected as different, while its required they are equal"
                )
            # this assumes state of 7 and 4 pars
            self.ic_var = (_np.hstack((_np.eye(7, 7), _np.zeros((7, 4))))).flatten()

        # On the rest
        if len(self.controls) % 4 > 0:
            raise ValueError(
                "In a zoh_leg controls must be multiple of 4: [T, ix, iy, iz] * nseg"
            )
        if len(tgrid) != self.nseg + 1:
            raise ValueError(
                "The t_grid and the controls have incompatible lenghts. It must be nseg*4 and nseg+1"
            )

        # We compile the function for the dynamics (this is used in the gradient computations)
        sys = self.ta.sys
        vars = [it[0] for it in sys]
        dyn = [it[1] for it in sys]
        self.dyn_cfunc = _hy.cfunc_dbl(dyn, vars, compact_mode=True)

    def compute_mismatch_constraints(self):
        """Propagates forward from *state0* and backward from *state1* and returns
        the 7-component state mismatch at the midpoint.

        Returns:
            :class:`list`: Mismatch vector of length 7. All entries are zero for a feasible transfer.
        """
        # Forward segments
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            # setting T, ix, iy, iz
            self.ta.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            if self.max_steps is not None:
                self.ta.propagate_until(self.tgrid[i + 1], max_steps = self.max_steps)
            else:
                self.ta.propagate_until(self.tgrid[i + 1])
        state_fwd = self.ta.state.copy()

        # Backward segments
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            # setting T, ix, iy, iz
            self.ta.pars[0] = self.controls[-4 * i - 4]
            self.ta.pars[1] = self.controls[-4 * i - 3]
            self.ta.pars[2] = self.controls[-4 * i - 2]
            self.ta.pars[3] = self.controls[-4 * i - 1]
            # propagating
            if self.max_steps is not None:
                self.ta.propagate_until(self.tgrid[-2 - i], max_steps = self.max_steps)
            else:   
                self.ta.propagate_until(self.tgrid[-2 - i])
        state_bck = self.ta.state
        return (state_fwd - state_bck).tolist()

    def compute_throttle_constraints(self):
        """Computes the throttle unit-norm constraints :math:`i_x^2 + i_y^2 + i_z^2 - 1` for every segment.

        Returns:
            :class:`list`: Constraint values of length nseg. All entries are zero when the direction
            vector is unit-norm on every segment.
        """
        retval = [0] * self.nseg
        for i in range(self.nseg):
            retval[i] = (
                self.controls[1 + 4 * i] ** 2
                + self.controls[2 + 4 * i] ** 2
                + self.controls[3 + 4 * i] ** 2
                - 1
            )
        return retval

    def compute_tc_grad(self):
        """Computes the gradients of the throttles constraints. Introducing the control vector as
        :math:`\\mathbf u = [T_0, i_{x0}, i_{y0}, i_{z0}, T_1, i_{x1}, i_{y1}, i_{z1}, ...]`, this method computes the following gradient:

        .. math::
        
           \\frac{\\partial \\mathbf {tc}}{\\partial \\mathbf u} \\rightarrow (\\mathbf{nseg} \\times 4\\mathbf{nseg}) 

        Returns:
            :class:`tuple` [:class:`numpy.ndarray`]: The gradient. Size will be (nseg,4nseg).
        """
        retval = _np.zeros((self.nseg, self.nseg * 4))
        for i in range(self.nseg):
            retval[i, 4 * i + 1] = 2 * self.controls[4 * i + 1]
            retval[i, 4 * i + 2] = 2 * self.controls[4 * i + 2]
            retval[i, 4 * i + 3] = 2 * self.controls[4 * i + 3]
        return retval

    def compute_mc_grad(self):
        """
        Computes the gradients of the mismatch constraints. Indicating the initial augmented state with :math:`\\mathbf x_s = [\\mathbf r_s, \\mathbf v_s, m_s]`, the
        final augmented state with :math:`\\mathbf x_f = [\\mathbf r_f, \\mathbf v_f, m_f]`, the time grid as :math:`T_{grid}` and the introducing the control vector
        :math:`\\mathbf u = [T_0, i_{x0}, i_{y0}, i_{z0}, T_1, i_{x1}, i_{y1}, i_{z1}, \\ldots]`, this method computes the following gradients:

        .. math::
        
           \\frac{\\partial \\mathbf {mc}}{\\partial \\mathbf x_s}  \\rightarrow (7\\times7)

        .. math::
        
           \\frac{\\partial \\mathbf {mc}}{\\partial \\mathbf x_f}  \\rightarrow (7\\times7)

        .. math::
        
           \\frac{\\partial \\mathbf {mc}}{\\partial \\mathbf u}  \\rightarrow (7\\times(4\\mathbf{nseg}))
        
        .. math::
        
           \\frac{\\partial \\mathbf {mc}}{\\partial \\mathbf T_{grid}}  \\rightarrow (7\\times(\\mathbf{nseg} + 1))

        Returns:
            :class:`tuple` [:class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`, :class:`numpy.ndarray`]: The four gradients. sizes will be (7,7), (7,7) (7,4nseg) and (7,nseg+1)
  """
        # We will return 4 matrices: dmc/dx0, dmc/dx1, dmc/dcontrol, dmc/dtgrid
        # STMs -> forward
        M_seg_fwd = []  # M10, M21, M32, ....
        M_fwd = []  # Mf0, Mf1, Mf2, ...
        C_seg_fwd = []  # dx_(i+1)/dcontrols_i
        C_fwd = _np.zeros((7, 4 * self.nseg_fwd))  # dx_f/dcontrols_i
        dyn_fwd = []
        self.ta_var.time = self.tgrid[0]
        self.ta_var.state[:7] = self.state0
        for i in range(self.nseg_fwd):
            self.ta_var.state[7:] = self.ic_var
            # setting T, ix, iy, iz
            self.ta_var.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            if self.max_steps is not None:
                self.ta_var.propagate_until(self.tgrid[i + 1], max_steps=self.max_steps)
            else:
                self.ta_var.propagate_until(self.tgrid[i + 1])
            # extracting the segment STMs
            M_seg_fwd.append(self.ta_var.state[7:].reshape(7, 11)[:, :7].copy())  # 7x7
            # extracting the control sensitivities in across the single segment
            C_seg_fwd.append(self.ta_var.state[7:].reshape(7, 11)[:, 7:].copy())  # 7x4
            # computing the dynamics for future usage
            dyn_fwd.append(
                self.dyn_cfunc(
                    self.ta_var.state[:7],
                    pars=[*self.controls[4 * i : 4 * i + 4], *self.pars_no_control],
                )
            )
        # We compute the STMs - Mf0, Mf1, Mf2, ...
        cur = _np.eye(7)
        for M in reversed(M_seg_fwd):
            cur = cur @ M
            M_fwd.append(cur)
        M_fwd = list(reversed(M_fwd)) + [_np.eye(7)]

        # 1 - dmc/dx0
        dmc_dx0 = M_fwd[0]

        # 2 - dmc/dcontrols
        i = 0
        for M, C in zip(M_fwd[1:], C_seg_fwd):
            C_fwd[:, 4 * i : 4 * i + 4] = M @ C
            i += 1

        # 3 - dmc/dtgrid
        dmcdtgrid = _np.zeros((7, self.nseg + 1))
        if self.nseg_fwd > 0:
            dmcdtgrid[:, 0] = -M_fwd[1] @ dyn_fwd[0]
            dmcdtgrid[:, self.nseg_fwd] = M_fwd[-1] @ dyn_fwd[-1]

            for i in range(1, self.nseg_fwd):
                dmcdtgrid[:, i] = M_fwd[i + 1] @ (
                    M_seg_fwd[i] @ dyn_fwd[i - 1] - dyn_fwd[i]
                )

        # STMs -> backward
        M_seg_bck = []  # M10, M21, M32, .... note: numbering reversed
        M_bck = []  # Mf0, Mf1, Mf2, .... note: numbering reversed
        C_seg_bck = []  # dxi/dcontrols
        C_bck = _np.zeros((7, 4 * self.nseg_bck))
        dyn_bck = []
        self.ta_var.time = self.tgrid[-1]
        self.ta_var.state[:7] = self.state1
        for i in range(self.nseg_bck):
            self.ta_var.state[7:] = self.ic_var
            # setting T, ix, iy, iz
            self.ta_var.pars[0] = self.controls[-4 * i - 4]
            self.ta_var.pars[1] = self.controls[-4 * i - 3]
            self.ta_var.pars[2] = self.controls[-4 * i - 2]
            self.ta_var.pars[3] = self.controls[-4 * i - 1]
            # propagating
            if self.max_steps is not None:
                self.ta_var.propagate_until(self.tgrid[-2 - i], max_steps=self.max_steps)
            else:
                self.ta_var.propagate_until(self.tgrid[-2 - i])
            # extracting the segment STMs
            M_seg_bck.append(self.ta_var.state[7:].reshape(7, 11)[:, :7].copy())
            # extracting the control sensitivities in across the single segment
            C_seg_bck.append(self.ta_var.state[7:].reshape(7, 11)[:, 7:].copy())  # 7x4
            # computing the dynamics for future usage
            dyn_bck.append(
                self.dyn_cfunc(
                    self.ta_var.state[:7],
                    pars=[self.controls[-4 * i - 4], self.controls[-4 * i - 3],
                          self.controls[-4 * i - 2], self.controls[-4 * i - 1]]
                    + self.pars_no_control,
                )
            )

        # We compute the STMs - Mf0, Mf1, Mf2, ...
        cur = _np.eye(7)
        for M in reversed(M_seg_bck):
            cur = cur @ M
            M_bck.append(cur)
        M_bck = list(reversed(M_bck)) + [_np.eye(7)]

        # 1 - dmc/dxf
        dmc_dx1 = -M_bck[0]  # mc = xfwd-x_bck hence the minus

        # 2 - dmc/dcontrols
        # NOTE: the controls are reversed in order and we need to account for it
        i = 0
        for M, C in zip(M_bck[1:], C_seg_bck):
            C_bck[:, 4 * self.nseg_bck - 4 - 4 * i : 4 * self.nseg_bck - 4 * i] = M @ C
            i += 1

        # 3 - dmc/dtgrid
        if (
            self.nseg_bck > 0
        ):  # we skip this in the corner case where no bck seg are there
            dmcdtgrid[:, -1] = M_bck[1] @ dyn_bck[0]
            dmcdtgrid[:, self.nseg_fwd] -= (
                M_bck[-1] @ dyn_bck[-1]
            )  # This is for the mid time point shared fwd and bck: we need to subtract to the existing

            for i in range(1, self.nseg_bck):
                dmcdtgrid[:, -1 - i] = -M_bck[i + 1] @ (
                    M_seg_bck[i] @ dyn_bck[i - 1] - dyn_bck[i]
                )

        # dmc/dcontrols (7x4*nseg)
        dmc_dcontrols = _np.hstack((C_fwd, -C_bck))  # mc = xfwd-x_bck hence the minus

        return dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid

    def get_state_info(self, N=2):
        """
        This method returns state histories sampled along each ZOH segment, for both the
        forward and backward propagation parts of the leg. The sampling is performed by
        calling :meth:`heyoka.taylor_adaptive.propagate_grid` on a uniformly-spaced grid
        of *N* points within each segment.

        Args:
            *N* (:class:`int`): Number of sampling points per segment (including the segment
            endpoints). The default (*N=2*) returns only the segment endpoints.

        Returns:
            :class:`tuple`: ``(state_fwd, state_bck)`` where:

                - ``state_fwd`` (:class:`list`): List of length ``nseg_fwd``. Each entry contains
                  the sampled 7D state history over the corresponding forward segment (from
                  ``tgrid[i]`` to ``tgrid[i+1]``).

                - ``state_bck`` (:class:`list`): List of length ``nseg_bck``. Each entry contains
                  the sampled 7D state history over the corresponding backward segment (from
                  ``tgrid[-1-i]`` to ``tgrid[-2-i]``).

            .. note::
               The backward propagation is carried out by integrating from the final time toward
               earlier times; depending on the backend conventions, the returned grids may thus be
               time-reversed with respect to a forward-time plot.

            .. note::
               This method uses the nominal integrator ``self.ta`` and overwrites its internal
               ``time``, ``state`` and the first four parameters (``T, i_x, i_y, i_z``). If the
               integrator state must be preserved, call this method on a dedicated copy.

        Examples:

        .. code-block:: python 
         
           ax = pk.plot.make_3Daxis()
           fwd, bck = leg.get_state_info(N=100)
           for segment in fwd:
               ax.scatter(segment[0,0], segment[0,1], segment[0,2], c='blue')
               ax.plot(segment[:,0], segment[:,1], segment[:,2], c='blue')
           for segment in bck:
               ax.scatter(segment[0,0], segment[0,1], segment[0,2], c='darkorange')
               ax.plot(segment[:,0], segment[:,1], segment[:,2], c='darkorange')
           ax.view_init(90, -90)
        """
        # Forward segments
        state_fwd = []
        self.ta.time = self.tgrid[0]
        self.ta.state[:] = self.state0
        for i in range(self.nseg_fwd):
            # setting T, ix, iy, iz
            self.ta.pars[0:4] = self.controls[4 * i : 4 * i + 4]
            # propagating
            plot_grid_fwd = _np.linspace(self.tgrid[i], self.tgrid[i + 1], N)
            if self.max_steps is not None:
                sol_fwd = self.ta.propagate_grid(plot_grid_fwd, max_steps=self.max_steps)[-1]
            else:
                sol_fwd = self.ta.propagate_grid(plot_grid_fwd)[-1]
            state_fwd.append(sol_fwd)
        
        # Backward segments
        state_bck = []
        self.ta.time = self.tgrid[-1]
        self.ta.state[:] = self.state1
        for i in range(self.nseg_bck):
            # setting T, ix, iy, iz
            self.ta.pars[0] = self.controls[-4 * i - 4]
            self.ta.pars[1] = self.controls[-4 * i - 3]
            self.ta.pars[2] = self.controls[-4 * i - 2]
            self.ta.pars[3] = self.controls[-4 * i - 1]
            # propagating
            plot_grid_bck = _np.linspace(self.tgrid[-1 - i], self.tgrid[-2 - i], N)
            if self.max_steps is not None:
                sol_bck = self.ta.propagate_grid(plot_grid_bck, max_steps=self.max_steps)[-1]
            else:
                sol_bck = self.ta.propagate_grid(plot_grid_bck)[-1]
            state_bck.append(sol_bck)
        return state_fwd,state_bck
