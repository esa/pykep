## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
##
## This file is part of the kep3 library.##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as _np
import pykep as _pk
from matplotlib import pyplot as _plt


class zoh_point2point:
    """Represents the optimal low-thrust transfer between two fixed points using `pykep`'s Zero Order Hold (direct) trajectory legs.

    This problem works internally using the :class:`~pykep.leg.zoh` and manipulates its transfer time T,
    final mass mf and the controls as to link two fixed points in space with a low-thrust trajectory.

    It can be used to better profile and understand performances of optimizers on this type of direct approach, but has a limited use
    in the design of interplanetary trajectories as per the fixed point limitation.

    The decision vector is::

        x = [mf] + controls + tof (+ [weights for softmax])

    where controls is a vector of control parameters :math:`[T, i_x, i_y, i_z] \\times n_\\text{seg}` representing magnitude and direction
    of the thrust applied in each segment.
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
        inequalities_for_tc = False,
        max_steps=None,
    ):
        """
        Initializes the zoh_point2point instance with given parameters.

        Args:
            *states* (:class:`list`): Initial state (only the first 6 states, i.e. no mass). Units as expected by the numerical integrator, defaults to [1.2, 0.0, -0.01, 0.01, 1.0, -0.01].

            *statef* (:class:`list`): Final state (only the first 6 states, i.e. no mass). Units as expected by the numerical integrator, defaults to [1.0, 0.0, -0.0, 0.01, 1.1, -0.0].

            *ms* (:class:`float`): Initial mass. Units as expected by the numerical integrator, defaults to 1.

            *max_thrust* (:class:`float`): Maximum thrust. Units as expected by the numerical integrator, defaults to 0.22.

            *tof_bounds* (:class:`list`): Bounds for time of flight. Units as expected by the numerical integrator, defaults to [3.4, 8.6].

            *mf_bounds* (:class:`list`): Bounds for final mass. Units as expected by the numerical integrator, defaults to [0.2, 1].

            *nseg* (:class:`int`): Number of segments for the trajectory. Defaults to 10.

            *cut* (:class:`float`): Cut parameter for the :class:`~pykep.leg.zoh`. Defaults to 0.6.

            *w_bounds_softmax* (:class:`list`): Bounds for the softmax weights (only used if time_encoding is 'softmax'). Defaults to [-1.0, 1.0].
            
            *inequalities_for_tc* (:class:`bool`): If True, the throttle constraints are formulated as inequalities (|i_u|**2 - 1 <= 0), otherwise they are formulated as equalities (|i_u| = 1). Defaults to False (equality formulation).

            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators

                - `ta`: Nominal dynamics (state dim 7, pars ≥ 4)

                - `ta_var`: Variational dynamics (state dim 84, same pars). When None, no gradients will be used.

            *max_steps* (:class:`int` or None): Maximum number of Taylor integrator steps per propagation call. When None, uses the default integrator behavior.

        """
        # Initialize defaults for optional constructor arguments.
        if states is None:
            states = [1.2, 0.0, -0.01, 0.01, 1.0, -0.01]
        if statef is None:
            statef = [1.0, 0.0, -0.0, 0.01, 1.1, -0.0]
        if tof_bounds is None:
            tof_bounds = [3.4, 8.6]
        if mf_bounds is None:
            mf_bounds = [0.2, 1]
        if w_bounds_softmax is None:
            w_bounds_softmax = [-1.0, 1.0]

        self.nseg = nseg
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.max_thrust = max_thrust
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = w_bounds_softmax
        self.max_steps = max_steps
        # Precompute nseg-dependent helpers reused in gradient() to avoid
        # rebuilding identical arrays at every optimizer call.
        self._dtgrid_dtof_uniform = _np.linspace(0.0, 1.0, self.nseg + 1)
        self._softmax_cumsum_matrix = _np.vstack(
            [_np.zeros(self.nseg), _np.tril(_np.ones((self.nseg, self.nseg)))]
        )
        self.inequalities_for_tc = inequalities_for_tc

        supported_time_encodings = [
            "uniform",
            "softmax",
        ]  # Variable-length encoding can be added here in the future.
        if self.time_encoding not in supported_time_encodings:
            raise NotImplementedError(
                f"Only {supported_time_encodings} time encodings are currently implemented"
            )

        # Build and store a ZOH leg instance.
        # Controls, time grid, and final mass are overwritten from the decision vector.
        controls = _np.random.uniform(-1, 1, (4 * nseg,))
        controls[0::4] = _np.abs(controls[0::4])  # force will be in [0, 1]
        controls[0::4] *= self.max_thrust  # force will be in [0, max_thrust]
        tgrid = _np.linspace(
            0, (self.tof_bounds[0] + self.tof_bounds[1]) / 2, self.nseg + 1
        )
        self.leg = _pk.leg.zoh(
            states + [ms],
            controls.tolist(),
            statef + [ms],
            tgrid,
            cut,
            tas,
            max_steps=self.max_steps,
        )

    def _compute_throttle_constraints(self):
        """Computes throttle constraints from the leg controls."""
        tc = _np.zeros(self.nseg)
        controls = self.leg.controls
        for i in range(self.nseg):
            base = 4 * i
            tc[i] = (
                controls[base + 1] ** 2
                + controls[base + 2] ** 2
                + controls[base + 3] ** 2
                - 1.0
            )
        return tc.tolist()

    def compute_tc_grad(self):
        """Computes throttle Jacobian w.r.t. leg controls."""
        dtc_dcontrols = _np.zeros((self.nseg, 4 * self.nseg))
        controls = self.leg.controls
        for i in range(self.nseg):
            base = 4 * i
            ix = controls[base + 1]
            iy = controls[base + 2]
            iz = controls[base + 3]
            dtc_dcontrols[i, base + 1] = 2.0 * ix
            dtc_dcontrols[i, base + 2] = 2.0 * iy
            dtc_dcontrols[i, base + 3] = 2.0 * iz
        return dtc_dcontrols

    def _set_leg_from_x(self, x):
        # Decode the decision vector into the leg tgrid, mf, and controls.
        # Here: x = [mf] + controls + tof
        state1 = list(self.leg.state1)
        state1[-1] = x[0]
        self.leg.state1 = state1
        controls = list(x[1 : 1 + 4 * self.nseg])
        controls[0::4] = [v * self.max_thrust for v in controls[0::4]]
        self.leg.controls = controls
        tof = x[1 + 4 * self.nseg]
        # Since we only have tof in the decision vector we assume a uniform epoch grid
        if self.time_encoding == "uniform":
            self.leg.tgrid = _np.linspace(0, tof, self.nseg + 1)
        elif self.time_encoding == "softmax":
            w = x[2 + 4 * self.nseg : 2 + 4 * self.nseg + self.nseg]
            softmax_weights, _ = _pk.compute_softmax_and_jacobian(w)
            segment_duration = tof * softmax_weights
            # Transform segment durations into a monotonically increasing time grid from 0 to tof.
            tgrid = _np.zeros(self.nseg + 1)
            tgrid[1:] = _np.cumsum(segment_duration)
            self.leg.tgrid = tgrid

    def get_bounds(self):
        if self.time_encoding == "uniform":
            lb = (
                [self.mf_bounds[0]]
                + [0, -1.0, -1.0, -1.0] * self.nseg
                + [self.tof_bounds[0]]
            )
            ub = (
                [self.mf_bounds[1]]
                + [1.0, 1.0, 1.0, 1.0] * self.nseg
                + [self.tof_bounds[1]]
            )
        elif self.time_encoding == "softmax":
            lb = (
                [self.mf_bounds[0]]
                + [0, -1.0, -1.0, -1.0] * self.nseg
                + [self.tof_bounds[0]]
                + [self.w_bounds_softmax[0]] * self.nseg
            )
            ub = (
                [self.mf_bounds[1]]
                + [1.0, 1.0, 1.0, 1.0] * self.nseg
                + [self.tof_bounds[1]]
                + [self.w_bounds_softmax[1]] * self.nseg
            )
        return (lb, ub)

    def fitness(self, x):
        # We set the leg using data in the decision vector
        self._set_leg_from_x(x)

        # Optimize for maximum final mass (minimum propellant).
        obj = -x[0]

        # We compute all constraints. The tc will be treated as equality or inequality
        # according to the value of self.inequalities_for_tc.
        cmc = self.leg.compute_mismatch_constraints()
        ctc = self._compute_throttle_constraints()
        
        return [obj] + cmc + ctc

    def gradient(self, x):
        """
        Computes the gradient of the fitness function with respect to the decision vector.

        The fitness function is:
            f(x) = [obj] + mismatch_constraints + throttle_constraints

        where:
            - obj = -mf (maximize final mass)
            - mismatch_constraints: 7 equality constraints (position, velocity, mass mismatch)
            - throttle_constraints: nseg equality constraints (|u|^2 - 1 = 0 for each segment)

        The decision vector is:
            x = [mf] + controls + tof  (+ weights for softmax)

        where controls are:
            - [T_1, ix_1, iy_1, iz_1, ..., T_nseg, ix_nseg, iy_nseg, iz_nseg]
        and weights = [w_1, w_2, ..., w_nseg] (only if time_encoding is 'softmax')

        Returns:
            gradient: numpy array of shape (nf, nx) where:
                nf = 1 + 7 + nseg (total fitness dimension)
                nx = 1 + 4*nseg + 1 (+ nseg for softmax) (decision vector dimension)
        """
        if not self.with_gradient:
            raise RuntimeError(
                "Gradient computation requires variational integrator (tas[1] must not be None)"
            )
        # let's start by setting the leg from the decision vector
        self._set_leg_from_x(x)
        # now we initialize the gradient matrix:
        nf = 1 + 7 + self.nseg  # total number of fitness components
        nx = 1 + 4 * self.nseg + 1  # total number of decision variables
        if self.time_encoding == "softmax":
            nx += (
                self.nseg
            )  # account for the additional decision variables for the softmax weights
        gradient = _np.zeros((nf, nx))

        # ============================================
        # Gradient of objective w.r.t decision vector:
        # ============================================
        gradient[0, 0] = -1.0  # d(-mf)/dmf = -1 .. all the rest is zeros

        # ============================================
        # Gradient of mismatch constraints w.r.t decision vector, made of three contributions:
        # ============================================
        ## First contribution -> partials of mismatch constraints w.r.t. final mass
        dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid = self.leg.compute_mc_grad()
        gradient[1:8, 0] = dmc_dx1[:, 6]  # mismatch constraints w.r.t mf
        ## Second contribution -> partials of mismatch constraints w.r.t. controls
        for i in range(self.nseg):
            dmc_block = dmc_dcontrols[:, 4 * i : 4 * i + 4]
            gradient[1:8, 1 + 4 * i : 1 + 4 * i + 4] = dmc_block
            gradient[1:8, 1 + 4 * i] *= self.max_thrust
        ## Third contribution -> partials of mismatch constraints w.r.t. time grid
        if self.time_encoding == "uniform":
            # note here that what we are after is:
            # dmc/dtof = dmc/dtgrid * dtgrid/dtof
            # and since for the uniform case tgrid = linspace(0,tof,nseg+1) = [0, tof/nseg, 2*tof/nseg, ..., tof], we have that dtgrid/dtof = [0, 1/nseg, 2/nseg, ..., 1]
            gradient[1:8, -1] = dmcdtgrid @ self._dtgrid_dtof_uniform
        elif self.time_encoding == "softmax":
            # let's start by extracting the softmax weights from the decision vector:
            w = x[2 + 4 * self.nseg : 2 + 4 * self.nseg + self.nseg]
            softmax_weights, J_softmax = _pk.compute_softmax_and_jacobian(w)
            tof = x[1 + 4 * self.nseg]

            # gradient w.r.t. tof:
            # tgrid = cumsum(tof*softmax_weights) => dtgrid/dtof = cumsum(softmax_weights)
            # note that also tgrid[0]=0 so dtgrid[0]/dtof = 0
            dtgrid_dtof = _np.zeros(self.nseg + 1)
            dtgrid_dtof[1:] = _np.cumsum(softmax_weights)
            gradient[1:8, 1 + 4 * self.nseg] = dmcdtgrid @ dtgrid_dtof

            # now the gradient w.r.t. softmax weights:
            # tgrid[i]=sum_{j=0}^{i-1}(tof*softmax_weights[j]) =>
            # dtgrid[i]/dw[k] = sum_{j=0}^{i-1} tof*dsoftmax_weights[j]/dw[k] =
            #                = tof*sum_{j=0}^{i-1} J_softmax[j,k]
            # now if we assume that: tril_matrix[i,j]=1 if j<i, else 0
            # then dtgrid/dw = tof*(tril_matrix @ J_softmax)
            dtgrid_dw = tof * (self._softmax_cumsum_matrix @ J_softmax)
            gradient[1:8, 2 + 4 * self.nseg : 2 + 4 * self.nseg + self.nseg] = (
                dmcdtgrid @ dtgrid_dw
            )
        else:
            raise NotImplementedError(
                f"Gradient w.r.t. tof not implemented for {self.time_encoding} time encoding"
            )
        # ============================================
        # Gradient of throttle constraints w.r.t. decision vector, made of two contributions:
        # ============================================
        ## First contribution -> partials of throttle constraints w.r.t. final mass
        gradient[8 : 8 + self.nseg, 0] = 0.0  # these are zeros
        ## Second contribution -> partials of throttle constraints w.r.t. controls
        dtc_dcontrols = self.compute_tc_grad()
        for i in range(self.nseg):
            dtc_block = dtc_dcontrols[i : i + 1, 4 * i : 4 * i + 4]
            gradient[
                8 + i : 8 + i + 1, 1 + 4 * i : 1 + 4 * i + 4
            ] = dtc_block
            gradient[8 + i : 8 + i + 1, 1 + 4 * i] *= self.max_thrust
        ## Third contribution -> partials of throttle constraints w.r.t. time grid
        # Throttle constraints are independent of the time grid, so this block is zero.
        ## Fourth contribution (softmax only) -> partials of throttle constraints w.r.t. softmax weights
        # This block is also zero because throttle constraints do not depend on the time encoding.
        return gradient.flatten()

    # If the variational integrator is passed at construction time (i.e., it is not None),
    # `gradient` should be provided (the pagmo UDP interface uses this method to check).
    def has_gradient(self):
        return self.with_gradient

    # Equality constraints only: mismatches and |i_u| = 1.
    def get_nec(self):
        if self.inequalities_for_tc:
            # We formulate the throttle constraints as inequalities:  |i_u|**2 - sss1 <= 0
            return 7
        else:
            return 7 + self.nseg
    
    def get_nic(self):
        if self.inequalities_for_tc:
            # We formulate the throttle constraints as inequalities:  |i_u|**2 -1 <= 0
            return self.nseg
        else:
            return 0

    def pretty(self, x):
        """
        Prints a detailed representation of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, time of flight, and (if time encoding is softmax) the softmax weights.
        """
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
        orbit_color = None,
        tof=None,
        **kwargs,
    ):
        """
        Plots the trajectory of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, time of flight, and (if time encoding is softmax) the softmax weights.
            
            *ax* (:class:`matplotlib.axes.Axes`): The matplotlib axes to plot on. If None, a new figure and axes will be created.
            
            *N* (:class:`int`): The number of points to plot along the trajectory.
            
            *to_cartesian* (:class:`~collections.abc.Callable`): A function that converts whatever state is used in the internal Taylor integrator to Cartesian (r,v).
            
            *mark_segments* (:class:`bool`): adds markers at each segment edge

            *mark_mismatch* (:class:`bool`): highlights mismatch constraints

            *orbit_color* (:class:`str`): when not None is the color of the border trajectories propagated for the same length as the trajectory encoded by x.
       
            *tof* (:class:`float`): time of flight for the border trajectories (only used if orbit_color is not None). If None, uses the tof encoded in x.
       
            **kwargs* (:class:`dict`): Additional keyword arguments to pass to the plotting functions (e.g., marker size, color, etc.)
       
        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the Lambert's problem trajectory added.
        """
        x_arr = _np.asarray(x)

        # We start with the boundaries ....
        if orbit_color is not None:
            if tof is None:
                tof = x_arr[1+4*self.nseg]

            # Start orbit
            ta = self.leg.ta
            ta.state[:] = self.leg.state0
            ta.pars[0] = 0.
            ta.time=0
            sol_s = ta.propagate_grid(_np.linspace(0, tof, N * self.nseg))[-1]
            sol_cart_s = _np.array([to_cartesian(it) for it in sol_s])
            ax.plot(sol_cart_s[:,0], sol_cart_s[:,1], sol_cart_s[:,2], orbit_color, alpha=0.5)

            # End orbit (backpropagated)
            ta.state[:] = self.leg.state1
            ta.time=0
            sol_f = ta.propagate_grid(_np.linspace(0, -tof, N * self.nseg))[-1]
            sol_cart_f = _np.array([to_cartesian(it) for it in sol_f])
            ax.plot(sol_cart_f[:,0], sol_cart_f[:,1], sol_cart_f[:,2], orbit_color, alpha=0.5)

        # And then the trajectory
        self._set_leg_from_x(x_arr)
        fwd, bck, _ = self.leg.get_state_info(N=N)
        # compute the color scheme
        throttles = x_arr[1 : 1 + 4 * self.nseg : 4]
        if ax is None:
            ax = _pk.plot.make_3Daxis()
        # plot
        last_point = None
        for i, segment in enumerate(fwd):
            color = (
                0.25 + (0.80 - 0.25) * throttles[i],
                0.41 + (0.36 - 0.41) * throttles[i],
                0.88 + (0.36 - 0.88) * throttles[i],
            )
            # We obtain the state in Cartesian
            segment_cart = _np.array([to_cartesian(it) for it in segment])
            if mark_segments:
                ax.scatter(
                    segment_cart[0, 0],
                    segment_cart[0, 1],
                    segment_cart[0, 2],
                    **kwargs,
                )
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c=color)
            last_point = segment_cart[-1, :3]
        if mark_mismatch and last_point is not None:
            ax.scatter(
                last_point[0],
                last_point[1],
                last_point[2],
                marker="^",
                **kwargs,
            )
        for i, segment in enumerate(bck):
            color = (
                0.25 + (0.80 - 0.25) * throttles[-1 - i],
                0.41 + (0.36 - 0.41) * throttles[-1 - i],
                0.88 + (0.36 - 0.88) * throttles[-1 - i],
            )
            # We obtain the state in Cartesian
            segment_cart = _np.array([to_cartesian(it) for it in segment])
            if mark_segments:
                ax.scatter(
                    segment_cart[0, 0],
                    segment_cart[0, 1],
                    segment_cart[0, 2],
                    **kwargs,
                )
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c=color)
            last_point = segment_cart[-1, :3]
        if mark_mismatch and last_point is not None:
            ax.scatter(
                last_point[0],
                last_point[1],
                last_point[2],
                marker="^",
                **kwargs
            )


        return ax
    
    def plot_throttle(self, x, ax=None, **kwargs):
        """
        Plots the throttle profile of the zero order hold point to point problem.

        The throttle values are assumed to be defined per segment (length ``nseg``)
        and are plotted as piecewise-constant values over the corresponding time
        intervals defined by ``self.leg.tgrid`` (length ``nseg+1``).

        Args:
            *x* (:class:`list`): The decision vector containing the segment throttles.
                The throttle values are read from ``x[1:4*self.nseg:4]``.
            *ax* (:class:`matplotlib.axes.Axes`): The matplotlib axes to plot on.
                If None, a new figure and axes will be created.

        Returns:
            The matplotlib axes with the throttle profile plotted.
        """
        if ax is None:
            _, ax = _plt.subplots()

        self._set_leg_from_x(x)
        tgrid = _np.asarray(self.leg.tgrid)              # length nseg+1
        throttles = _np.asarray(x[1 : 1 + 4 * self.nseg : 4])  # length nseg

        # Repeat the last throttle so that the last tgrid point is included in the step plot
        u_plot = _np.r_[throttles, throttles[-1]]
        ax.step(tgrid, u_plot, where="post", **kwargs)

        ax.set_xlabel("time grid")
        ax.set_ylabel("throttle value (nd)")
        return ax

