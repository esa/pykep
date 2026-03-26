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

    This problem works internally using the :class:`~pykep.leg.sims_flanagan` and manipulates its transfer time T,
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
        states=[1.2, 0.0, -0.01, 0.01, 1.0, -0.01],
        statef=[1.0, 0.0, -0.0, 0.01, 1.1, -0.0],
        ms=1.0,
        max_thrust=0.22,
        tof_bounds=[3.4, 8.6],
        mf_bounds=[0.2, 1],
        nseg=10,
        cut=0.6,
        tas=(_pk.ta.get_zoh_kep(1e-10), None),
        time_encoding="uniform",
        w_bounds_softmax=[-1.0, 1.0],
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

            *cut* (:class:`float`): Cut parameter for the :class:`~pykep.leg.sims_flanagan`. Defaults to 0.6.

            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators

                - `ta`: Nominal dynamics (state dim 7, pars ≥ 4)

                - `ta_var`: Variational dynamics (state dim 84, same pars). When None, no gradients will be used.

        """
        # We define some additional datamembers useful later-on
        self.nseg = nseg
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.max_thrust = max_thrust
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = w_bounds_softmax

        supported_time_encodings = [
            "uniform",
            "softmax",
        ]  # we will add variable length here in the future
        if self.time_encoding not in supported_time_encodings:
            raise NotImplementedError(
                f"Only {supported_time_encodings} time encodings are currently implemented"
            )

        # We build and store a ZOH leg as data member.
        # We will change controls, tgrid and mf as those are encoded in the decision vector,
        # but to construct we need some values to construct the first instance...

        controls = _np.random.uniform(-1, 1, (4 * nseg,))
        controls[0::4] = _np.abs(controls[0::4])  # force will be in [0, 1]
        controls[0::4] *= self.max_thrust  # force will be in [0, max_thrust]
        self.controls = controls
        tgrid = _np.linspace(
            0, (self.tof_bounds[0] + self.tof_bounds[1]) / 2, self.nseg + 1
        )
        self.leg = _pk.leg.zoh(
            states + [ms], list(controls), statef + [ms], tgrid, cut, tas
        )

    def _set_leg_from_x(self, x):
        # Here is where the decision vector gets decoded into the leg tgrid, mf and controls
        # Lets do this case: x = [mf] + controls + tof
        self.leg.state1[-1] = x[0]
        self.leg.controls = x[1 : 1 + 4 * self.nseg].copy()
        self.leg.controls[0::4] *= self.max_thrust
        self.leg.controls = list(self.leg.controls)
        tof = x[1 + 4 * self.nseg]
        # Since we only have tof in the decision vector we assume a uniform epoch grid
        if self.time_encoding == "uniform":
            self.leg.tgrid = _np.linspace(0, tof, self.nseg + 1)
        elif self.time_encoding == "softmax":
            w = x[2 + 4 * self.nseg : 2 + 4 * self.nseg + self.nseg]
            softmax_weights, _ = _pk.compute_softmax_and_jacobian(w)
            segment_duration = tof * softmax_weights
            # now we need to transform the segment durations into a time grid increasing monotically from 0 to tof
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

        # We optimize for maximum final mass (minimum propellent)
        obj = -x[0]

        # We compute the equality constraints
        ceq = self.leg.compute_mismatch_constraints()
        ceq += self.leg.compute_throttle_constraints()
        retval = _np.array([obj] + ceq)
        return retval

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

        where controls = [T_1, ix_1, iy_1, iz_1, ..., T_nseg, ix_nseg, iy_nseg, iz_nseg]
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
            gradient[1:8, 1 + 4 * i : 1 + 4 * i + 4] = dmc_dcontrols[
                :, 4 * i : 4 * i + 4
            ]
            gradient[
                1:8, 1 + 4 * i
            ] *= (
                self.max_thrust
            )  # account for the multiplication by max_thrust for the T component of the controls
        ## Third contribution -> partials of mismatch constraints w.r.t. time grid
        if self.time_encoding == "uniform":
            # note here that what we are after is:
            # dmc/dtof = dmc/dtgrid * dtgrid/dtof
            # and since for the uniform case tgrid = linspace(0,tof,nseg+1) = [0, tof/nseg, 2*tof/nseg, ..., tof], we have that dtgrid/dtof = [0, 1/nseg, 2/nseg, ..., 1]
            gradient[1:8, -1] = dmcdtgrid @ _np.linspace(
                0, 1, self.nseg + 1
            )  # dmc/dtof = dmc/dtgrid * dtgrid/dtof
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
            tril_matrix = _np.zeros((self.nseg + 1, self.nseg))
            for i in range(1, self.nseg + 1):
                tril_matrix[i, :i] = 1.0
            dtgrid_dw = tof * (tril_matrix @ J_softmax)
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
        dtc_dcontrols = self.leg.compute_tc_grad()
        gradient[8 : 8 + self.nseg, 1 : 1 + 4 * self.nseg] = (
            dtc_dcontrols  # the partials of the throttle constraints w.r.t. the controls of all segments
        )
        ## Third contribution -> partials of throttle constraints w.r.t. controls
        # the throttle constrain does not have dependence on the time grid, so these are zeros, nothing to do here
        ## Fourth contribution (softmax only) -> partials of throttle constraints w.r.t. softmax weights: they are zero so nothing to do here
        return gradient.flatten()

    # If the variational integator is also passed in construction (i.e. its not None)
    # gradient should be provided (pagmo UDP interface used this method to know)
    def has_gradient(self):
        return self.with_gradient

    # Only equality contraints. Mismatches and the |i_u| = 1.
    def get_nec(self):
        return 7 + self.nseg

    def pretty(self, x):
        """
        Prints a detailed representation of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containins: final mass, thrust direction, time of flight and (if time encoding is softmax) the weights for the softmax time grid.
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
        **kwargs,
    ):
        """
        Plots the trajectory of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containins: final mass, thrust direction, time of flight and (if time encoding is softmax) the weights for the softmax time grid.
            
            *ax* (:class:`matplotlib.axes.Axes`): The matplotlib axes to plot on. If None, a new figure and axes will be created.
            
            *N* (:class:`int`): The number of points to plot along the trajectory.
            
            *to_cartesian* (:class:`~collections.abc.Callable`): A function that converts whatever state is used in the internal Taylor integrator to Cartesian (r,v).
            
            *mark_segments* (:class:`bool`): adds markers ath each segment edge

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the Lambert's problem trajectory added.
        """
        x_arr = _np.asarray(x)
        # to be replaced with a plot method akin the sims-flanagan point2point one
        self._set_leg_from_x(x_arr)
        fwd, bck = self.leg.get_state_info(N=N)
        # compute the color scheme
        throttles = x_arr[1 : 1 + 4 * self.nseg : 4]
        if ax is None:
            ax = _pk.plot.make_3Daxis()
        # plot
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
        if mark_mismatch:
            ax.scatter(
                segment_cart[-1, 0],
                segment_cart[-1, 1],
                segment_cart[-1, 2],
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
        if mark_mismatch:
            ax.scatter(
                segment_cart[-1, 0],
                segment_cart[-1, 1],
                segment_cart[-1, 2],
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

