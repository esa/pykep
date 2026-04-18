# Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
#                          Advanced Concepts Team, European Space Agency (ESA)
#
# This file is part of the pykep library.
#
# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as _np
import pykep as _pk

class zoh_ss_point2point:
    """Represents the time optimal Solar Sail transfer between two fixed points using Zero Order Hold (direct) trajectory legs.

    This problem works internally using :class:`~pykep.leg.zoh_py` with
    ``dim_dynamics=6`` and ``dim_controls=2`` and manipulates its transfer time tof
    and the sail clock and cone angles (controls) as to link two fixed points in space with a Solar Sail trajectory.

    It can be used to better profile and understand performances of optimizers on this type of direct approach, but has a limited use
    in the design of interplanetary trajectories as per the fixed point limitation.

    The decision vector is::

        x = controls + [tof] (+ [weights for softmax])

    where controls is a vector of control parameters :math:`[\\alpha, \\beta] \\times n_\\text{seg}` representing the sail cone and clock angles, assumed
    fixed in the RTN frame along each segment.
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
    ):
        """
        Initializes the zoh_ss_point2point instance with given parameters.

        Args:
            *states* (:class:`list`): Initial state. Units as expected by the numerical integrator, defaults to [1.2, 0.0, -0.01, 0.01, 1.0, -0.01].

            *statef* (:class:`list`): Final state. Units as expected by the numerical integrator, defaults to [1.0, 0.0, -0.0, 0.01, 1.1, -0.0].

            *tof_bounds* (:class:`list`): Bounds for time of flight. Units as expected by the numerical integrator, defaults to [3.4, 8.6].

            *nseg* (:class:`int`): Number of segments for the trajectory. Defaults to 10.

            *cut* (:class:`float`): Cut parameter for :class:`~pykep.leg.zoh_py`. Defaults to 0.6.

            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators

                - `ta`: Nominal dynamics (state dim 6, pars ≥ 2)

                - `ta_var`: Variational dynamics (54, same pars). When None, no gradients will be used.

            *max_steps* (:class:`int` or None): Maximum number of Taylor integrator steps per propagation call. When None, uses the default integrator behavior.

        """
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")

        # Initialize defaults for optional constructor arguments.
        if states is None:
            states = [1.2, 0.0, -0.01, 0.01, 1.0, -0.01]
        if statef is None:
            statef = [1.0, 0.0, -0.0, 0.01, 1.1, -0.0]
        if tof_bounds is None:
            tof_bounds = [3.4, 8.6]
        if w_bounds_softmax is None:
            w_bounds_softmax = [-1.0, 1.0]

        self.nseg = nseg
        self.tof_bounds = tof_bounds
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

        supported_time_encodings = [
            "uniform",
            "softmax",
        ]  # Variable-length encoding can be added here in the future.
        if self.time_encoding not in supported_time_encodings:
            raise NotImplementedError(
                f"Only {supported_time_encodings} time encodings are currently implemented"
            )

        # Build and store a ZOH leg instance.
        # The random controls only initialize the leg object once and are
        # overwritten before fitness/gradient/plot evaluations.
        controls = _np.random.uniform(-_np.pi/2, _np.pi/2, (2 * nseg,))
        self.controls = controls
        tgrid = _np.linspace(
            0, (self.tof_bounds[0] + self.tof_bounds[1]) / 2, self.nseg + 1
        )
        self.leg = _pk.leg.zoh_py(
            states,
            list(controls),
            statef,
            tgrid,
            cut,
            tas,
            max_steps=self.max_steps,
            dim_dynamics=6,
            dim_controls=2,
        )

    def _expected_nx(self):
        nx = 2 * self.nseg + 1
        if self.time_encoding == "softmax":
            nx += self.nseg
        return nx

    def _set_leg_from_x(self, x):
        if len(x) != self._expected_nx():
            raise ValueError(
                f"Invalid decision vector length: got {len(x)}, expected {self._expected_nx()}"
            )
        # Decode the decision vector into the leg tgrid and controls.
        # The following is valid for all options of the time encoding assuming x = controls + [tof] + (....)
        controls = list(x[: 2 * self.nseg])
        self.leg.controls=controls
        tof = x[2 * self.nseg]
        # Handle each time encoding separately.
        if self.time_encoding == "uniform":
            # Uniform epoch grid
            self.leg.tgrid = _np.linspace(0, tof, self.nseg + 1)
        elif self.time_encoding == "softmax":
            w = x[1 + 2 * self.nseg : 1 + 2 * self.nseg + self.nseg]
            softmax_weights, _ = _pk.compute_softmax_and_jacobian(w)
            segment_duration = tof * softmax_weights
            # Transform segment durations into a monotonically increasing time grid from 0 to tof.
            tgrid = _np.zeros(self.nseg + 1)
            tgrid[1:] = _np.cumsum(segment_duration)
            self.leg.tgrid = tgrid

    def get_bounds(self):
        if self.time_encoding == "uniform":
            lb = (
                [-_np.pi/2+1e-3, 0.] * self.nseg
                + [self.tof_bounds[0]]
            )
            ub = (
                [_np.pi/2-1e-3, 2*_np.pi] * self.nseg
                + [self.tof_bounds[1]]
            )
        elif self.time_encoding == "softmax":
            lb = (
                [-_np.pi/2+1e-3, 0.] * self.nseg
                + [self.tof_bounds[0]]
                + [self.w_bounds_softmax[0]] * self.nseg
            )
            ub = (
                [_np.pi/2-1e-3, 2*_np.pi] * self.nseg
                + [self.tof_bounds[1]]
                + [self.w_bounds_softmax[1]] * self.nseg
            )
        return (lb, ub)

    def fitness(self, x):
        # We set the leg using data in the decision vector
        self._set_leg_from_x(x)

        # We optimize for minimum time
        obj = x[2 * self.nseg]

        # We compute the equality constraints
        ceq = self.leg.compute_mismatch_constraints()
        retval = _np.array([obj] + ceq)
        return retval

    def gradient(self, x):
        """
        Computes the gradient of the fitness function with respect to the decision vector.

        The fitness function is:
            f(x) = [obj] + mismatch_constraints

        where:
            - obj = tof (minimize time of flight)
            - mismatch_constraints: 6 equality constraints (position and velocity mismatch)

        The decision vector is:
            x = controls + tof (+ weights for softmax)

        where controls = [alpha_1, beta_1, ..., alpha_nseg, beta_nseg]
        and weights = [w_1, w_2, ..., w_nseg] (only if time_encoding is 'softmax')

        Returns:
            gradient: numpy array of shape (nf, nx) where:
                nf = 1 + 6 (total fitness dimension)
                nx = 2*nseg + 1 (+ nseg for softmax) (decision vector dimension)
        """
        if not self.with_gradient:
            raise RuntimeError(
                "Gradient computation requires variational integrator (tas[1] must not be None)"
            )
        # let's start by setting the leg from the decision vector
        self._set_leg_from_x(x)
        # Initialize the gradient matrix.
        nf = 1 + 6 # total number of fitness components
        nx = 2 * self.nseg + 1  # total number of decision variables
        if self.time_encoding == "softmax":
            nx += (
                self.nseg
            )  # account for the additional decision variables for the softmax weights
        gradient = _np.zeros((nf, nx))

        # ============================================
        # Gradient of objective w.r.t decision vector:
        # ============================================
        gradient[0, 2 * self.nseg] = 1.0  # d(tof)/dtof = 1 .. all the rest is zeros

        # ============================================
        # Gradient of mismatch constraints w.r.t decision vector, made of three contributions:
        # ============================================
        ## First contribution -> partials of mismatch constraints w.r.t. controls
        dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid = self.leg.compute_mc_grad()
        for i in range(self.nseg):
            gradient[1:7, 2 * i : 2 * i + 2] = dmc_dcontrols[
                :, 2 * i : 2 * i + 2
            ]
        ## Second contribution -> partials of mismatch constraints w.r.t. time grid
        if self.time_encoding == "uniform":
            # note here that what we are after is:
            # dmc/dtof = dmc/dtgrid * dtgrid/dtof
            # and since for the uniform case tgrid = linspace(0,tof,nseg+1) = [0, tof/nseg, 2*tof/nseg, ..., tof], we have that dtgrid/dtof = [0, 1/nseg, 2/nseg, ..., 1]
            gradient[1:7, -1] = dmcdtgrid @ self._dtgrid_dtof_uniform
        elif self.time_encoding == "softmax":
            # let's start by extracting the softmax weights from the decision vector:
            w = x[1 + 2 * self.nseg : 1 + 2 * self.nseg + self.nseg]
            softmax_weights, J_softmax = _pk.compute_softmax_and_jacobian(w)
            tof = x[2 * self.nseg]

            # gradient w.r.t. tof:
            # tgrid = cumsum(tof*softmax_weights) => dtgrid/dtof = cumsum(softmax_weights)
            # note that also tgrid[0]=0 so dtgrid[0]/dtof = 0
            dtgrid_dtof = _np.zeros(self.nseg + 1)
            dtgrid_dtof[1:] = _np.cumsum(softmax_weights)
            gradient[1:7, 2 * self.nseg] = dmcdtgrid @ dtgrid_dtof

            # now the gradient w.r.t. softmax weights:
            # tgrid[i]=sum_{j=0}^{i-1}(tof*softmax_weights[j]) =>
            # dtgrid[i]/dw[k] = sum_{j=0}^{i-1} tof*dsoftmax_weights[j]/dw[k] =
            #                = tof*sum_{j=0}^{i-1} J_softmax[j,k]
            # now if we assume that: tril_matrix[i,j]=1 if j<i, else 0
            # then dtgrid/dw = tof*(tril_matrix @ J_softmax)
            dtgrid_dw = tof * (self._softmax_cumsum_matrix @ J_softmax)
            gradient[1:7, 1 + 2 * self.nseg : 1 + 2 * self.nseg + self.nseg] = (
                dmcdtgrid @ dtgrid_dw
            )
        else:
            raise NotImplementedError(
                f"Gradient w.r.t. tof not implemented for {self.time_encoding} time encoding"
            )
        return gradient.flatten()

    # If the variational integrator is passed at construction time (i.e., it is not None),
    # `gradient` should be provided (the pagmo UDP interface uses this method to check).
    def has_gradient(self):
        return self.with_gradient

    # Equality constraints only: trajectory mismatches.
    def get_nec(self):
        return 6

    def pretty(self, x):
        """
        Prints a detailed representation of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containing controls, time of flight, and (if time encoding is softmax) the softmax weights.
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
        plot_sail=True,
        sail_size=0.05,
        **kwargs,
    ):
        """
        Plots the trajectory of the zero order hold point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containing controls,
                time of flight, and (if time encoding is softmax) the softmax weights.
                
            *ax* (:class:`matplotlib.axes.Axes`): The matplotlib axes to plot on.
                If None, a new figure and axes will be created.
                
            *N* (:class:`int`): The number of points to plot along the trajectory.
            
            *to_cartesian* (:class:`~collections.abc.Callable`): A function that converts whatever
                state is used in the internal Taylor integrator to Cartesian (r,v).
                
            *mark_segments* (:class:`bool`): Adds markers at each segment edge.
            
            *mark_mismatch* (:class:`bool`): Marks the terminal mismatch point.
            
            *plot_sail* (:class:`bool`): Adds a visualization of a rectangular sail.
            
            *sail_size* (:class:`float`): Half-side length of the rendered sail square,
                in the same units as the trajectory positions.

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the trajectory added.

        Notes:
            If propagation fails on one or more segments, this method still plots the
            successfully propagated forward/backward segments. This is intentional so
            users can visually inspect where and how the trajectory integration breaks.
        """
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        x_arr = _np.asarray(x)
        self._set_leg_from_x(x_arr)
        # Intentionally keep plotting even when propagation is only partial.
        # This helps diagnose failing trajectory segments visually.
        fwd, bck, success = self.leg.get_state_info(N=N)

        nseg = self.nseg
        sail_angles = x_arr[: 2 * nseg].reshape(nseg, 2)  # (nseg, 2): [[alpha, beta], ...]

        if ax is None:
            ax = _pk.plot.make_3Daxis()

        def _set_axes_equal(ax):
            """Force equal aspect ratio on a 3D matplotlib axis."""
            limits = _np.array([
                ax.get_xlim3d(),
                ax.get_ylim3d(),
                ax.get_zlim3d(),
            ])
            origin = _np.mean(limits, axis=1)
            radius = 0.5 * _np.max(_np.abs(limits[:, 1] - limits[:, 0]))
            ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
            ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
            ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

        def _rtn_basis(state_cartesian):
            """Returns (R̂, T̂, N̂) unit vectors from a Cartesian state [rx,ry,rz,vx,vy,vz,...]."""
            r = _np.array(state_cartesian[:3], dtype=float)
            v = _np.array(state_cartesian[3:6], dtype=float)
            R_hat = r / _np.linalg.norm(r)
            h = _np.cross(r, v)
            N_hat = h / _np.linalg.norm(h)
            T_hat = _np.cross(N_hat, R_hat)
            return R_hat, T_hat, N_hat

        def _sail_normal_cartesian(alpha, beta, R_hat, T_hat, N_hat):
            """
            Cone angle alpha: angle between sail normal and R̂ (radial).
            Clock angle beta: azimuth in the T̂-N̂ plane.
            Normal/acceleration direction in RTN (matching zoh.ta._zoh_ss):
            n = [cos(alpha), sin(alpha)*sin(beta), sin(alpha)*cos(beta)]
            """
            n_rtn = _np.array([
                _np.cos(alpha),
                _np.sin(alpha) * _np.sin(beta),
                _np.sin(alpha) * _np.cos(beta),
            ])
            M = _np.column_stack([R_hat, T_hat, N_hat])
            return M @ n_rtn

        def _sail_patch_vertices(center, normal, size):
            """
            Build 4 corners of a square sail of half-side `size` centered at `center`,
            lying in the plane perpendicular to `normal`.
            """
            ref = _np.array([0.0, 0.0, 1.0])
            if abs(_np.dot(normal, ref)) > 0.9:
                ref = _np.array([0.0, 1.0, 0.0])
            u = _np.cross(normal, ref)
            u /= _np.linalg.norm(u)
            v = _np.cross(normal, u)
            v /= _np.linalg.norm(v)
            corners = _np.array([
                center + size * ( u + v),
                center + size * (-u + v),
                center + size * (-u - v),
                center + size * ( u - v),
            ])
            return corners

        def _draw_sail(segment, alpha, beta):
            """Draw sail patch and normal quiver at the midpoint of a segment."""
            mid_idx = len(segment) // 2
            state_mid = to_cartesian(segment[mid_idx])
            center = _np.array(state_mid[:3], dtype=float)
            R_hat, T_hat, N_hat = _rtn_basis(state_mid)
            n_cart = _sail_normal_cartesian(alpha, beta, R_hat, T_hat, N_hat)
            corners = _sail_patch_vertices(center, n_cart, sail_size)
            poly = Poly3DCollection(
                [corners],
                alpha=0.45,
                facecolor='silver',
                edgecolor='dimgray',
                linewidth=0.8,
                zorder=5,
            )
            ax.add_collection3d(poly)
            ax.quiver(
                center[0], center[1], center[2],
                n_cart[0] * sail_size * 2,
                n_cart[1] * sail_size * 2,
                n_cart[2] * sail_size * 2,
                color='tab:blue', linewidth=1.2, arrow_length_ratio=0.3,
            )

        def _sail_active(alpha):
            """Returns True if the sail is delivering a non-zero force."""
            return _np.abs(_np.abs(alpha) - _np.pi / 2) > 1e-3

        # --- Plot forward segments (sail indices 0 … nseg_fwd-1)
        last_point = None
        for i, segment in enumerate(fwd):
            segment_cart = _np.array([to_cartesian(it) for it in segment])
            if mark_segments:
                ax.scatter(
                    segment_cart[0, 0],
                    segment_cart[0, 1],
                    segment_cart[0, 2],
                    **kwargs,
                )
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c='tab:orange')
            last_point = segment_cart[-1, :3]

            if plot_sail:
                alpha, beta = sail_angles[i]
                if _sail_active(alpha):
                    _draw_sail(segment, alpha, beta)

        if mark_mismatch and last_point is not None:
            ax.scatter(
                last_point[0],
                last_point[1],
                last_point[2],
                marker="^",
                **kwargs,
            )

        # --- Plot backward segments (sail indices nseg-1 … nseg_fwd)
        for i, segment in enumerate(bck):
            segment_cart = _np.array([to_cartesian(it) for it in segment])
            if mark_segments:
                ax.scatter(
                    segment_cart[0, 0],
                    segment_cart[0, 1],
                    segment_cart[0, 2],
                    **kwargs,
                )
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c='tab:orange')
            last_point = segment_cart[-1, :3]

            if plot_sail:
                # bck[0] shoots from the terminal state → last sail segment
                seg_idx = nseg - 1 - i
                alpha, beta = sail_angles[seg_idx]
                if _sail_active(alpha):
                    _draw_sail(segment, alpha, beta)

        if mark_mismatch and last_point is not None:
            ax.scatter(
                last_point[0],
                last_point[1],
                last_point[2],
                marker="^",
                **kwargs,
            )

        ax.set_box_aspect([1, 1, 1])
        _set_axes_equal(ax)
        return ax
