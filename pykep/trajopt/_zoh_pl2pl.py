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
from matplotlib import pyplot as _plt


class zoh_pl2pl:
    """Represents the optimal low-thrust transfer between two :class:`~pykep.planet` using the Zero Order Hold (direct) method with free
    departure and arrival velocities and under a generic (low-thrust) dynamics.

    This problem works internally using the :class:`~pykep.leg.zoh` and manipulates its initial and final states, as well as its transfer time T, 
    final mass mf and the controls as to link two orbits described via planets to a low-thrust trajectory. The spacecraft can also be allowed to depart
    and arrive with residual relative velocities.
  
    This class manages two distinct and independent unit systems:

    1. **Ephemeris Units** (from :class:`~pykep.planet`)
        - Time input: MJD2000 days
        - Position/velocity output: user-defined (e.g., SI: meters, meters/second)
        - **Critical contract:** velocities and accelerations must be per *second*

    2. **Integrator Units** (from the Taylor adaptive integrator ``tas``)
        - These are free, typically assuming μ = 1 as well as a length scale

    The transition between systems is controlled by three user-supplied parameters:

    - ``L`` (default: AU): length scale mapping ephemeris position to Integrator Units
    - ``V`` (default: EARTH_VELOCITY): velocity scale mapping ephemeris velocity to Integrator Units
    - ``TIME = L / V``: derived time scale (in SI seconds) converting seconds to Integrator Units
    - ``ACC = V² / L``: derived acceleration scale converting ephemeris acceleration to Integrator Units

    As a consequence, the class applies these standardized conversions to the units used in the backbone numerical integration:

    - Position: $r_{tas} = r_{pla} / L$
    - Velocity: $v_{tas} = v_{pla} / V$
    - Acceleration: $a_{tas} = a_{pla} / \\text{ACC}$
    - Time: $t_{tas} = t_{pla} \\cdot \\text{DAY2SEC} / \\text{TIME}$

    where $t_{pla}$ is in days.

    In a nuthsell: choose L and V to match the units returned by pls/plf, but those ephemeris
    rates must be per second (not per day or something else).

        The decision vector is::

        x = [t0_fractional, mf, vinf_dep, idep_x, idep_y, idep_z, vinf_arr, iarr_x, iarr_y, iarr_z] + controls + [tof] (+ [weights] for softmax)

    where:

    - ``t0_fractional``: non-dimensional, fraction of the ``t0_bounds`` range
    - ``mf``: non-dimensional final mass
    - ``vinf_dep``, ``vinf_arr``: non-dimensional excess velocity magnitudes (scaled by ``V``)
    - ``idep``, ``iarr``: unit Cartesian direction vectors of the relative velocities.
    - ``controls``: ``[T, i_x, i_y, i_z] * nseg`` where ``T`` is non-dimensional throttle.
        The interpretation of ``i_x, i_y, i_z`` is dynamics-dependent: any frame is admissible,
        as long as they are constrained compatibly with throttle constraints.
    - ``tof``: non-dimensional time of flight (scaled by ``TIME``)
    - ``weights``: softmax weights (only if ``time_encoding='softmax'``)

    .. note::

        The API differs from the :class:`~pykep.trajopt.point2point` problem: unit scales (``L``, ``V``) must be explicitly 
        provided to correctly bridge ephemeris and integrator unit systems. Mismatched scales will silently 
        produce incorrect trajectories and gradients.
    """

    def __init__(
        self,
        pls=_pk.planet(_pk.udpla.jpl_lp(body="EARTH")),
        plf=_pk.planet(_pk.udpla.jpl_lp(body="MARS")),
        ms=1.0,
        max_thrust=0.22,
        t0_bounds=[6700.0, 6800.0],
        tof_bounds=[3.4, 8.6],
        mf_bounds=[0.2, 1.0],
        vinf_dep_bounds=[0.0, 0.2],
        vinf_arr_bounds=[0.0, 0.2],
        nseg=10,
        cut=0.6,
        tas=(_pk.ta.get_zoh_kep(1e-10), None),
        cart2state = None,
        state2cart = None,
        inequalities_for_tc = False,
        time_encoding="uniform",
        w_bounds_softmax=[-1.0, 1.0],
        L=_pk.AU,
        V=_pk.EARTH_VELOCITY,
    ):
        """
        Initializes a planet-to-planet ZOH low-thrust transfer problem with free departure and arrival excess velocities.

        Args:
            *pls* (:class:`~pykep.planet`): Departure planet or ephemeris provider.
                It must implement ``eph(t_mjd2000)`` and, if gradients are requested,
                also ``acc(t_mjd2000)``. The returned position, velocity and acceleration
                units must be consistent with ``L`` and ``V``. Defaults to JPL low-precision Earth.

            *plf* (:class:`~pykep.planet`): Arrival planet or ephemeris provider.
                It follows the same contract as ``pls``. Defaults to JPL low-precision Mars.

            *ms* (:class:`float`): Initial spacecraft mass in integrator units.
                This is the mass used in ``self.leg.state0`` and therefore must match the
                mass normalization assumed by the chosen dynamics in ``tas``. Defaults to 1.0.

            *max_thrust* (:class:`float`): Maximum thrust magnitude in integrator units.
                The throttle component stored in each segment of the decision vector is scaled
                by this value before being assigned to the ZOH leg. Defaults to 0.22.

            *t0_bounds* (:class:`list`): Lower and upper bounds for the departure epoch,
                expressed in MJD2000 days. The first decision variable is a fractional value in
                ``[0, 1]`` that is mapped into this interval. Defaults to ``[6700.0, 6800.0]``.

            *tof_bounds* (:class:`list`): Bounds for the time of flight in integrator time units.
                Conversion to days is performed internally using ``TIME / DAY2SEC`` when querying
                the arrival ephemeris. Defaults to ``[3.4, 8.6]``.

            *mf_bounds* (:class:`list`): Bounds for the final spacecraft mass in integrator units.
                The optimization objective is ``-mf``. Defaults to ``[0.2, 1.0]``.

            *vinf_dep_bounds* (:class:`list`): Bounds for the magnitude of the departure excess
                velocity in integrator velocity units, that is, after division by ``V``.
                Defaults to ``[0.0, 0.2]``.

            *vinf_arr_bounds* (:class:`list`): Bounds for the magnitude of the arrival excess
                velocity in integrator velocity units, that is, after division by ``V``.
                Defaults to ``[0.0, 0.2]``.

            *nseg* (:class:`int`): Number of constant-control segments used by the ZOH transcription.
                Each segment contributes four decision variables: throttle magnitude and a direction vector.
                Defaults to 10.

            *cut* (:class:`float`): Cut parameter passed to :class:`~pykep.leg.zoh`.
                It defines the internal forward/backward split used by the transcription.
                Defaults to 0.6.

            *inequalities_for_tc* (:class:`bool`): If True, the throttle constraints are
                formulated as inequalities (``|i_u|**2 - 1 <= 0``), otherwise they are formulated
                as equalities (``|i_u| = 1``). Defaults to False (equality formulation).

            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators

                - `ta`: Nominal dynamics (state dim 7, pars ≥ 4)

                - `ta_var`: Variational dynamics (state dim 84, same pars). When None, no gradients will be used.

            *cart2state* (:class:`list`): Optional list ``[transform, jacobian]`` for non-Cartesian dynamics.
                ``transform`` is a callable that maps a 6D Cartesian state (in integrator units) to the
                6D state expected by the dynamics. ``jacobian`` is a callable that returns the 6×6 derivative
                of ``transform`` w.r.t. the Cartesian state, as a flattened 36-element array (row-major order).
                When provided, the class applies this transformation to boundary states before propagation,
                and uses the Jacobian for gradient computation. Both functions receive the state as a 1D array.
                Defaults to ``None`` (Cartesian dynamics).

            *state2cart* (:class:`callable`): Optional callable for converting dynamics states back to Cartesian.
                It maps a 6D state (in the dynamics' representation) to 6D Cartesian positions and velocities
                (in integrator units). Used by ``plot()`` and ``pretty()`` methods for visualization.
                Ignored if ``cart2state`` is not provided. Defaults to ``None``.

            *time_encoding* (:class:`str`): Time-grid encoding scheme.
                ``'uniform'`` uses equally spaced segment boundaries over the total time of flight.
                ``'softmax'`` adds ``nseg`` weights to the decision vector and maps them to positive
                segment durations summing to the total time of flight. Defaults to ``'uniform'``.

            *w_bounds_softmax* (:class:`list`): Lower and upper bounds for the softmax weights.
                These bounds are used only when ``time_encoding='softmax'``. Defaults to ``[-1.0, 1.0]``.

            *L* (:class:`float`): Length scale used to map ephemeris positions to integrator units.
                If ``pls`` and ``plf`` return positions in meters, a typical choice is ``pykep.AU``
                or another length unit consistent with the integrator dynamics.

            *V* (:class:`float`): Velocity scale used to map ephemeris velocities to integrator units.
                If ephemeris velocities are returned in meters per second, ``V`` must also be expressed
                in meters per second.
        """
        # We define some additional datamembers useful later-on
        self.pls = pls
        self.plf = plf
        self.ms = ms
        self.nseg = nseg
        self.t0_bounds = t0_bounds
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.vinf_dep_bounds = vinf_dep_bounds
        self.vinf_arr_bounds = vinf_arr_bounds
        self.max_thrust = max_thrust
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = w_bounds_softmax
        self.inequalities_for_tc = inequalities_for_tc
        self.cart2state = cart2state
        self.state2cart = state2cart

        # Check user provided functions for cartesian-state conversions if provided
        if self.cart2state is not None:
            if not isinstance(self.cart2state, list):
                raise ValueError("cart2state must be a list containing the conversion and optionally its Jacobian")
            if len(self.cart2state) == 0 or len(self.cart2state) > 2:
                raise ValueError("cart2state must contain at least the conversion function and at most its Jacobian")       
            if not callable(self.cart2state[0]):
                raise ValueError("The first element of cart2state must be a callable function")
            if tas[1] is not None and not callable(self.cart2state[1]):
                raise ValueError("cart2state Jacobian function must be provided when using gradients")
            
        if self.state2cart is not None:
            if not callable(self.state2cart):
                raise ValueError("state2cart must be a callable function")
            
        # Non dimensional units
        # (these are needed to bridge to the calls to the eph method: mjd2000 -> SI, and the dynamics which may be written in non dimensional units)
        self.L = L
        self.V = V
        self.TIME = self.L / self.V
        self.ACC = self.V / self.TIME

        supported_time_encodings = ["uniform", "softmax"]
        if self.time_encoding not in supported_time_encodings:
            raise NotImplementedError(
                f"Only {supported_time_encodings} time encodings are currently implemented"
            )

        # We build and store a ZOH leg as data member.
        controls = _np.random.uniform(-1, 1, (4 * nseg,))
        controls[0::4] = _np.abs(controls[0::4])
        controls[0::4] *= self.max_thrust

        # Get initial states from planets at mid t0 (in SI units)
        t0_mid = (self.t0_bounds[0] + self.t0_bounds[1]) / 2
        tof_mid = (self.tof_bounds[0] + self.tof_bounds[1]) / 2
        rs, vs = self.pls.eph(t0_mid)
        rf, vf = self.plf.eph(t0_mid + self.TIME * tof_mid / _pk.DAY2SEC)

        # Convert to non-dimensional units
        rs_nd = [it / self.L for it in rs]
        vs_nd = [it / self.V for it in vs]
        rf_nd = [it / self.L for it in rf]
        vf_nd = [it / self.V for it in vf]

        states = rs_nd + vs_nd
        statef = rf_nd + vf_nd

        tgrid = _np.linspace(0, tof_mid, self.nseg + 1)

        self.leg = _pk.leg.zoh(
            states + [ms],
            list(controls),
            statef + [ms],
            tgrid,
            cut,
            tas,
        )

    def _compute_throttle_constraints(self):
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

    def _compute_throttle_jacobian_leg(self):
        dtc_dcontrols = _np.zeros((self.nseg, 4 * self.nseg))
        controls = self.leg.controls
        for i in range(self.nseg):
            base = 4 * i
            dtc_dcontrols[i, base + 1] = 2.0 * controls[base + 1]
            dtc_dcontrols[i, base + 2] = 2.0 * controls[base + 2]
            dtc_dcontrols[i, base + 3] = 2.0 * controls[base + 3]
        return dtc_dcontrols

    def compute_t0(self, x):
        return self.t0_bounds[0] + x[0] * (self.t0_bounds[1] - self.t0_bounds[0])

    def _set_leg_from_x(self, x):
        """
        Here is where the decision vector gets decoded into the leg initial/final states, tgrid, mf and controls
        x = [t0, mf, vinf_dep, idep_x, idep_y, idep_z, vinf_arr, iarr_x, iarr_y, iarr_z] + controls + [tof] (+ [weights] for softmax)

        All quantities are in non-dimensional units.
        """
        t0 = self.compute_t0(x)
        state1 = list(self.leg.state1)
        state1[-1] = x[1]  # mf (non-dimensional)

        # Extract velocity magnitudes and directions (non-dimensional)
        vinf_dep_mag = x[2]
        idep = x[3:6]
        vinf_arr_mag = x[6]
        iarr = x[7:10]

        # Compute velocity differences in non-dimensional units
        dv_dep_nd = vinf_dep_mag * _np.array(idep)
        dv_arr_nd = vinf_arr_mag * _np.array(iarr)

        # Set controls (scale throttle entries before assigning to the C++ list-backed property)
        controls = _np.array(x[10 : 10 + 4 * self.nseg], copy=True)
        controls[0::4] *= self.max_thrust
        self.leg.controls = controls.tolist()

        tof = x[10 + 4 * self.nseg]  # non-dimensional

        # Get planet states at departure and arrival (SI units)
        t0_arr = t0 + self.TIME * tof / _pk.DAY2SEC  # Convert nd tof to days
        rs, vs = self.pls.eph(t0)
        rf, vf = self.plf.eph(t0_arr)

        # Convert planet states to non-dimensional
        rs_nd = _np.array([it / self.L for it in rs])
        vs_nd = _np.array([it / self.V for it in vs])
        rf_nd = _np.array([it / self.L for it in rf])
        vf_nd = _np.array([it / self.V for it in vf])

        # Update initial and final states with velocity differences (all non-dimensional)
        v_sc_dep_nd = vs_nd + dv_dep_nd
        v_sc_arr_nd = vf_nd + dv_arr_nd

        # Apply coordinate transformation if provided
        if self.cart2state is not None:
            cart_state0 = list(rs_nd) + list(v_sc_dep_nd)
            cart_state1 = list(rf_nd) + list(v_sc_arr_nd)
            state0_transformed = self.cart2state[0](cart_state0)
            state1_transformed = self.cart2state[0](cart_state1)
            self.leg.state0 = list(state0_transformed) + [self.ms]
            state1[:6] = list(state1_transformed)
        else:
            self.leg.state0 = list(rs_nd) + list(v_sc_dep_nd) + [self.ms]
            state1[:6] = list(rf_nd) + list(v_sc_arr_nd)
        
        self.leg.state1 = state1

        # Set time grid based on encoding (non-dimensional time)
        if self.time_encoding == "uniform":
            self.leg.tgrid = _np.linspace(0, tof, self.nseg + 1)
        elif self.time_encoding == "softmax":
            w = x[11 + 4 * self.nseg : 11 + 4 * self.nseg + self.nseg]
            softmax_weights, _ = _pk.compute_softmax_and_jacobian(w)
            segment_duration = tof * softmax_weights
            tgrid = _np.zeros(self.nseg + 1)
            tgrid[1:] = _np.cumsum(segment_duration)
            self.leg.tgrid = tgrid

        # Return ephemerides in SI and non-dimensional for gradient computation
        return rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd

    def get_bounds(self):
        if self.time_encoding == "uniform":
            lb = (
                [0.0, self.mf_bounds[0]]
                + [self.vinf_dep_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [self.vinf_arr_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [0, -1.0, -1.0, -1.0] * self.nseg
                + [self.tof_bounds[0]]
            )
            ub = (
                [1.0, self.mf_bounds[1]]
                + [self.vinf_dep_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [self.vinf_arr_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [1.0, 1.0, 1.0, 1.0] * self.nseg
                + [self.tof_bounds[1]]
            )
        elif self.time_encoding == "softmax":
            lb = (
                [0.0, self.mf_bounds[0]]
                + [self.vinf_dep_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [self.vinf_arr_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [0, -1.0, -1.0, -1.0] * self.nseg
                + [self.tof_bounds[0]]
                + [self.w_bounds_softmax[0]] * self.nseg
            )
            ub = (
                [1.0, self.mf_bounds[1]]
                + [self.vinf_dep_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [self.vinf_arr_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [1.0, 1.0, 1.0, 1.0] * self.nseg
                + [self.tof_bounds[1]]
                + [self.w_bounds_softmax[1]] * self.nseg
            )
        return (lb, ub)

    def fitness(self, x):
        # We set the leg using data in the decision vector
        self._set_leg_from_x(x)

        # We optimize for maximum final mass (minimum propellant)
        obj = -x[1]

        # We compute the equality constraints
        ceq = self.leg.compute_mismatch_constraints()
        ceq += self._compute_throttle_constraints()

        # Add direction normalization constraints
        idep = x[3:6]
        iarr = x[7:10]
        ceq += [_np.dot(idep, idep) - 1.0]
        ceq += [_np.dot(iarr, iarr) - 1.0]

        retval = _np.array([obj] + ceq)
        return retval

    def gradient(self, x):
        """
        Computes the gradient of the fitness function with respect to the decision vector.

        All computations are in non-dimensional units except t0 which is in MJD2000.
        """
        if not self.with_gradient:
            raise RuntimeError(
                "Gradient computation requires variational integrator (tas[1] must not be None)"
            )

        # if .acc method of eph not implemented, we cannot compute gradients
        if not hasattr(self.pls, "acc") or not hasattr(self.plf, "acc"):
            raise NotImplementedError(
                "Gradient computation requires .acc method in planet ephemerides for computing accelerations"
            )

        # Set the leg from the decision vector
        rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd = self._set_leg_from_x(x)

        # Extract components for clarity
        vinf_dep_mag = x[2]
        idep = x[3:6]
        vinf_arr_mag = x[6]
        iarr = x[7:10]

        # Initialize the gradient matrix
        nf = 1 + 7 + self.nseg + 2
        nx = 10 + 4 * self.nseg + 1
        if self.time_encoding == "softmax":
            nx += self.nseg

        gradient = _np.zeros((nf, nx))

        # Gradient of objective w.r.t decision vector
        gradient[0, 1] = -1.0

        # Gradient of mismatch constraints w.r.t decision vector
        dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid = self.leg.compute_mc_grad()

        t0 = self.compute_t0(x)
        tof = x[10 + 4 * self.nseg]

        # Compute planet accelerations in SI units
        as_si = _np.array(self.pls.acc(t0))
        af_si = _np.array(self.plf.acc(t0 + self.TIME * tof / _pk.DAY2SEC))

        # Convert accelerations to non-dimensional
        as_nd = as_si / self.ACC
        af_nd = af_si / self.ACC

        # Partials of mismatch constraints w.r.t. t0
        # d(rs_nd)/dt0 = d(rs/L)/dt0_days = vs * DAY2SEC / L = vs_nd * (V/L) * DAY2SEC
        # d(vs_nd)/dt0 = (1/V) * as * DAY2SEC = (as_nd * V/TIME) / V * DAY2SEC = as_nd * DAY2SEC / TIME
        dx0_dt0_cart = _np.concatenate(
            [
                vs_nd * _pk.DAY2SEC / self.TIME,
                as_nd * _pk.DAY2SEC / self.TIME,
                [0.0],
            ]
        )
        
        # Apply Jacobian of cart2state transformation if provided
        if self.cart2state is not None:
            cart_state0 = _np.concatenate([rs_nd, vs_nd])
            J0 = self.cart2state[1](cart_state0).reshape((6, 6))
            dx0_dt0_nd = J0 @ dx0_dt0_cart[:6]
            dx0_dt0_nd = _np.concatenate([dx0_dt0_nd, [0.0]])
        else:
            dx0_dt0_nd = dx0_dt0_cart
        
        gradient[1:8, 0] = dmc_dx0 @ dx0_dt0_nd

        # Also account for final state change through t0 (since tf = t0 + tof * TIME / DAY2SEC)
        # dtf/dt0 = 1, so d(rf_nd)/dt0 = vf_nd * (V/L) * DAY2SEC, d(vf_nd)/dt0 = af_nd * DAY2SEC/TIME
        dx1_dt0_cart = _np.concatenate(
            [
                vf_nd * _pk.DAY2SEC / self.TIME,
                af_nd * _pk.DAY2SEC / self.TIME,
                [0.0],
            ]
        )
        
        # Apply Jacobian of cart2state transformation if provided
        if self.cart2state is not None:
            cart_state1 = _np.concatenate([rf_nd, vf_nd])
            J1 = self.cart2state[1](cart_state1).reshape((6, 6))
            dx1_dt0_nd = J1 @ dx1_dt0_cart[:6]
            dx1_dt0_nd = _np.concatenate([dx1_dt0_nd, [0.0]])
        else:
            dx1_dt0_nd = dx1_dt0_cart
        
        gradient[1:8, 0] += dmc_dx1 @ dx1_dt0_nd

        # Partials w.r.t. final mass
        gradient[1:8, 1] = dmc_dx1[:, 6]

        # Partials w.r.t. departure velocity magnitude
        if self.cart2state is not None:
            cart_state0 = _np.concatenate([rs_nd, vs_nd])
            J0 = self.cart2state[1](cart_state0).reshape((6, 6))
            gradient[1:8, 2] = (dmc_dx0 @ _np.vstack([_np.zeros((3, 6)), J0 @ _np.eye(6)[:, 3:6]])) @ idep
        else:
            gradient[1:8, 2] = dmc_dx0[:, 3:6] @ idep

        # Partials w.r.t. departure direction
        if self.cart2state is not None:
            cart_state0 = _np.concatenate([rs_nd, vs_nd])
            J0 = self.cart2state[1](cart_state0).reshape((6, 6))
            gradient[1:8, 3:6] = (dmc_dx0 @ _np.vstack([_np.zeros((3, 6)), J0 @ _np.eye(6)[:, 3:6]])) * vinf_dep_mag
        else:
            gradient[1:8, 3:6] = dmc_dx0[:, 3:6] * vinf_dep_mag

        # Partials w.r.t. arrival velocity magnitude
        if self.cart2state is not None:
            cart_state1 = _np.concatenate([rf_nd, vf_nd])
            J1 = self.cart2state[1](cart_state1).reshape((6, 6))
            gradient[1:8, 6] = (dmc_dx1 @ _np.vstack([_np.zeros((3, 6)), J1 @ _np.eye(6)[:, 3:6]])) @ iarr
        else:
            gradient[1:8, 6] = dmc_dx1[:, 3:6] @ iarr

        # Partials w.r.t. arrival direction
        gradient[1:8, 7:10] = dmc_dx1[:, 3:6] * vinf_arr_mag

        # Partials w.r.t. controls
        for i in range(self.nseg):
            gradient[1:8, 10 + 4 * i : 10 + 4 * i + 4] = dmc_dcontrols[
                :, 4 * i : 4 * i + 4
            ]
            gradient[1:8, 10 + 4 * i] *= self.max_thrust

        # Partials w.r.t. time of flight
        if self.time_encoding == "uniform":
            # dtgrid/dtof in non-dimensional time
            dtgrid_dtof = _np.linspace(0, 1, self.nseg + 1)
            gradient[1:8, -1] = dmcdtgrid @ dtgrid_dtof

            # Add contribution from final state change w.r.t tof (non-dimensional)
            # d(rf_nd)/d(tof_nd) = vf_nd * (V/L) * TIME  [simplifies to vf_nd for V = L/TIME]
            # d(vf_nd)/d(tof_nd) = af_nd
            dx1_dtof_nd = _np.concatenate(
                [vf_nd, af_nd, [0.0]]
            )
            gradient[1:8, -1] += dmc_dx1 @ dx1_dtof_nd

        elif self.time_encoding == "softmax":
            w = x[11 + 4 * self.nseg : 11 + 4 * self.nseg + self.nseg]
            softmax_weights, J_softmax = _pk.compute_softmax_and_jacobian(w)

            # Gradient w.r.t. tof (non-dimensional)
            dtgrid_dtof = _np.zeros(self.nseg + 1)
            dtgrid_dtof[1:] = _np.cumsum(softmax_weights)
            gradient[1:8, 10 + 4 * self.nseg] = dmcdtgrid @ dtgrid_dtof

            # Add contribution from final state change
            dx1_dtof_nd = _np.concatenate(
                [vf_nd, af_nd, [0.0]]
            )
            gradient[1:8, 10 + 4 * self.nseg] += dmc_dx1 @ dx1_dtof_nd

            # Gradient w.r.t. softmax weights
            tril_matrix = _np.tril(_np.ones((self.nseg, self.nseg)))
            dtgrid_dw = tof * tril_matrix @ J_softmax
            dtgrid_dw = _np.vstack([_np.zeros(self.nseg), dtgrid_dw])
            gradient[1:8, 11 + 4 * self.nseg : 11 + 4 * self.nseg + self.nseg] = (
                dmcdtgrid @ dtgrid_dw
            )

        # Gradient of throttle constraints
        gradient[8 : 8 + self.nseg, 0:10] = 0.0
        dtc_dcontrols = self._compute_throttle_jacobian_leg()
        gradient[8 : 8 + self.nseg, 10 : 10 + 4 * self.nseg] = dtc_dcontrols

        # Gradient of direction constraints
        gradient[8 + self.nseg, 3:6] = 2.0 * idep
        gradient[9 + self.nseg, 7:10] = 2.0 * iarr

        # Scaling the gradient to t_fractional d/dt_frac = d/dt0 * dt0/dt_frac
        gradient[:, 0] *= self.t0_bounds[1] - self.t0_bounds[0]

        return gradient.flatten()

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        if self.inequalities_for_tc:
            # Throttle constraints and boundary impulses are inequalities: |i_u|^2 - 1 <= 0
            return 7 
        else:
            return 7 + self.nseg + 2
        
    def get_nic(self):
        if self.inequalities_for_tc:
            # Throttle constraints and boundary impulses are inequalities: |i_u|^2 - 1 <= 0
            return self.nseg + 2
        else:
            return 0

    def pretty(self, x):
        """
        Prints a detailed representation of the zero order hold planet to planet problem with free velocities.
        """
        self._set_leg_from_x(x)

        # Extract velocity components (non-dimensional)
        vinf_dep_mag = x[2]
        idep = x[3:6]
        vinf_arr_mag = x[6]
        iarr = x[7:10]

        dv_dep_nd = vinf_dep_mag * idep
        dv_arr_nd = vinf_arr_mag * iarr

        # Convert to SI for display
        dv_dep_si = dv_dep_nd * self.V
        dv_arr_si = dv_arr_nd * self.V
        vinf_dep_mag_si = vinf_dep_mag * self.V
        vinf_arr_mag_si = vinf_arr_mag * self.V

        tof_nd = x[10 + 4 * self.nseg]
        tof_days = tof_nd * self.TIME / _pk.DAY2SEC

        t0 = self.compute_t0(x)

        print(
            f"\nLow-thrust ZOH transfer (free velocities, non-dimensional formulation)"
        )
        print(f"Departure: {self.pls.name}\nArrival: {self.plf.name}")
        print(
            f"\nLaunch epoch: {t0:.5f} MJD2000, a.k.a. {_pk.epoch(t0, _pk.epoch.julian_type.MJD2000)}"
        )
        print(
            f"Arrival epoch: {t0+tof_days:.5f} MJD2000, a.k.a. {_pk.epoch(t0+tof_days, _pk.epoch.julian_type.MJD2000)}"
        )
        print(f"Time of flight: {tof_nd:.5f} (nd), {tof_days:.5f} days")
        print(f"\nDeparture excess velocity:")
        print(f"  Magnitude: {vinf_dep_mag:.6f} (nd), {vinf_dep_mag_si/1000:.6f} km/s")
        print(
            f"  Vector (nd): [{dv_dep_nd[0]:.6f}, {dv_dep_nd[1]:.6f}, {dv_dep_nd[2]:.6f}]"
        )
        print(
            f"  Vector (km/s): [{dv_dep_si[0]/1000:.6f}, {dv_dep_si[1]/1000:.6f}, {dv_dep_si[2]/1000:.6f}]"
        )
        print(
            f"  Direction: [{idep[0]:.6f}, {idep[1]:.6f}, {idep[2]:.6f}], norm: {_np.linalg.norm(idep):.6f}"
        )
        print(f"\nArrival excess velocity:")
        print(f"  Magnitude: {vinf_arr_mag:.6f} (nd), {vinf_arr_mag_si/1000:.6f} km/s")
        print(
            f"  Vector (nd): [{dv_arr_nd[0]:.6f}, {dv_arr_nd[1]:.6f}, {dv_arr_nd[2]:.6f}]"
        )
        print(
            f"  Vector (km/s): [{dv_arr_si[0]/1000:.6f}, {dv_arr_si[1]/1000:.6f}, {dv_arr_si[2]/1000:.6f}]"
        )
        print(
            f"  Direction: [{iarr[0]:.6f}, {iarr[1]:.6f}, {iarr[2]:.6f}], norm: {_np.linalg.norm(iarr):.6f}"
        )
        print(f"\nFinal mass: {x[1]:.6f} (nd)")
        print(
            f"\nScaling factors: L={self.L/1e9:.3f} Gm, V={self.V/1000:.3f} km/s, TIME={self.TIME/86400:.3f} days"
        )
        print(f"\nDetails on the ZOH leg:")
        print(self.leg)

    def plot(
        self,
        x,
        ax=None,
        N=30,
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

            *mark_segments* (:class:`bool`): adds markers ath each segment edge

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The modified Axes object with the Lambert's problem trajectory added.
        """
        x_arr = _np.asarray(x)
        # to be replaced with a plot method akin the sims-flanagan point2point one
        self._set_leg_from_x(x_arr)
        fwd, bck, _ = self.leg.get_state_info(N=N)
        
        # Convert trajectories back to Cartesian if needed
        if self.state2cart is not None:
            fwd_cart = []
            for segment in fwd:
                segment_cart = _np.array([_np.concatenate([self.state2cart(state[:6]), state[6:]]) for state in segment])
                fwd_cart.append(segment_cart)
            fwd = fwd_cart
            
            bck_cart = []
            for segment in bck:
                segment_cart = _np.array([_np.concatenate([self.state2cart(state[:6]), state[6:]]) for state in segment])
                bck_cart.append(segment_cart)
            bck = bck_cart
        
        # compute the color scheme
        throttles = x_arr[10 : 10 + 4 * self.nseg : 4]
        if ax is None:
            ax = _pk.plot.make_3Daxis()
        # plot
        for i, segment in enumerate(fwd):
            color = (
                0.25 + (0.80 - 0.25) * throttles[i],
                0.41 + (0.36 - 0.41) * throttles[i],
                0.88 + (0.36 - 0.88) * throttles[i],
            )
            if mark_segments:
                ax.scatter(
                    segment[0, 0] * self.L,
                    segment[0, 1] * self.L,
                    segment[0, 2] * self.L,
                    **kwargs,
                )
            ax.plot(
                segment[:, 0] * self.L,
                segment[:, 1] * self.L,
                segment[:, 2] * self.L,
                c=color,
            )
        if mark_mismatch:
            ax.scatter(
                segment[-1, 0] * self.L,
                segment[-1, 1] * self.L,
                segment[-1, 2] * self.L,
                marker="^",
                **kwargs,
            )
        for i, segment in enumerate(bck):
            color = (
                0.25 + (0.80 - 0.25) * throttles[-1 - i],
                0.41 + (0.36 - 0.41) * throttles[-1 - i],
                0.88 + (0.36 - 0.88) * throttles[-1 - i],
            )
            if mark_segments:
                ax.scatter(
                    segment[0, 0] * self.L,
                    segment[0, 1] * self.L,
                    segment[0, 2] * self.L,
                    **kwargs,
                )
            ax.plot(segment[:, 0] * self.L, segment[:, 1] * self.L, segment[:, 2] * self.L, c=color)
        if mark_mismatch:
            ax.scatter(
                segment[-1, 0] * self.L, segment[-1, 1] * self.L, segment[-1, 2] * self.L, marker="^", **kwargs
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
        tgrid = _np.asarray(self.leg.tgrid)  # length nseg+1
        throttles = _np.asarray(x[10 : 10 + 4 * self.nseg : 4])  # length nseg

        # Repeat the last throttle so that the last tgrid point is included in the step plot
        u_plot = _np.r_[throttles, throttles[-1]]

        ax.step(tgrid, u_plot, where="post", **kwargs)

        ax.set_xlabel("time grid")
        ax.set_ylabel("throttle value (nd)")
        return ax
