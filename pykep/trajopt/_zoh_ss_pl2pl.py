# Copyright (c) 2023-2026 Dario Izzo
# This file is part of the pykep library.
# SPDX-License-Identifier: MPL-2.0
import numpy as _np
import pykep as _pk
from matplotlib import pyplot as _plt


class zoh_ss_pl2pl:
    """Represents the time optimal Solar Sail transfer between two planets using Zero Order Hold (direct) trajectory legs.

    This class mirrors `zoh_pl2pl` but for Solar Sail dynamics. It works
    internally using :class:`~pykep.leg.zoh` with ``dim_dynamics=6`` and
    ``dim_controls=2`` and manipulates the departure epoch, the transfer
    time (tof) and the sail clock and cone angles (controls) so as to link
    two ephemeris-provided planetary states with a Solar Sail trajectory.

    The decision vector follows the planet-to-planet layout::

        x = [t0_fractional, vinf_dep, idep_x, idep_y, idep_z,
             vinf_arr, iarr_x, iarr_y, iarr_z] + controls + [tof] (+ [weights])

    where:
    - ``t0_fractional``: non-dimensional departure time offset in integrator time units (``TIME``); maps to MJD2000 via ``t0 = t0_bounds[0] + x[0] * TIME * SEC2DAY``
    - ``vinf_dep``, ``vinf_arr``: excess velocity magnitudes (non-dimensional, scaled by ``V``)
    - ``idep``, ``iarr``: unit direction vectors of the excess velocities
    - ``controls``: ``[alpha, beta] * nseg`` sail cone/clock angles per segment
    - ``tof``: non-dimensional time of flight (scaled by ``TIME``)
    """

    def __init__(
        self,
        pls=None,
        plf=None,
        t0_bounds=None,
        tof_bounds=None,
        vinf_dep_bounds=None,
        vinf_arr_bounds=None,
        nseg=10,
        cut=0.6,
        tas=None,
        inequalities_for_tc=False,
        time_encoding="uniform",
        w_bounds_softmax=None,
        max_steps=None,
        L=_pk.AU,
        V=_pk.EARTH_VELOCITY,
    ):
        """
        Initializes a planet-to-planet ZOH solar-sail transfer problem with bounded departure and arrival excess velocities.

        Args:
            *pls* (:class:`~pykep.planet`): Departure planet or ephemeris provider.
                It must implement, if gradients are requested,
                ``acc(t_mjd2000)``. The returned position, velocity and acceleration
                units must be consistent with ``L`` and ``V``. Defaults to JPL low-precision Earth.

            *plf* (:class:`~pykep.planet`): Arrival planet or ephemeris provider.
                It follows the same contract as ``pls``. Defaults to JPL low-precision Mars.

            *t0_bounds* (:class:`list`): Lower and upper bounds for the departure epoch,
                expressed in MJD2000 days. Defaults to ``[6700.0, 6800.0]``.

            *tof_bounds* (:class:`list`): Bounds for the time of flight in integrator time units.
               Defaults to ``[3.4, 8.6]``.

            *vinf_dep_bounds* (:class:`list`): Bounds for the magnitude of the departure excess
                velocity in integrator velocity units. Defaults to ``[0.0, 0.2]``.

            *vinf_arr_bounds* (:class:`list`): Bounds for the magnitude of the arrival excess
                velocity in integrator velocity units. Defaults to ``[0.0, 0.2]``.

            *nseg* (:class:`int`): Number of constant-control segments used by the ZOH transcription.
                Defaults to 10.

            *cut* (:class:`float`): Cut parameter passed to :class:`~pykep.leg.zoh`.
                It defines the internal forward/backward split used by the transcription.
                Defaults to 0.6.

            *tas* (:class:`tuple`): `(ta, ta_var)` Taylor-adaptive integrators for solar sail dynamics

                - `ta`: Nominal dynamics (state dim 6, pars ≥ 2)

                - `ta_var`: Variational dynamics (state dim 54, same pars). When None, no gradients will be used.
        

            *inequalities_for_tc* (:class:`bool`): If True, direction-normalization constraints
                are exposed as inequalities; otherwise as equalities.

            *time_encoding* (:class:`str`): Time-grid encoding scheme.
                ``'uniform'`` uses equally spaced segment boundaries over the total time of flight.
                ``'softmax'`` adds ``nseg`` weights to the decision vector and maps them to positive
                segment durations summing to the total time of flight. Defaults to ``'uniform'``.

            *w_bounds_softmax* (:class:`list`): Lower and upper bounds for the softmax weights.
                These bounds are used only when ``time_encoding='softmax'``. Defaults to ``[-1.0, 1.0]``.

            *max_steps* (:class:`int` or ``None``): Maximum number of internal integrator steps
                allowed by :class:`~pykep.leg.zoh`. If ``None``, the integrator default is used.

            *L* (:class:`float`): Length scale used to map ephemeris positions to integrator units.

            *V* (:class:`float`): Velocity scale used to map ephemeris velocities to integrator units.
        """
        if not isinstance(nseg, (int, _np.integer)) or nseg <= 0:
            raise ValueError("nseg must be a positive integer")

        # Initialize defaults for optional constructor arguments.
        if pls is None:
            pls = _pk.planet(_pk.udpla.jpl_lp(body="EARTH"))
        if plf is None:
            plf = _pk.planet(_pk.udpla.jpl_lp(body="MARS"))
        if t0_bounds is None:
            t0_bounds = [6700.0, 6800.0]
        if tof_bounds is None:
            tof_bounds = [3.4, 8.6]
        if vinf_dep_bounds is None:
            vinf_dep_bounds = [0.0, 0.2]
        if vinf_arr_bounds is None:
            vinf_arr_bounds = [0.0, 0.2]
        if w_bounds_softmax is None:
            w_bounds_softmax = [-1.0, 1.0]

        if tas is None:
            tas = (_pk.ta.get_zoh_ss(1e-10), None)

        # Store inputs
        self.pls = pls
        self.plf = plf
        self.nseg = nseg
        self.t0_bounds = t0_bounds
        self.tof_bounds = tof_bounds
        self.vinf_dep_bounds = vinf_dep_bounds
        self.vinf_arr_bounds = vinf_arr_bounds
        self.with_gradient = tas[1] is not None
        self.time_encoding = time_encoding
        self.w_bounds_softmax = w_bounds_softmax
        self.max_steps = max_steps
        self.inequalities_for_tc = inequalities_for_tc

        # Non dimensional units (bridge ephemeris ↔ integrator units)
        self.L = L
        self.V = V
        self.TIME = self.L / self.V
        self.ACC = self.V / self.TIME

        supported_time_encodings = ["uniform", "softmax"]
        if self.time_encoding not in supported_time_encodings:
            raise NotImplementedError(
                f"Only {supported_time_encodings} time encodings are currently implemented"
            )

        # Build and store a ZOH leg instance. Random initial sail angles only initialize
        # the leg object once and are overwritten before fitness/gradient/plot evaluations.
        controls = _np.random.uniform(-_np.pi / 2, _np.pi / 2, (2 * nseg,))

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
        # [t0_frac, vinf_dep, idep_x, idep_y, idep_z, vinf_arr, iarr_x, iarr_y, iarr_z] + controls(2*nseg) + [tof] (+ [weights])
        nx = 9 + 2 * self.nseg + 1
        if self.time_encoding == "softmax":
            nx += self.nseg
        return nx

    def compute_t0(self, x):
        return self.t0_bounds[0] + x[0] * self.TIME * _pk.SEC2DAY

    def _set_leg_from_x(self, x):
        if len(x) != self._expected_nx():
            raise ValueError(
                f"Invalid decision vector length: got {len(x)}, expected {self._expected_nx()}"
            )

        # Parse decision vector (pl2pl layout)
        t0 = self.compute_t0(x)
        vinf_dep_mag = x[1]
        idep = x[2:5]
        vinf_arr_mag = x[5]
        iarr = x[6:9]

        # Controls for sail: 2 per segment
        controls = list(x[9 : 9 + 2 * self.nseg])
        self.leg.controls = controls

        tof = x[9 + 2 * self.nseg]

        # get planet states at departure and arrival (SI units)
        t0_arr = t0 + self.TIME * tof / _pk.DAY2SEC
        rs, vs = self.pls.eph(t0)
        rf, vf = self.plf.eph(t0_arr)

        # convert to non-dimensional
        rs_nd = _np.array([it / self.L for it in rs])
        vs_nd = _np.array([it / self.V for it in vs])
        rf_nd = _np.array([it / self.L for it in rf])
        vf_nd = _np.array([it / self.V for it in vf])

        # apply vinf deltas (non-dimensional)
        dv_dep_nd = [vinf_dep_mag * it for it in idep]
        dv_arr_nd = [vinf_arr_mag * it for it in iarr]

        v_sc_dep_nd = [a + b for a,b in zip(vs_nd, dv_dep_nd)]
        v_sc_arr_nd = [a + b for a,b in zip(vf_nd, dv_arr_nd)]

        # set leg boundary states (solar sail dynamics use 6D states)
        self.leg.state0 = list(rs_nd) + v_sc_dep_nd
        self.leg.state1 = list(rf_nd) + v_sc_arr_nd

        # set time grid
        if self.time_encoding == "uniform":
            self.leg.tgrid = _np.linspace(0, tof, self.nseg + 1)
        elif self.time_encoding == "softmax":
            w = x[10 + 2 * self.nseg : 10 + 2 * self.nseg + self.nseg]
            softmax_weights, _ = _pk.compute_softmax_and_jacobian(w)
            segment_duration = tof * softmax_weights
            tgrid = _np.zeros(self.nseg + 1)
            tgrid[1:] = _np.cumsum(segment_duration)
            self.leg.tgrid = tgrid

        # return ephemerides for possible gradient use
        return rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd

    def get_bounds(self):
        # Bounds follow pl2pl decision vector layout
        if self.time_encoding == "uniform":
            lb = (
                [0.0]
                + [self.vinf_dep_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [self.vinf_arr_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + ([-_np.pi / 2 + 1e-3, 0.0] * self.nseg)
                + [self.tof_bounds[0]]
            )
            ub = (
                [(self.t0_bounds[1]-self.t0_bounds[0]) * _pk.DAY2SEC / self.TIME]
                + [self.vinf_dep_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [self.vinf_arr_bounds[1]]
                + [1.0, 1.0, 1.0]
                + ([_np.pi / 2 - 1e-3, 2 * _np.pi] * self.nseg)
                + [self.tof_bounds[1]]
            )
        elif self.time_encoding == "softmax":
            lb = (
                [0.0]
                + [self.vinf_dep_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + [self.vinf_arr_bounds[0]]
                + [-1.0, -1.0, -1.0]
                + ([-_np.pi / 2 + 1e-3, 0.0] * self.nseg)
                + [self.tof_bounds[0]]
                + [self.w_bounds_softmax[0]] * self.nseg
            )
            ub = (
                [(self.t0_bounds[1]-self.t0_bounds[0]) * _pk.DAY2SEC / self.TIME]
                + [self.vinf_dep_bounds[1]]
                + [1.0, 1.0, 1.0]
                + [self.vinf_arr_bounds[1]]
                + [1.0, 1.0, 1.0]
                + ([_np.pi / 2 - 1e-3, 2 * _np.pi] * self.nseg)
                + [self.tof_bounds[1]]
                + [self.w_bounds_softmax[1]] * self.nseg
            )
        return (lb, ub)

    def fitness(self, x):
        # We set the leg using data in the decision vector
        self._set_leg_from_x(x)

        # Objective: minimize time of flight
        obj = x[9 + 2 * self.nseg]

        # Mismatch constraints (6 equality constraints from the ZOH leg)
        c = self.leg.compute_mismatch_constraints()

        # Direction normalization constraints (always present in the return vector;
        # whether they are equalities or inequalities is decided by get_nec/get_nic)
        idep = x[2:5]
        iarr = x[6:9]
        c += [idep[0] * idep[0] + idep[1] * idep[1] + idep[2] * idep[2] - 1.0]
        c += [iarr[0] * iarr[0] + iarr[1] * iarr[1] + iarr[2] * iarr[2] - 1.0]

        return [obj] + c

    def gradient(self, x):
        if not self.with_gradient:
            raise RuntimeError(
                "Gradient computation requires variational integrator (tas[1] must not be None)"
            )

        # if .acc method of eph not implemented, we cannot compute gradients
        if not hasattr(self.pls, "acc") or not hasattr(self.plf, "acc"):
            raise NotImplementedError(
                "Gradient computation requires .acc method in planet ephemerides for computing accelerations"
            )

        # Set leg and collect ephemerides (SI and non-dimensional)
        rs, vs, rf, vf, rs_nd, vs_nd, rf_nd, vf_nd = self._set_leg_from_x(x)

        vinf_dep_mag = x[1]
        idep = x[2:5]
        vinf_arr_mag = x[5]
        iarr = x[6:9]

        # objective + mismatch(6) + 2 direction normalization constraints
        nf = 1 + 6 + 2
        nx = self._expected_nx()
        gradient = _np.zeros((nf, nx))

        # Objective gradient: d(tof)/dtof = 1
        gradient[0, 9 + 2 * self.nseg] = 1.0

        dmc_dx0, dmc_dx1, dmc_dcontrols, dmcdtgrid = self.leg.compute_mc_grad()

        t0 = self.compute_t0(x)
        tof = x[9 + 2 * self.nseg]

        # Planet accelerations in SI, converted to non-dimensional
        t0_arr = t0 + self.TIME * tof / _pk.DAY2SEC
        as_nd = _np.array(self.pls.acc(t0)) / self.ACC
        af_nd = _np.array(self.plf.acc(t0_arr)) / self.ACC

        # Gradient of mismatch w.r.t. t0:
        #   d(rs_nd)/dt0 = vs_nd * DAY2SEC / TIME
        #   d(vs_nd)/dt0 = as_nd * DAY2SEC / TIME
        dx0_dt0 = _np.concatenate([
            vs_nd * _pk.DAY2SEC / self.TIME,
            as_nd * _pk.DAY2SEC / self.TIME,
        ])
        gradient[1:7, 0] = dmc_dx0 @ dx0_dt0

        # Arrival state also depends on t0 (tf = t0 + tof*TIME/DAY2SEC, so dtf/dt0 = 1)
        dx1_dt0 = _np.concatenate([
            vf_nd * _pk.DAY2SEC / self.TIME,
            af_nd * _pk.DAY2SEC / self.TIME,
        ])
        gradient[1:7, 0] += dmc_dx1 @ dx1_dt0

        # Gradient w.r.t. departure excess velocity magnitude
        # d(v_sc_dep_nd)/d(vinf_dep_mag) = idep  →  d(state0)/d(vinf_dep_mag) only in velocity block
        gradient[1:7, 1] = dmc_dx0[:, 3:6] @ idep

        # Gradient w.r.t. departure direction unit vector
        # d(v_sc_dep_nd)/d(idep) = vinf_dep_mag * I
        gradient[1:7, 2:5] = dmc_dx0[:, 3:6] * vinf_dep_mag

        # Gradient w.r.t. arrival excess velocity magnitude
        gradient[1:7, 5] = dmc_dx1[:, 3:6] @ iarr

        # Gradient w.r.t. arrival direction unit vector
        gradient[1:7, 6:9] = dmc_dx1[:, 3:6] * vinf_arr_mag

        # Controls contribution (controls start at index 9, 2 per segment)
        for i in range(self.nseg):
            gradient[1:7, 9 + 2 * i : 9 + 2 * i + 2] = dmc_dcontrols[:, 2 * i : 2 * i + 2]

        # Time-of-flight contribution (tgrid + moving arrival state)
        # d(rf_nd)/d(tof) = vf_nd,  d(vf_nd)/d(tof) = af_nd  (in non-dimensional time)
        dx1_dtof = _np.concatenate([vf_nd, af_nd])

        if self.time_encoding == "uniform":
            dtgrid_dtof = _np.linspace(0.0, 1.0, self.nseg + 1)
            gradient[1:7, 9 + 2 * self.nseg] = dmcdtgrid @ dtgrid_dtof
            gradient[1:7, 9 + 2 * self.nseg] += dmc_dx1 @ dx1_dtof
        elif self.time_encoding == "softmax":
            w = x[10 + 2 * self.nseg : 10 + 2 * self.nseg + self.nseg]
            softmax_weights, J_softmax = _pk.compute_softmax_and_jacobian(w)
            dtgrid_dtof = _np.zeros(self.nseg + 1)
            dtgrid_dtof[1:] = _np.cumsum(softmax_weights)
            gradient[1:7, 9 + 2 * self.nseg] = dmcdtgrid @ dtgrid_dtof
            gradient[1:7, 9 + 2 * self.nseg] += dmc_dx1 @ dx1_dtof
            tril_matrix = _np.tril(_np.ones((self.nseg, self.nseg)))
            dtgrid_dw = tof * tril_matrix @ J_softmax
            dtgrid_dw = _np.vstack([_np.zeros(self.nseg), dtgrid_dw])
            gradient[1:7, 10 + 2 * self.nseg : 10 + 2 * self.nseg + self.nseg] = (
                dmcdtgrid @ dtgrid_dw
            )

        # Direction normalization constraints gradients
        idep = x[2:5]
        iarr = x[6:9]
        gradient[7, 2:5] = 2.0 * idep
        gradient[8, 6:9] = 2.0 * iarr

        # Scale t0_fractional gradient: d/d(t0_frac) = d/dt0 * dt0/dt0_frac
        gradient[:, 0] *= self.TIME * _pk.SEC2DAY

        return gradient.flatten()

    def has_gradient(self):
        return self.with_gradient

    def get_nec(self):
        # When inequalities_for_tc=False: mismatch(6) + direction_norm(2) are all equalities.
        # When inequalities_for_tc=True:  only mismatch(6) are equalities; direction constraints become inequalities.
        if self.inequalities_for_tc:
            return 6
        else:
            return 6 + 2

    def get_nic(self):
        # When inequalities_for_tc=True: direction normalization constraints (2) are inequalities.
        if self.inequalities_for_tc:
            return 2
        return 0

    def pretty(self, x):
        """
        Prints a detailed representation of the solar-sail planet-to-planet ZOH problem.
        """
        self._set_leg_from_x(x)

        vinf_dep_mag = x[1]
        idep = x[2:5]
        vinf_arr_mag = x[5]
        iarr = x[6:9]

        dv_dep_nd = vinf_dep_mag * _np.array(idep)
        dv_arr_nd = vinf_arr_mag * _np.array(iarr)
        dv_dep_si = dv_dep_nd * self.V
        dv_arr_si = dv_arr_nd * self.V
        vinf_dep_mag_si = vinf_dep_mag * self.V
        vinf_arr_mag_si = vinf_arr_mag * self.V

        tof_nd = x[9 + 2 * self.nseg]
        tof_days = tof_nd * self.TIME / _pk.DAY2SEC
        t0 = self.compute_t0(x)

        print(f"\nSolar Sail ZOH transfer (free velocities, non-dimensional formulation)")
        print(f"Departure: {self.pls.name}\nArrival: {self.plf.name}")
        print(
            f"\nLaunch epoch: {t0:.5f} MJD2000, a.k.a. {_pk.epoch(t0, _pk.epoch.julian_type.MJD2000)}"
        )
        print(
            f"Arrival epoch: {t0 + tof_days:.5f} MJD2000, a.k.a. {_pk.epoch(t0 + tof_days, _pk.epoch.julian_type.MJD2000)}"
        )
        print(f"Time of flight: {tof_nd:.5f} (nd), {tof_days:.5f} days")
        print(f"\nDeparture excess velocity:")
        print(f"  Magnitude: {vinf_dep_mag:.6f} (nd), {vinf_dep_mag_si / 1000:.6f} km/s")
        print(
            f"  Vector (nd): [{dv_dep_nd[0]:.6f}, {dv_dep_nd[1]:.6f}, {dv_dep_nd[2]:.6f}]"
        )
        print(
            f"  Vector (km/s): [{dv_dep_si[0] / 1000:.6f}, {dv_dep_si[1] / 1000:.6f}, {dv_dep_si[2] / 1000:.6f}]"
        )
        print(
            f"  Direction: [{idep[0]:.6f}, {idep[1]:.6f}, {idep[2]:.6f}], norm: {_np.linalg.norm(idep):.6f}"
        )
        print(f"\nArrival excess velocity:")
        print(f"  Magnitude: {vinf_arr_mag:.6f} (nd), {vinf_arr_mag_si / 1000:.6f} km/s")
        print(
            f"  Vector (nd): [{dv_arr_nd[0]:.6f}, {dv_arr_nd[1]:.6f}, {dv_arr_nd[2]:.6f}]"
        )
        print(
            f"  Vector (km/s): [{dv_arr_si[0] / 1000:.6f}, {dv_arr_si[1] / 1000:.6f}, {dv_arr_si[2] / 1000:.6f}]"
        )
        print(
            f"  Direction: [{iarr[0]:.6f}, {iarr[1]:.6f}, {iarr[2]:.6f}], norm: {_np.linalg.norm(iarr):.6f}"
        )
        print(
            f"\nScaling factors: L={self.L / 1e9:.3f} Gm, V={self.V / 1000:.3f} km/s, TIME={self.TIME / 86400:.3f} days"
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
        plot_sail=True,
        sail_size=0.05,
        **kwargs,
    ):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        x_arr = _np.asarray(x)
        self._set_leg_from_x(x_arr)
        fwd, bck, success = self.leg.get_state_info(N=N)

        nseg = self.nseg
        sail_angles = x_arr[9 : 9 + 2 * nseg].reshape(nseg, 2)

        if ax is None:
            ax = _pk.plot.make_3Daxis()

        def _set_axes_equal(ax):
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
            r = _np.array(state_cartesian[:3], dtype=float)
            v = _np.array(state_cartesian[3:6], dtype=float)
            R_hat = r / _np.linalg.norm(r)
            h = _np.cross(r, v)
            N_hat = h / _np.linalg.norm(h)
            T_hat = _np.cross(N_hat, R_hat)
            return R_hat, T_hat, N_hat

        def _sail_normal_cartesian(alpha, beta, R_hat, T_hat, N_hat):
            n_rtn = _np.array([
                _np.cos(alpha),
                _np.sin(alpha) * _np.sin(beta),
                _np.sin(alpha) * _np.cos(beta),
            ])
            M = _np.column_stack([R_hat, T_hat, N_hat])
            return M @ n_rtn

        def _sail_patch_vertices(center, normal, size):
            ref = _np.array([0.0, 0.0, 1.0])
            if abs(_np.dot(normal, ref)) > 0.9:
                ref = _np.array([0.0, 1.0, 0.0])
            u = _np.cross(normal, ref)
            u /= _np.linalg.norm(u)
            v = _np.cross(normal, u)
            v /= _np.linalg.norm(v)
            corners = _np.array([
                center + size * (u + v),
                center + size * (-u + v),
                center + size * (-u - v),
                center + size * (u - v),
            ])
            return corners

        def _draw_sail(segment, alpha, beta):
            mid_idx = len(segment) // 2
            state_mid = segment[mid_idx]
            center = _np.array(state_mid[:3], dtype=float)
            R_hat, T_hat, N_hat = _rtn_basis(state_mid)
            n_cart = _sail_normal_cartesian(alpha, beta, R_hat, T_hat, N_hat)
            corners = _sail_patch_vertices(center, n_cart, sail_size)
            poly = Poly3DCollection(
                [corners],
                alpha=0.45,
                facecolor="silver",
                edgecolor="dimgray",
                linewidth=0.8,
                zorder=5,
            )
            ax.add_collection3d(poly)
            ax.quiver(
                center[0],
                center[1],
                center[2],
                n_cart[0] * sail_size * 2,
                n_cart[1] * sail_size * 2,
                n_cart[2] * sail_size * 2,
                color="tab:blue",
                linewidth=1.2,
                arrow_length_ratio=0.3,
            )

        def _sail_active(alpha):
            return _np.abs(_np.abs(alpha) - _np.pi / 2) > 1e-3

        last_point = None
        for i, segment in enumerate(fwd):
            segment_cart = _np.array(segment)
            if mark_segments:
                ax.scatter(segment_cart[0, 0], segment_cart[0, 1], segment_cart[0, 2], **kwargs)
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c="tab:orange")
            last_point = segment_cart[-1, :3]
            if plot_sail:
                alpha, beta = sail_angles[i]
                if _sail_active(alpha):
                    _draw_sail(segment, alpha, beta)

        if mark_mismatch and last_point is not None:
            ax.scatter(last_point[0], last_point[1], last_point[2], marker="^", **kwargs)

        for i, segment in enumerate(bck):
            segment_cart = _np.array(segment)
            if mark_segments:
                ax.scatter(segment_cart[0, 0], segment_cart[0, 1], segment_cart[0, 2], **kwargs)
            ax.plot(segment_cart[:, 0], segment_cart[:, 1], segment_cart[:, 2], c="tab:orange")
            last_point = segment_cart[-1, :3]
            if plot_sail:
                seg_idx = nseg - 1 - i
                alpha, beta = sail_angles[seg_idx]
                if _sail_active(alpha):
                    _draw_sail(segment, alpha, beta)

        if mark_mismatch and last_point is not None:
            ax.scatter(last_point[0], last_point[1], last_point[2], marker="^", **kwargs)

        ax.set_box_aspect([1, 1, 1])
        _set_axes_equal(ax)
        return ax
