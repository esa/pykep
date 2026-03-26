## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)##
## This file is part of the kep3 library.##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as _np
import pykep as pk


class sf_pl2pl:
    """Represents the optimal low-thrust transfer between two :class:`~pykep.planet` using the Sims-Flanagan (direct) method.

    This problem works internally using the :class:`~pykep.leg.sims_flanagan` and manipulates its initial and final states, as well as its transfer time T, final mass mf
    and the controls as to link the two planets with a low-thrust trajectory.

    The particular transcription used is suitable only for few revolutions, after which convergence will start to be problematic.

    The decision vector is::

        z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, throttles, tof] - all in S.I. units except t0 and tof in days

    where throttles is a vector of throttles structured as [u0x, u0y,u0z, ...]. By throttles we intend non dimensiona thrust levels in [0,1].

    """

    def __init__(
        self,
        pls=pk.planet(pk.udpla.jpl_lp(body="EARTH")),
        plf=pk.planet(pk.udpla.jpl_lp(body="MARS")),
        ms=1500,
        mu=pk.MU_SUN,
        max_thrust=0.12,
        veff=3000 * pk.G0,
        t0_bounds=[6700.0, 6800.0],
        tof_bounds=[200.0, 300.0],
        mf_bounds=[1300.0, 1500.0],
        vinfs=3.0,
        vinff=0.0,
        nseg=10,
        cut=0.6,
        mass_scaling=1500,
        r_scaling=pk.AU,
        v_scaling=pk.EARTH_VELOCITY,
        with_gradient=True,
    ):
        """

        Args:
            *pls* (:class:`~pykep.planet`): Initial planet. Defaults to jpl_lp Earth.

            *plf* (:class:`~pykep.planet`): Final planet. Defaults to jpl_lp Mars.

            *ms* (:class:`float`): Initial spacecraft mass in kg. Defaults to 1000 kg.

            *mu* (:class:`float`): Gravitational parameter, default is for the Sun (:class:`~pykep.MU_SUN`).

            *max_thrust* (:class:`float`): Maximum thrust in Newtons. Defaults to 0.12 N.

            *veff* (:class:`float`): Exhaust velocity (Specific impulse * _pk.G0) in m/s. Defaults to 3000 s * _pk.G0.

            *t0_bounds* (:class:`list`): Bounds for departure epoch in MJD2000. Defaults to [6700.0, 6800.0].

            *tof_bounds* (:class:`list`): Bounds for time of flight in days. Defaults to [200, 300] days.

            *mf_bounds* (:class:`list`): Bounds for final mass in kg. Defaults to [1300.0, 1500.0] kg.

            *vinfs* (:class:`float`): Allowed magnitude for the departure's relative velocity in km/s. Defaults to 3.

            *vinff* (:class:`float`): Allowed magnitude for the arrival's relative velocity in km/s. Defaults to 0.

            *nseg* (:class:`int`): Number of segments for the trajectory. Defaults to 10.

            *cut* (:class:`float`): Cut parameter for the :class:`~pykep.leg.sims_flanagan`. Defaults to 0.6.

            *mass_scaling* (:class:`float`): Scaling factor for mass (used to scale constraints). Defaults to 1500.

            *r_scaling* (:class:`float`): Scaling factor for distance, (used to scale constraints). Defaults AU (:class:`~pykep.AU`).

            *v_scaling* (:class:`float`): Scaling factor for velocity (used to scale constraints). Defaults the Earth's velocity (:class:`~pykep.EARTH_VELOCITY`).

            *with_gradient* (:class:`bool`): Indicates if gradient information should be used. Defaults True.

        """
        # We add as data member one single Sims-Flanagan leg and set it using problem data
        self.leg = pk.leg.sims_flanagan()

        self.leg.ms = ms
        self.leg.max_thrust = max_thrust
        self.leg.veff = veff
        self.leg.mu = mu
        self.leg.cut = cut

        # We define some additional datamembers useful later-on
        self.pls = pls
        self.plf = plf
        self.t0_bounds = t0_bounds
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.vinfs = vinfs * 1000  # now in m/s
        self.vinff = vinff * 1000  # now in m/s
        self.nseg = nseg
        self.mass_scaling = mass_scaling
        self.r_scaling = r_scaling
        self.v_scaling = v_scaling
        self.with_gradient = with_gradient

    # z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, throttles, tof]
    def get_bounds(self):
        lb = (
            [self.t0_bounds[0], self.mf_bounds[0]]
            + [-self.vinfs] * 3  # in m/s.
            + [-self.vinff] * 3  # in m/s.
            + [-1, -1, -1] * self.nseg
            + [self.tof_bounds[0]]
        )
        ub = (
            [self.t0_bounds[1], self.mf_bounds[1]]
            + [self.vinfs] * 3  # in m/s.
            + [self.vinff] * 3  # in m/s.
            + [1, 1, 1] * self.nseg
            + [self.tof_bounds[1]]
        )
        return (lb, ub)

    def _set_leg_from_x(self, x):
        # We set the leg using data in the decision vector
        rs, vs = self.pls.eph(x[0])
        rf, vf = self.plf.eph(x[0] + x[-1])
        self.leg.rvs = [rs, [a + b for a, b in zip(vs, x[2:5])]]  # we add vinfs
        self.leg.rvf = [rf, [a + b for a, b in zip(vf, x[5:8])]]  # we add vinff
        self.leg.tof = x[-1] * pk.DAY2SEC
        self.leg.mf = x[1]
        self.leg.throttles = x[8:-1]

        # We return the eph as to avoid having to recompute them later on if needed.
        return rs, vs, rf, vf

    # z = [t0, mf, Vsx, Vsy, Vsz, Vfx, Vfy, Vfz, throttles, tof]
    def fitness(self, x):
        # 1 - We set the optimality principle
        mf = x[1]

        # 2 - We compute the constraints violations (mismatch+throttle)
        self._set_leg_from_x(x)  # set the leg
        ceq = self.leg.compute_mismatch_constraints()
        cineq = self.leg.compute_throttle_constraints()

        # 3 - We add the departure vinfs constraint (quadratic)
        cineq = cineq + [
            (x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - self.vinfs**2) / (self.v_scaling**2)
        ]
        # We add the departure vinff constraint (quadratic)
        cineq = cineq + [
            (x[5] ** 2 + x[6] ** 2 + x[7] ** 2 - self.vinff**2) / (self.v_scaling**2)
        ]

        # 4 - We assemble the return fitness
        retval = _np.array([-mf] + ceq + cineq)  # here we can sum lists

        # 5 - We scale the values in nd units (numerical solvers are sensitive to well-scaled values)
        retval[0] /= self.mass_scaling
        retval[1:4] /= self.r_scaling
        retval[4:7] /= self.v_scaling
        retval[7] /= self.mass_scaling

        return retval

    def has_gradient(self):
        return self.with_gradient

    def gradient(self, x):
        rs, vs, rf, vf = self._set_leg_from_x(x)
        mcg_xs, mcg_xf, mcg_th_tof = self.leg.compute_mc_grad()
        tcg_th = self.leg.compute_tc_grad()

        # 1 - The gradient of the objective function (obj = -mf)
        retval = [-1.0 / self.mass_scaling]

        # 2 - The gradient of the mismatch contraints (mcg). We divide them in pos,vel and mass as
        # they have a different sparsity structure
        # pos-vel - mismatch constraints gradients
        for i in range(6):
            if i < 3:
                scaling = self.r_scaling
            else:
                scaling = self.v_scaling
            # First w.r.t. t0 (it is in days in the decision vector, so we will need to convert)
            tmp = (
                +_np.dot(mcg_xs[i, :3], vs)
                - _np.dot(mcg_xs[i, 3:6], rs) * self.leg.mu / (_np.linalg.norm(rs) ** 3)
                + _np.dot(mcg_xf[i, :3], vf)
                - _np.dot(mcg_xf[i, 3:6], rf) * self.leg.mu / (_np.linalg.norm(rf) ** 3)
            )  # here we assumed that the planet is keplerian (with mu), else there will be a small error in this gradient computation as the acceleration will not exactly be central and mu/r^2.
            retval.append(tmp / scaling / pk.SEC2DAY)
            # Then w.r.t. mf
            retval.append(mcg_xf[i, -1] / scaling)
            # Then w.r.t vinfs
            retval.extend(mcg_xs[i, 3:6] / scaling)
            # Then w.r.t vinff
            retval.extend(mcg_xf[i, 3:6] / scaling)
            # Then the [throttles, tof]
            retval.extend(mcg_th_tof[i, :] / scaling)
            retval[-1] += (
                _np.dot(mcg_xf[i, :3], vf)
                - _np.dot(mcg_xf[i, 3:6], rf) * self.leg.mu / (_np.linalg.norm(rf) ** 3)
            ) / scaling
            retval[-1] *= pk.DAY2SEC
        ## mass - mismatch constraint, here there is no dependency
        ## from t0 nor vinfs, vinff (see sparsity structure) so the code is slightly different
        for i in range(6, 7):
            # First w.r.t. mf
            retval.append(mcg_xf[i, -1] / self.mass_scaling)
            # Then the [throttles, tof]
            retval.extend(mcg_th_tof[i, :] / self.mass_scaling)
            retval[-1] += (
                _np.dot(mcg_xf[i, :3], vf)
                - _np.dot(mcg_xf[i, 3:6], rf) * self.leg.mu / (_np.linalg.norm(rf) ** 3)
            ) / self.mass_scaling
            retval[-1] *= pk.DAY2SEC

        ## 3 -  The gradient of the throttle constraints
        for i in range(self.leg.nseg):
            retval.extend(tcg_th[i, 3 * i : 3 * i + 3])

        ## 4 - The gradient of the vinfs, vinf constraints
        retval.extend([2 * x[2], 2 * x[3], 2 * x[4], 2 * x[5], 2 * x[6], 2 * x[7]])
        retval[-6:] = [a / self.v_scaling**2 for a in retval[-6:]]

        return retval

    def gradient_sparsity(self):
        dim = 9 + 3 * self.nseg
        # The objective function only depends on the final mass, which is in the chromosome.
        retval = [[0, 1]]
        # The mismatch constraints on x,y,z,vx,vy,vz depend on all variables.
        for i in range(1, 7):
            for j in range(dim):
                retval.append([i, j])
        # The mismatch constraints on m depend on all variables except the vinfs, vinff, t0
        retval.append([7, 1])
        for j in range(8, dim):
            retval.append([7, j])
        # The throttle constraints only depend on the specific throttles (3).
        for i in range(self.nseg):
            retval.append([8 + i, 3 * i + 8])
            retval.append([8 + i, 3 * i + 9])
            retval.append([8 + i, 3 * i + 10])
        # The constraints on vinfs, vinff only depend on the vinfs, vinff
        for j in range(2, 5):
            retval.append([8 + self.nseg, j])
        for j in range(5, 8):
            retval.append([9 + self.nseg, j])
        # We return the sparsity pattern
        return retval

    def get_nec(self):
        return 7

    def get_nic(self):
        return self.nseg + 2

    def pretty(self, x):
        """
        Prints a human readable representation of the transfer.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, and time of flight.
        """
        self._set_leg_from_x(x)
        print(f"\nLow-thrust NEP transfer")
        print(f"Departure: {self.pls.get_name()}\nArrival: {self.plf.get_name()}")
        print(
            f"\nLaunch epoch: {x[0]:.5f} MJD2000, a.k.a. {pk.epoch(x[0], pk.epoch.julian_type.MJD2000)}"
        )
        print(
            f"Arrival epoch: {x[0]+x[-1]:.5f} MJD2000, a.k.a. {pk.epoch(x[0]+x[-1], pk.epoch.julian_type.MJD2000)}"
        )

        print(f"Time of flight (days): {x[-1]:.5f} ")
        print(
            f"\nLaunch DV (km/s) {_np.sqrt(x[2] ** 2 + x[3] ** 2 + x[4] ** 2) / 1000:.8f} - [{x[2]/1000.},{x[3]/1000.},{x[4]/1000.}]"
        )
        print(
            f"Arrival DV (km/s) {_np.sqrt(x[5] ** 2 + x[6] ** 2 + x[7] ** 2) / 1000:.8f} - [{x[5]/1000.},{x[6]/1000.},{x[7]/1000.}]"
        )
        print(f"Final mass (kg): {x[1]}")
        print(f"\nDetails on the low-thrust leg: ")
        print(self.leg)

    def plot(
        self,
        x,
        ax=None,
        units=pk.AU,
        show_midpoints=False,
        show_gridpoints=False,
        show_throttles=False,
        length=0.1,
        arrow_length_ratio=0.05,
        **kwargs,
    ):
        """
        Plots the trajectory leg  3D axes.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, and time of flight.

            *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): The 3D axis to plot on. Defaults to None.

            *units* (:class:`float`, optional): The unit scale for the plot. Defaults to _pk.AU.

            *show_midpoints* (:class:`bool`, optional): Whether to show midpoints on the trajectory. Defaults to False.

            *show_gridpoints* (:class:`bool`, optional): Whether to show grid points on the trajectory. Defaults to False.

            *show_throttles* (:class:`bool`, optional): Whether to show throttle vectors. Defaults to False.

            *length* (:class:`float`, optional): Length of the throttle vectors. Defaults to 0.1.

            *arrow_length_ratio* (:class:`float`, optional): Arrow length ratio for the throttle vectors. Defaults to 0.05.

            *\\*\\*kwargs*: Additional keyword arguments for the plot.

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The 3D axis with the plotted trajectory.
        """
        self._set_leg_from_x(x)
        sf = self.leg
        # Making the axis
        if ax is None:
            ax = pk.plot.make_3Daxis(figsize=(7, 7))

        rs, _ = sf.rvs
        rf, _ = sf.rvf
        ax.scatter(rs[0] / pk.AU, rs[1] / units, rs[2] / units, c="k", s=20)
        ax.scatter(rf[0] / pk.AU, rf[1] / units, rf[2] / units, c="k", s=20)

        # Plotting planets
        ax = pk.plot.add_planet(ax, self.pls, when=x[0])
        ax = pk.plot.add_planet_orbit(ax, self.pls, c="gray", alpha=0.5)
        ax = pk.plot.add_planet(ax, self.plf, when=x[0] + x[-1])
        ax = pk.plot.add_planet_orbit(ax, self.plf, c="gray", alpha=0.5)

        # Plotting the trajctory leg
        ax = pk.plot.add_sf_leg(
            ax,
            sf,
            units=units,
            show_throttles=show_throttles,
            length=length,
            show_gridpoints=show_gridpoints,
            show_midpoints=show_midpoints,
            arrow_length_ratio=arrow_length_ratio,
            **kwargs,
        )
        return ax
