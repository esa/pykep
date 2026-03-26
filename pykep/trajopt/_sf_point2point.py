## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)## 
## This file is part of the kep3 library.## 
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as _np
import pykep as pk


class sf_point2point:
    """Represents the optimal low-thrust transfer between two fixed points using the Sims-Flanagan (direct) method.

    This problem works internally using the :class:`~pykep.leg.sims_flanagan` and manipulates its transfer time T, final mass mf and the controls as to
    link two fixed points in space with a low-thrust trajectory.

    It can be used to better profile and understand performances of optimizers on this type of direct approach, but has a limited use
    in the design of interplanetary trajectories as per the fixed point limitation.

    The decision vector is::

        z = [mf, throttles, tof]

    where throttles is a vector of throttles structures as [u0x, u0y,u0z, ...]. By throttles we intend non dimensiona thrust levels in [0,1].
    """

    def __init__(
        self,
        rvs=[
            _np.array([1, 0.1, -0.1]) * pk.AU,
            _np.array([0.2, 1, -0.2]) * pk.EARTH_VELOCITY,
        ],
        rvf=[
            _np.array([-1.2, -0.1, 0.1]) * pk.AU,
            _np.array([0.2, -1.023, 0.44]) * pk.EARTH_VELOCITY,
        ],
        ms=1000.,
        mu=pk.MU_SUN,
        max_thrust=0.12,
        veff=3000*pk.G0,
        tof_bounds=[80., 400.],
        mf_bounds=[200.0, 1000.0],
        nseg=10,
        cut=0.6,
        mass_scaling=1000,
        r_scaling=pk.AU,
        v_scaling=pk.EARTH_VELOCITY,
        with_gradient=True,
    ):
        """
        Initializes the sf_point2point instance with given parameters.

        Args:
            *rvs* (:class:`list`): Initial position and velocity vectors. Defaults to two vectors scaled by :class:`~pykep.AU` and Earth's velocity.

            *rvf* (:class:`list`): Final position and velocity vectors. Defaults to two vectors scaled by :class:`~pykep.AU` and Earth's velocity.

            *ms* (:class:`float`): Initial spacecraft mass in kg. Defaults to 1000 kg.

            *mu* (:class:`float`): Gravitational parameter, default is for the Sun (:class:`~pykep.MU_SUN`).

            *max_thrust* (:class:`float`): Maximum thrust in Newtons. Defaults to 0.12 N.

            *isp* (:class:`float`): Specific impulse in seconds. Defaults to 3000 s.

            *tof_bounds* (:class:`list`): Bounds for time of flight in days. Defaults to [0, 400] days.

            *mf_bounds* (:class:`list`): Bounds for final mass in kg. Defaults to [200.0, 1000.0] kg.

            *nseg* (:class:`int`): Number of segments for the trajectory. Defaults to 10.

            *cut* (:class:`float`): Cut parameter for the :class:`~pykep.leg.sims_flanagan`. Defaults to 0.6.

            *mass_scaling* (:class:`float`): Scaling factor for mass (used to scale constraints). Defaults to 1000.

            *r_scaling* (:class:`float`): Scaling factor for distance, (used to scale constraints). Defaults AU (:class:`~pykep.AU`).

            *v_scaling* (:class:`float`): Scaling factor for velocity (used to scale constraints). Defaults the Earth's velocity (:class:`~pykep.EARTH_VELOCITY`).

            *with_gradient* (:class:`bool`): Indicates if gradient information should be used. Defaults True.

        """
        # We add as data member one single Sims-Flanagan leg and set it using problem data
        self.leg = pk.leg.sims_flanagan()
        self.leg.rvs = rvs
        self.leg.ms = ms
        self.leg.rvf = rvf
        self.leg.max_thrust = max_thrust
        self.leg.veff = veff
        self.leg.mu = mu
        self.leg.cut = cut
        
        # We define some additional datamembers useful later-on
        self.nseg = nseg
        self.tof_bounds = tof_bounds
        self.mf_bounds = mf_bounds
        self.mass_scaling = mass_scaling
        self.r_scaling = r_scaling
        self.v_scaling = v_scaling
        self.with_gradient = with_gradient
    
    def _set_leg_from_x(self, x):
        # We set the data in the leg using the decision vector
        self.leg.tof = x[-1] * pk.DAY2SEC
        self.leg.mf = x[0]
        self.leg.throttles = x[1:-1]

    def get_bounds(self):
        lb = [self.mf_bounds[0]] + [-1, -1, -1] * self.nseg + [self.tof_bounds[0]]
        ub = [self.mf_bounds[1]] + [1, 1, 1] * self.nseg + [self.tof_bounds[1]]
        return (lb, ub)


    def fitness(self, x):
        # 1 - We set the leg using data in the decision vector
        self._set_leg_from_x(x)
        obj = -x[0] / self.mass_scaling

        # 2 - We compute the constraints violations (mismatch+throttle)
        ceq = self.leg.compute_mismatch_constraints()
        cineq = self.leg.compute_throttle_constraints()
        retval = _np.array([obj] + ceq + cineq)  # here we can sum lists

        # 3 - We scale the values in nd units (numerical solvers are sensitive to well-scaled values)
        retval[1:4] /= self.r_scaling
        retval[4:7] /= self.v_scaling
        retval[7] /= self.mass_scaling

        return retval

    def has_gradient(self):
        return self.with_gradient
    
    def gradient(self, x):
        self._set_leg_from_x(x)
        _, mcg_xf, mcg_th_tof = self.leg.compute_mc_grad()
        tcg_th = self.leg.compute_tc_grad()

        # 1 - The gradient of the objective function (obj = -mf)
        retval = [-1.0 / self.mass_scaling]
        # 2 - The gradient of the mismatch contraints (mcg). We divide them in pos, vel mass as they have different scaling units
        # pos
        for i in range(3):
            # First w.r.t. mf
            retval.append(mcg_xf[i, -1] / self.r_scaling)
            # Then the [throttles, tof]
            retval.extend(mcg_th_tof[i, :] / self.r_scaling)
            retval[-1] *= pk.DAY2SEC
        # vel
        for i in range(3, 6):
            # First w.r.t. mf
            retval.append(mcg_xf[i, -1] / self.v_scaling)
            # Then the [throttles, tof]
            retval.extend(mcg_th_tof[i, :] / self.v_scaling)
            retval[-1] *= pk.DAY2SEC
        # mass
        for i in range(6, 7):
            # First w.r.t. mf
            retval.append(mcg_xf[i, -1] / self.mass_scaling)
            # Then the [throttles, tof]
            retval.extend(mcg_th_tof[i, :] / self.mass_scaling)
            retval[-1] *= pk.DAY2SEC
        # 3 -  The gradient of the throttle constraints
        for i in range(self.nseg):
            retval.extend(tcg_th[i, 3 * i : 3 * i + 3])

        return retval

    def gradient_sparsity(self):
        dim = 2 + 3 * self.nseg
        # The objective function only depends on the final mass, which is in the chromosome.
        retval = [[0, 0]]
        # The mismatch constraints depend on all variables.
        for i in range(1, 8):
            for j in range(dim):
                retval.append([i, j])
        # The throttle constraints only depend on the specific throttles (3).
        for i in range(self.nseg):
            retval.append([8 + i, 3 * i + 1])
            retval.append([8 + i, 3 * i + 2])
            retval.append([8 + i, 3 * i + 3])
        # We return the sparsity pattern
        return retval

    def get_nec(self):
        return 7

    def get_nic(self):
        return self.nseg

    def pretty(self, x):
        """
        Prints a detailed representation of the Point to point problem.

        Args:
            *x* (:class:`list`): The decision vector containing final mass, thrust direction, and time of flight.
        """
        self._set_leg_from_x(x)
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
        **kwargs
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
            **kwargs
        )

        return ax
