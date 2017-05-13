from PyGMO.problem import base as base_problem
from PyKEP.core import epoch, propagate_taylor, EARTH_VELOCITY, AU, MU_SUN, MU_EARTH, G0, DAY2SEC
from PyKEP.planet import jpl_lp
from PyKEP.sims_flanagan._sims_flanagan import spacecraft
import numpy as np

class lt_margo(base_problem):
    """
    This class is a PyGMO problem representing a low-thrust interplanetary trajectory from the Earth (or from the
    Sun-Earth L1 or L2) to a target NEO.
    The trajectory is modeled using the Sims-Flanagan model, extended to include the Earth's gravity (which is
    assumed constant along each segment). It is possible to use a non-uniform partition with a denser grid of
    segments in the initial part of the trajectory, in order to have a better sampling there.

    This problem was developed during the M-ARGO CDF study.

    The decision vector (chromosome) is::

      [t0, tof, mf] + [throttles1] + [throttles2] + ...
    """
    def __init__(self,
                 target = jpl_lp('mars'),
                 n_seg = 20,
                 grid_type = "uniform",
                 t0 = [epoch(0), epoch(1000)],
                 tof = [200, 500],
                 m0 = 50.0,
                 Tmax = 0.001,
                 Isp = 2000.0,
                 earth_gravity = False,
                 start = "earth"):
        """
        prob = lt_margo(target = jpl_lp('mars'), n_seg = 20, grid_type = "uniform", t0 = [epoch(0), epoch(1000)],
        tof = [200, 500], m0 = 50.0, Tmax = 0.001, Isp = 2000.0, earth_gravity = False, start = "earth")

        - target: target PyKEP.planet
        - n_seg: number of segments
        - grid_type: "uniform" for uniform segments, "nonuniform" for a denser grid in the first part of the trajectory
        - t0: list of two epochs defining the launch window
        - tof: list of two floats definind the minimum and maximum time of flight (days)
        - m0: initial mass of the spacecraft
        - Tmax: maximum thrust
        - Isp: engine specific impulse
        - earth_gravity: boolean specifying whether to take Earth's gravity into account
        - start: starting point ("earth", "l1", or "l2").

        .. note::

        L1 and L2 are approximated as the points on the line connecting the Sun and the Earth at a distance of, respectively, 0.99 and 1.01 AU from the Sun.

        .. note::

        If the Earth's gravity is enabled, the starting point cannot be the Earth
        """

        # Various checks
        if start not in ["earth", "l1", "l2"]:
            raise ValueError("start must be either 'earth', 'l1' or 'l2'")
        if grid_type not in ["uniform", "nonuniform"]:
            raise ValueError("grid_type must be either 'uniform' or 'nonuniform'")
        if earth_gravity and start == "earth":
            raise ValueError("If Earth gravity is enabled the starting point cannot be the Earth")

        # 1a) The decision vector length ([t0, tof, mf, [Tx,Ty,Tz] * n_seg])
        dim = 3 + n_seg * 3
        # 1b) The total number of constraints (mismatch + throttles)
        c_dim = 7 + n_seg
        # 1c) The number of inequality constraints (throttles)
        c_ineq_dim = n_seg
        # First we call the constructor for the base PyGMO problem
        super(lt_margo, self).__init__(dim, 0, 1, c_dim, c_ineq_dim, 1e-4)

        # 2) We then define some class data members
        # public:
        self.target = target
        # private:
        self.__n_seg = n_seg
        self.__grid_type = grid_type # gridding function (must be 0 in 0, 1 in 1, and strictly increasing)
        self.__sc = spacecraft(m0, Tmax, Isp)
        self.__earth = jpl_lp('earth')
        self.__earth_gravity = earth_gravity
        self.__start = start
        # grid construction
        if grid_type == "uniform":
            grid = np.array([i/n_seg for i in range(n_seg + 1)])
        elif grid_type == "nonuniform":
            grid_f = lambda x : x**2 if x<0.5 else 0.25+1.5*(x-0.5) # quadratic in [0,0.5], linear in [0.5,1]
            grid = np.array([grid_f(i/n_seg) for i in range(n_seg + 1)])
        fwd_seg = np.searchsorted(grid, 0.5, side='right') # index corresponding to the middle of the transfer
        bwd_seg = n_seg - fwd_seg
        fwd_grid = grid[:fwd_seg + 1]
        bwd_grid = grid[fwd_seg:]
        self.__fwd_seg = fwd_seg
        self.__fwd_grid = fwd_grid
        self.__fwd_dt = np.array([(fwd_grid[i+1] - fwd_grid[i]) for i in range(fwd_seg)]) * DAY2SEC
        self.__bwd_seg = bwd_seg
        self.__bwd_grid = bwd_grid
        self.__bwd_dt = np.array([(bwd_grid[i+1] - bwd_grid[i]) for i in range(bwd_seg)]) * DAY2SEC

        # 3) We compute the bounds
        lb = [t0[0].mjd2000] + [tof[0]] + [0] + [-1, -1, -1] * n_seg
        ub = [t0[1].mjd2000] + [tof[1]] + [m0] + [1, 1, 1] * n_seg

        # 4) And we set the bounds
        self.set_bounds(lb, ub)

    # Objective function
    def _objfun_impl(self, x):
        return (-x[2],)

    # Constraints function
    def _compute_constraints_impl(self, x):
        ceq = list()
        cineq = list()
        throttles_constraints = []
        mismatch_constraints = []

        throttles = [x[3 + 3 * i : 6 + 3 * i] for i in range(self.__n_seg)]
        for t in throttles:
            throttles_constraints.append(t[0]**2 + t[1]**2 + t[2]**2 - 1.)

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd = self._propagate(x)[0:6]

        mismatch_constraints.extend([a - b for a,b in zip(rfwd[-1], rbwd[0])])
        mismatch_constraints.extend([a - b for a,b in zip(vfwd[-1], vbwd[0])])
        mismatch_constraints.append(mfwd[-1] - mbwd[0])

        ceq.extend(mismatch_constraints)
        cineq.extend(throttles_constraints)

        # Making the mismatches non dimensional
        ceq[0] /= AU
        ceq[1] /= AU
        ceq[2] /= AU
        ceq[3] /= EARTH_VELOCITY
        ceq[4] /= EARTH_VELOCITY
        ceq[5] /= EARTH_VELOCITY
        ceq[6] /= self.__sc.mass

        # We assemble the constraint vector
        retval = list()
        retval.extend(ceq)
        retval.extend(cineq)

        return retval

    # Propagates the trajectory
    def _propagate(self, x):
        # 1 - We decode the chromosome
        t0 = x[0]
        T = x[1]
        m_f = x[2]
        # We extract the number of segments for forward and backward propagation
        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        # We extract information on the spacecraft
        m_i = self.__sc.mass
        max_thrust = self.__sc.thrust
        isp = self.__sc.isp
        veff = isp * G0
        # And on the leg
        throttles = [x[3 + 3 * i : 6 + 3 * i] for i in range(n_seg)]
        # Return lists
        n_points_fwd = fwd_seg + 1
        n_points_bwd = bwd_seg + 1
        rfwd = [[0.0] * 3] * (n_points_fwd)
        vfwd = [[0.0] * 3] * (n_points_fwd)
        mfwd = [0.0] * (n_points_fwd)
        ufwd = [[0.0] * 3] * (fwd_seg)
        rbwd = [[0.0] * 3] * (n_points_bwd)
        vbwd = [[0.0] * 3] * (n_points_bwd)
        mbwd = [0.0] * (n_points_bwd)
        ubwd = [[0.0] * 3] * (bwd_seg)

        # 2 - We compute the initial and final epochs and ephemerides
        t_i = epoch(t0)
        r_i, v_i = self.__earth.eph(t_i)
        if self.__start == 'l1':
            r_i = [r * 0.99 for r in r_i]
            v_i = [v * 0.99 for v in v_i]
        elif self.__start == 'l2':
            r_i = [r * 1.01 for r in r_i]
            v_i = [v * 1.01 for v in v_i]
        t_f = epoch(t0 + T)
        r_f, v_f = self.target.eph(t_f)

        # 3 - Forward propagation
        fwd_grid = t0 + T * self.__fwd_grid # days
        fwd_dt = T * self.__fwd_dt # seconds
        # Initial conditions
        rfwd[0] = r_i
        vfwd[0] = v_i
        mfwd[0] = m_i
        # Propagate
        for i, t in enumerate(throttles[:fwd_seg]):
            ufwd[i] = [max_thrust * thr for thr in t]
            if self.__earth_gravity:
                r_E, v_E = self.__earth.eph(epoch(fwd_grid[i]))
                r_delta = [a - b for a,b in zip(r_E, rfwd[i])]
                r3 = sum([r**2 for r in r_delta])**(3/2)
                ufwd[i] = [ui + mfwd[i] * MU_EARTH / r3 * ri for ui,ri in zip(ufwd[i], r_delta)]
            rfwd[i+1], vfwd[i+1], mfwd[i+1] = propagate_taylor(rfwd[i],vfwd[i],mfwd[i],ufwd[i],fwd_dt[i],MU_SUN,veff,-10,-10)

        # 4 - Backward propagation
        bwd_grid = t0 + T * self.__bwd_grid # days
        bwd_dt = T * self.__bwd_dt # seconds
        # Final conditions
        rbwd[-1] = r_f
        vbwd[-1] = v_f
        mbwd[-1] = m_f
        # Propagate
        for i, t in enumerate(throttles[-1:-bwd_seg - 1:-1]):
            ubwd[-i-1] = [max_thrust * thr for thr in t]
            if self.__earth_gravity:
                r_E, v_E = self.__earth.eph(epoch(bwd_grid[-i-1]))
                r_delta = [a - b for a,b in zip(r_E, rbwd[-i-1])]
                r3 = sum([r**2 for r in r_delta])**(3/2)
                ubwd[-i-1] = [ui + mbwd[-i-1] * MU_EARTH / r3 * ri for ui,ri in zip(ubwd[-i-1], r_delta)]
            rbwd[-i-2], vbwd[-i-2], mbwd[-i-2] = propagate_taylor(rbwd[-i-1],vbwd[-i-1],mbwd[-i-1],ubwd[-i-1],-bwd_dt[-i-1],MU_SUN,veff,-10,-10)

        return rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt

    # Visualizes the trajectory
    def plot_trajectory(self, x, units=AU, plot_segments=True, ax=None):
        """
        ax = prob.plot_trajectory(self, x, units=AU, plot_segments=True, ax=None)

        - x: encoded trajectory
        - units: the length unit to be used in the plot
        - plot_segments: when true plots also the segments boundaries
        - [out] ax: matplotlib axis where to plot

        Plots the trajectory
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from PyKEP.orbit_plots import plot_planet, plot_taylor

        # Creating the axis if necessary
        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.gca(projection='3d')
        else:
            axis = ax

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = x[1]
        isp = self.__sc.isp
        veff = isp * G0
        fwd_grid = t0 + T * self.__fwd_grid # days
        bwd_grid = t0 + T * self.__bwd_grid # days

        throttles = [x[3 + 3 * i : 6 + 3 * i] for i in range(n_seg)]
        alphas = [min(1., np.linalg.norm(t)) for t in throttles]

        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt = self._propagate(x)

        # Plotting the Sun, the Earth and the target
        axis.scatter([0], [0], [0], color='y')
        plot_planet(self.__earth, epoch(t0), units=units, legend=True, color=(0.7, 0.7, 1), ax=axis)
        plot_planet(self.target, epoch(t0 + T), units=units, legend=True, color=(0.7, 0.7, 1), ax=axis)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / units
        yfwd[0] = rfwd[0][1] / units
        zfwd[0] = rfwd[0][2] / units

        for i in range(fwd_seg):
            plot_taylor(rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[i], MU_SUN, veff, N=10, units=units, color=(alphas[i], 0, 1-alphas[i]), ax=axis)
            xfwd[i+1] = rfwd[i+1][0] / units
            yfwd[i+1] = rfwd[i+1][1] / units
            zfwd[i+1] = rfwd[i+1][2] / units
        if plot_segments:
            axis.scatter(xfwd[:-1], yfwd[:-1], zfwd[:-1], label='nodes', marker='o', s=1)

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / units
        ybwd[-1] = rbwd[-1][1] / units
        zbwd[-1] = rbwd[-1][2] / units

        for i in range(bwd_seg):
            plot_taylor(rbwd[-i-1], vbwd[-i-1], mbwd[-i-1], ubwd[-i-1], -bwd_dt[-i-1], MU_SUN, veff, N=10, units=units, color=(alphas[-i-1], 0, 1-alphas[-i-1]), ax=axis)
            xbwd[-i-2] = rbwd[-i-2][0] / units
            ybwd[-i-2] = rbwd[-i-2][1] / units
            zbwd[-i-2] = rbwd[-i-2][2] / units
        if plot_segments:
            axis.scatter(xbwd[1:], ybwd[1:], zbwd[1:], marker='o', s=1)

        if ax is None:  # show only if axis is not set
            plt.show()

        return axis

    def plot_dists_thrust(self, x, ax=None):
        """
        ax = prob.plot_dists_thrust(self, x, ax=None)

        - x: encoded trajectory
        - [out] ax: matplotlib axis where to plot

        Plots the distance of the spacecraft from the Earth/Sun and the thrust profile
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        # Creating the axis if necessary
        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.add_subplot(111)
        else:
            axis = ax

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = x[1]
        fwd_grid = t0 + T * self.__fwd_grid # days
        bwd_grid = t0 + T * self.__bwd_grid # days

        throttles = [x[3 + 3 * i : 6 + 3 * i] for i in range(n_seg)]
        thrusts = [min(1., np.linalg.norm(t)) for t in throttles]
        thrusts.append(thrusts[-1]) # duplicate the last for plotting

        dist_earth = [0.0] * (n_seg + 2) # distances spacecraft - Earth
        dist_sun = [0.0] * (n_seg + 2) # distances spacecraft - Sun
        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt = self._propagate(x)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / AU
        yfwd[0] = rfwd[0][1] / AU
        zfwd[0] = rfwd[0][2] / AU
        r_E = [ri/AU for ri in self.__earth.eph(epoch(fwd_grid[0]))[0]]
        dist_earth[0] = np.linalg.norm([r_E[0]-xfwd[0], r_E[1]-yfwd[0], r_E[2]-zfwd[0]])
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])

        for i in range(fwd_seg):
            xfwd[i+1] = rfwd[i+1][0] / AU
            yfwd[i+1] = rfwd[i+1][1] / AU
            zfwd[i+1] = rfwd[i+1][2] / AU
            r_E = [ri/AU for ri in self.__earth.eph(epoch(fwd_grid[i+1]))[0]]
            dist_earth[i+1] = np.linalg.norm([r_E[0]-xfwd[i+1], r_E[1]-yfwd[i+1], r_E[2]-zfwd[i+1]])
            dist_sun[i+1] = np.linalg.norm([xfwd[i+1], yfwd[i+1], zfwd[i+1]])

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / AU
        ybwd[-1] = rbwd[-1][1] / AU
        zbwd[-1] = rbwd[-1][2] / AU
        r_E = [ri/AU for ri in self.__earth.eph(epoch(bwd_grid[-1]))[0]]
        dist_earth[-1] = np.linalg.norm([r_E[0]-xbwd[-1], r_E[1]-ybwd[-1], r_E[2]-zbwd[-1]])
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])

        for i in range(bwd_seg):
            xbwd[-i-2] = rbwd[-i-2][0] / AU
            ybwd[-i-2] = rbwd[-i-2][1] / AU
            zbwd[-i-2] = rbwd[-i-2][2] / AU
            r_E = [ri/AU for ri in self.__earth.eph(epoch(bwd_grid[-i-2]))[0]]
            dist_earth[-i-2] = np.linalg.norm([r_E[0]-xbwd[-i-2], r_E[1]-ybwd[-i-2], r_E[2]-zbwd[-i-2]])
            dist_sun[-i-2] = np.linalg.norm([xbwd[-i-2], ybwd[-i-2], zbwd[-i-2]])

        axis.set_xlabel("t [mjd2000]")
        # Plot Earth distance
        axis.plot(times,dist_earth, c = 'b', label = "sc-Earth")
        # Plot Sun distance
        axis.plot(times,dist_sun, c = 'y', label = "sc-Sun")
        axis.set_ylabel("distance [AU]", color = 'k')
        axis.set_ylim(bottom = 0.)
        axis.tick_params('y', colors = 'k')
        axis.legend(loc=2)
        # draw threshold where Earth gravity equals 0.1*Tmax
        axis.axhline(y = np.sqrt(MU_EARTH * self.__sc.mass / (self.__sc.thrust * 0.1)) / AU, c = 'g', ls = ":", lw = 1)
        # Plot thrust profile
        axis = axis.twinx()
        axis.step(np.concatenate((fwd_grid, bwd_grid[1:])), thrusts, where = "post", c = 'r', linestyle = '--')
        axis.set_ylabel("T/Tmax",color='r')
        axis.tick_params('y',colors='r')
        axis.set_xlim([times[0],times[-1]])
        axis.set_ylim([0,1.2])

        if ax is None:  # show only if axis is not set
            plt.show()

        return axis

    def double_segments(self,x):
        """
        new_prob, new_x = prob.double_segments(self,x)

        - x: encoded trajectory

        Returns the decision vector encoding a low trust trajectory having double the number of segments with respect to x
        and a 'similar' throttle history. In case high fidelity is True, and x is a feasible trajectory, the returned decision vector
        also encodes a feasible trajectory that can be further optimized
        Returns the problem and the decision vector encoding a low-thrust trajectory having double the number of
        segments with respect to x and the same throttle history. If x is a feasible trajectory, the new chromosome is "almost"
        feasible, due to the new refined Earth gravity that is now different in the 2 halves of each segment).
        """
        new_x = list(x[:3])
        for i in range(self.__n_seg):
            new_x.extend(x[3 + 3 * i : 6 + 3 * i] * 2)

        new_prob = earth_gravity(
             target = self.target,
             n_seg = 2 * self.__n_seg,
             grid_type = self.__grid_type,
             t0 = [epoch(self.lb[0]), epoch(self.ub[0])],
             tof = [self.lb[1], self.ub[1]],
             m0 = self.__sc.mass,
             Tmax = self.__sc.thrust,
             Isp = self.__sc.isp,
             earth_gravity = self.__earth_gravity,
             start = self.__start
        )

        return new_prob, new_x

    def pretty(self, x):
        """
        prob.pretty(x)

        - x: encoded trajectory

        Prints human readable information on the trajectory represented by the decision vector x
        """
        n_seg = self.__n_seg
        m_i = self.__sc.mass
        t0 = x[0]
        T = x[1]
        m_f = x[2]
        thrusts = [np.linalg.norm(x[3 + 3 * i : 6 + 3 * i]) for i in range(n_seg)]

        tf = t0 + T
        mP = m_i - m_f
        deltaV = self.__sc.isp * G0 * np.log(m_i / m_f)

        dt = np.append(self.__fwd_dt, self.__bwd_dt) * T / DAY2SEC
        time_thrusts_on = sum(dt[i] for i in range(len(thrusts)) if thrusts[i] > 0.1)

        print("Departure:", epoch(t0), "(", t0, "mjd2000 )")
        print("Time of flight:", T, "days")
        print("Arrival:", epoch(tf), "(", tf, "mjd2000 )")
        print("Delta-v:", deltaV, "m/s")
        print("Propellant consumption:", mP, "kg")
        print("Thrust-on time:", time_thrusts_on, "days")
