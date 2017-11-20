import pykep as pk
import numpy as np
import math

mpcorbline = "K14Y00D 24.3   0.15 K1794 105.44160   34.12337  117.64264    1.73560  0.0865962  0.88781021   1.0721510  2 MPO369254   104   1  194 days 0.24 M-v 3Eh MPCALB     2803          2014 YD            20150618"


class lt_margo:
    """
    This class can be used as a User Defined Problem (UDP) in the pygmo2 software and if successfully solved,
    represents a low-thrust interplanetary trajectory from the Earth (or from the
    Sun-Earth L1 or L2 Lagrangian point) to a target NEO. The trajectory is modeled using the Sims-Flanagan model, 
    extended to include the Earth's gravity (assumed constant along each segment).
    The propulsion model can be both nuclear (NEP) or solar (SEP).

    This problem was developed and used during the European Space Agency interplanetary cubesat M-ARGO study.

    The decision vector (chromosome) is::

      [t0, tof, mf] + [throttles1] + [throttles2] + ...
    """

    def __init__(self,
                 target=pk.planet.mpcorb(mpcorbline),
                 n_seg=30,
                 grid_type="uniform",
                 t0=[pk.epoch(8000), pk.epoch(9000)],
                 tof=[200, 365.25 * 3],
                 m0=20.0,
                 Tmax=0.0017,
                 Isp=3000.0,
                 earth_gravity=False,
                 sep=False,
                 start="earth"):
        """
        prob = lt_margo(target = pk.planet.mpcorb(mpcorbline), n_seg = 30, grid_type = "uniform", t0 = [epoch(8000), epoch(9000)], tof = [200, 365.25*3], m0 = 20.0, Tmax = 0.0017, Isp = 3000.0, earth_gravity = False, sep = False, start = "earth")

        Args:
            - target (``pykep.planet``): target planet
            - n_seg (``int``): number of segments to use in the problem transcription (time grid)
            - grid_type (``string``): "uniform" for uniform segments, "nonuniform" toi use a denser grid in the first part of the trajectory
            - t0 (``list``): list of two pykep.epoch defining the bounds on the launch epoch
            - tof (``list``): list of two floats defining the bounds on the time of flight (days)
            - m0 (``float``): initial mass of the spacecraft (kg)
            - Tmax (``float``): maximum thrust at 1 AU (N)
            - Isp (``float``): engine specific impulse (s)
            - earth_gravity (``bool``): activates the Earth gravity effect in the dynamics
            - sep (``bool``): Activates a Solar Electric Propulsion model for the thrust - distance dependency.
            - start(``string``): Starting point ("earth", "l1", or "l2").

        .. note::

           L1 and L2 are approximated as the points on the line connecting the Sun and the Earth at a distance of, respectively, 0.99 and 1.01 AU from the Sun.

        .. note::

           If the Earth's gravity is enabled, the starting point cannot be the Earth
        """

        # 1) Various checks
        if start not in ["earth", "l1", "l2"]:
            raise ValueError("start must be either 'earth', 'l1' or 'l2'")
        if grid_type not in ["uniform", "nonuniform"]:
            raise ValueError(
                "grid_type must be either 'uniform' or 'nonuniform'")
        if earth_gravity and start == "earth":
            raise ValueError(
                "If Earth gravity is enabled the starting point cannot be the Earth")

        # 2) Class data members
        # public:
        self.target = target
        # private:
        self.__n_seg = n_seg
        self.__grid_type = grid_type
        self.__sc = pk.sims_flanagan._sims_flanagan.spacecraft(m0, Tmax, Isp)
        self.__earth = pk.planet.jpl_lp('earth')
        self.__earth_gravity = earth_gravity
        self.__sep = sep
        self.__start = start
        # grid construction
        if grid_type == "uniform":
            grid = np.array([i / n_seg for i in range(n_seg + 1)])
        elif grid_type == "nonuniform":
            grid_f = lambda x: x**2 if x < 0.5 else 0.25 + 1.5 * \
                (x - 0.5)  # quadratic in [0,0.5], linear in [0.5,1]
            grid = np.array([grid_f(i / n_seg) for i in range(n_seg + 1)])
        # index corresponding to the middle of the transfer
        fwd_seg = int(np.searchsorted(grid, 0.5, side='right'))
        bwd_seg = n_seg - fwd_seg
        fwd_grid = grid[:fwd_seg + 1]
        bwd_grid = grid[fwd_seg:]
        self.__fwd_seg = fwd_seg
        self.__fwd_grid = fwd_grid
        self.__fwd_dt = np.array([(fwd_grid[i + 1] - fwd_grid[i])
                                  for i in range(fwd_seg)]) * pk.DAY2SEC
        self.__bwd_seg = bwd_seg
        self.__bwd_grid = bwd_grid
        self.__bwd_dt = np.array([(bwd_grid[i + 1] - bwd_grid[i])
                                  for i in range(bwd_seg)]) * pk.DAY2SEC

        # 3) Bounds
        lb = [t0[0].mjd2000] + [tof[0]] + [0] + [-1, -1, -1] * n_seg
        ub = [t0[1].mjd2000] + [tof[1]] + [m0] + [1, 1, 1] * n_seg
        self.__lb = lb
        self.__ub = ub

    # Fitness function
    def fitness(self, x):
        obj = -x[2]

        ceq = list()
        cineq = list()
        throttles_constraints = []
        mismatch_constraints = []

        throttles = [x[3 + 3 * i: 6 + 3 * i] for i in range(self.__n_seg)]
        for t in throttles:
            throttles_constraints.append(t[0]**2 + t[1]**2 + t[2]**2 - 1.)
        cineq.extend(throttles_constraints)

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, _, _, _, _, dfwd, dbwd = self._propagate(
            x)

        if self.__earth_gravity:
            flyby_constraints = []
            for d in dfwd + dbwd:
                d2 = d[0]**2 + d[1]**2 + d[2]**2
                l22 = (0.01 * pk.AU)**2
                flyby_constraints.append(1. - d2 / l22)
            cineq.extend(flyby_constraints)

        mismatch_constraints.extend([a - b for a, b in zip(rfwd[-1], rbwd[0])])
        mismatch_constraints.extend([a - b for a, b in zip(vfwd[-1], vbwd[0])])
        mismatch_constraints.append(mfwd[-1] - mbwd[0])
        ceq.extend(mismatch_constraints)

        # Making the mismatches non dimensional
        ceq[0] /= pk.AU
        ceq[1] /= pk.AU
        ceq[2] /= pk.AU
        ceq[3] /= pk.EARTH_VELOCITY
        ceq[4] /= pk.EARTH_VELOCITY
        ceq[5] /= pk.EARTH_VELOCITY
        ceq[6] /= self.__sc.mass

        # We assemble the constraint vector
        retval = [obj]
        retval.extend(ceq)
        retval.extend(cineq)

        return retval

    # Get bounds
    def get_bounds(self):
        return (self.__lb, self.__ub)

    # Get number of inequality contraints
    def get_nic(self):
        if self.__earth_gravity:
            return 2 * self.__n_seg
        else:
            return self.__n_seg

    # Get number of equality contraints
    def get_nec(self):
        return 7

    # Gradient sparsity
    def gradient_sparsity(self):
        # Objective
        retval = [[0, 2]]
        # Mismatches
        retval += [[i, j] for i in range(1, 7)
                   for j in range(3 * (1 + self.__n_seg))]
        retval += [[7, j] for j in range(1, 3 * (1 + self.__n_seg))]
        # Throttles
        retval += [[8 + i, 3 + 3 * i + j]
                   for i in range(self.__n_seg) for j in range(3)]
        # Prevent flyby
        if self.__earth_gravity:
            retval += [[8 + self.__n_seg + i, j]
                       for i in range(self.__fwd_seg) for j in [0, 1] + list(range(3, 6 + 3 * i))]
            retval += [[8 + self.__n_seg + self.__fwd_seg + i, j] for i in range(self.__bwd_seg) for j in [
                0, 1, 2] + list(range(3 + 3 * (self.__n_seg - self.__bwd_seg + i), 3 + 3 * self.__n_seg))]
        return retval

    def get_name(self):
        return "MARGO cubesat transfer to " + self.target.name

    def get_extra_info(self):
        retval = "\tTarget planet: " + self.target.name
        retval += "\n\tEarth gravity: " + str(self.__earth_gravity)
        retval += "\n\tSolar Electric Propulsion: " + str(self.__sep)
        retval += "\n\tStart mass: " + str(self.__sc.mass) + " kg"
        retval += "\n\tMaximum thrust as 1AU: " + str(self.__sc.thrust) + " N"
        retval += "\n\tSpecific impulse: " + str(self.__sc.isp) + " s"
        retval += "\n\n\tLaunch window: [" + \
            str(self.__lb[0]) + ", " + str(self.__ub[0]) + "] - MJD2000"
        retval += "\n\tBounds on time of flight: [" + str(
            self.__lb[1]) + ", " + str(self.__ub[1]) + "] - days"
        retval += "\n\n\tNumber of segments: " + str(self.__n_seg)
        retval += "\n\tGrid type: " + self.__grid_type

        return retval

    # SEP model
    def _sep_model(self, r):
        SAA = 0
        Psupply = 13.75
        eff = 0.92

        Pbmp = (-40.558 * r**3 + 173.49 * r**2 -
                259.19 * r + 141.86) * math.cos(SAA)
        P = -146.26 * r**3 + 658.52 * r**2 - 1059.2 * r + 648.24  # 6 panels
        # P = -195.02 * r**3 + 878.03 * r**2 - 1412.3 * r + 864.32 # 8 panels
        if Pbmp < Psupply:
            P -= (Psupply - Pbmp)
        Pin = eff * P
        if Pin > 120:
            Pin = 120  # thermal max 120W

        Tmax = (26.27127 * Pin - 708.973) / 1000000
        if Tmax < 0:
            Tmax = 0

        Isp = -0.0011 * Pin**3 + 0.175971 * Pin**2 + 4.193797 * Pin + 2037.213

        return Tmax, Isp

    # Propagates the trajectory
    def _propagate(self, x):
        # 1 - We decode the chromosome
        t0 = x[0]
        T = x[1]
        m_f = x[2]
        # We extract the number of segments for forward and backward
        # propagation
        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        # We extract information on the spacecraft
        m_i = self.__sc.mass
        max_thrust = self.__sc.thrust
        isp = self.__sc.isp
        veff = isp * pk.G0
        # And on the leg
        throttles = [x[3 + 3 * i: 6 + 3 * i] for i in range(n_seg)]
        # Return lists
        n_points_fwd = fwd_seg + 1
        n_points_bwd = bwd_seg + 1
        rfwd = [[0.0] * 3] * (n_points_fwd)
        vfwd = [[0.0] * 3] * (n_points_fwd)
        mfwd = [0.0] * (n_points_fwd)
        ufwd = [[0.0] * 3] * (fwd_seg)
        dfwd = [[0.0] * 3] * (fwd_seg)  # distances E/S
        rbwd = [[0.0] * 3] * (n_points_bwd)
        vbwd = [[0.0] * 3] * (n_points_bwd)
        mbwd = [0.0] * (n_points_bwd)
        ubwd = [[0.0] * 3] * (bwd_seg)
        dbwd = [[0.0] * 3] * (bwd_seg)

        # 2 - We compute the initial and final epochs and ephemerides
        t_i = pk.epoch(t0)
        r_i, v_i = self.__earth.eph(t_i)
        if self.__start == 'l1':
            r_i = [r * 0.99 for r in r_i]
            v_i = [v * 0.99 for v in v_i]
        elif self.__start == 'l2':
            r_i = [r * 1.01 for r in r_i]
            v_i = [v * 1.01 for v in v_i]
        t_f = pk.epoch(t0 + T)
        r_f, v_f = self.target.eph(t_f)

        # 3 - Forward propagation
        fwd_grid = t0 + T * self.__fwd_grid  # days
        fwd_dt = T * self.__fwd_dt  # seconds
        # Initial conditions
        rfwd[0] = r_i
        vfwd[0] = v_i
        mfwd[0] = m_i
        # Propagate
        for i, t in enumerate(throttles[:fwd_seg]):
            if self.__sep:
                r = math.sqrt(rfwd[i][0]**2 + rfwd[i][1]
                              ** 2 + rfwd[i][2]**2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ufwd[i] = [max_thrust * thr for thr in t]
            if self.__earth_gravity:
                r_E, v_E = self.__earth.eph(pk.epoch(fwd_grid[i]))
                dfwd[i] = [a - b for a, b in zip(r_E, rfwd[i])]
                r3 = sum([r**2 for r in dfwd[i]])**(3 / 2)
                disturbance = [mfwd[i] * pk.MU_EARTH /
                               r3 * ri for ri in dfwd[i]]
                rfwd[i + 1], vfwd[i + 1], mfwd[i + 1] = pk.propagate_taylor_disturbance(
                    rfwd[i], vfwd[i], mfwd[i], ufwd[i], disturbance, fwd_dt[i], pk.MU_SUN, veff, -10, -10)
            else:
                rfwd[i + 1], vfwd[i + 1], mfwd[i + 1] = pk.propagate_taylor(
                    rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[i], pk.MU_SUN, veff, -10, -10)

        # 4 - Backward propagation
        bwd_grid = t0 + T * self.__bwd_grid  # days
        bwd_dt = T * self.__bwd_dt  # seconds
        # Final conditions
        rbwd[-1] = r_f
        vbwd[-1] = v_f
        mbwd[-1] = m_f
        # Propagate
        for i, t in enumerate(throttles[-1:-bwd_seg - 1:-1]):
            if self.__sep:
                r = math.sqrt(rbwd[-i - 1][0]**2 + rbwd[-i - 1]
                              [1]**2 + rbwd[-i - 1][2]**2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ubwd[-i - 1] = [max_thrust * thr for thr in t]
            if self.__earth_gravity:
                r_E, v_E = self.__earth.eph(pk.epoch(bwd_grid[-i - 1]))
                dbwd[-i - 1] = [a - b for a, b in zip(r_E, rbwd[-i - 1])]
                r3 = sum([r**2 for r in dbwd[-i - 1]])**(3 / 2)
                disturbance = [mfwd[i] * pk.MU_EARTH /
                               r3 * ri for ri in dbwd[-i - 1]]
                rbwd[-i - 2], vbwd[-i - 2], mbwd[-i - 2] = pk.propagate_taylor_disturbance(
                    rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], disturbance, -bwd_dt[-i - 1], pk.MU_SUN, veff, -10, -10)
            else:
                rbwd[-i - 2], vbwd[-i - 2], mbwd[-i - 2] = pk.propagate_taylor(
                    rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1], pk.MU_SUN, veff, -10, -10)

        return rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd

    # Visualizes the trajectory
    def plot_traj(self, x, units=pk.AU, plot_segments=True, plot_thrusts=False, axes=None):
        """
        ax = prob.plot_traj(self, x, units=AU, plot_segments=True, plot_thrusts=False, axes=None)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - units (``float``): the length unit to be used in the plot
            - plot_segments (``bool``): when True plots also the segments boundaries
            - plot_thrusts (``bool``): when True plots also the thrust vectors

        Returns:
            matplotlib.axes: axes where to plot

        Visualizes the the trajectory in a 3D plot
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = x[1]
        isp = self.__sc.isp
        veff = isp * pk.G0
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [x[3 + 3 * i: 6 + 3 * i] for i in range(n_seg)]
        alphas = [min(1., np.linalg.norm(t)) for t in throttles]

        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd = self._propagate(
            x)

        # Plotting the Sun, the Earth and the target
        axes.scatter([0], [0], [0], color='y')
        pk.orbit_plots.plot_planet(self.__earth, pk.epoch(
            t0), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)
        pk.orbit_plots.plot_planet(self.target, pk.epoch(
            t0 + T), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / units
        yfwd[0] = rfwd[0][1] / units
        zfwd[0] = rfwd[0][2] / units

        for i in range(fwd_seg):
            if self.__sep:
                r = math.sqrt(rfwd[i][0]**2 + rfwd[i][1]
                              ** 2 + rfwd[i][2]**2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            if self.__earth_gravity:
                r3 = sum([r**2 for r in dfwd[i]])**(3 / 2)
                disturbance = [mfwd[i] * pk.MU_EARTH /
                               r3 * ri for ri in dfwd[i]]
                pk.orbit_plots.plot_taylor_disturbance(rfwd[i], vfwd[i], mfwd[i], ufwd[i], disturbance, fwd_dt[
                                                       i], pk.MU_SUN, veff, N=10, units=units, color=(alphas[i], 0, 1 - alphas[i]), ax=axes)
            else:
                pk.orbit_plots.plot_taylor(rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[
                                           i], pk.MU_SUN, veff, N=10, units=units, color=(alphas[i], 0, 1 - alphas[i]), ax=axes)
            xfwd[i + 1] = rfwd[i + 1][0] / units
            yfwd[i + 1] = rfwd[i + 1][1] / units
            zfwd[i + 1] = rfwd[i + 1][2] / units
        if plot_segments:
            axes.scatter(xfwd[:-1], yfwd[:-1], zfwd[:-1],
                         label='nodes', marker='o', s=5, c='k')

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / units
        ybwd[-1] = rbwd[-1][1] / units
        zbwd[-1] = rbwd[-1][2] / units

        for i in range(bwd_seg):
            if self.__sep:
                r = math.sqrt(rbwd[-i - 1][0]**2 + rbwd[-i - 1]
                              [1]**2 + rbwd[-i - 1][2]**2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            if self.__earth_gravity:
                r3 = sum([r**2 for r in dbwd[-i - 1]])**(3 / 2)
                disturbance = [mfwd[i] * pk.MU_EARTH /
                               r3 * ri for ri in dbwd[-i - 1]]
                pk.orbit_plots.plot_taylor_disturbance(rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], disturbance, -bwd_dt[
                                                       -i - 1], pk.MU_SUN, veff, N=10, units=units, color=(alphas[-i - 1], 0, 1 - alphas[-i - 1]), ax=axes)
            else:
                pk.orbit_plots.plot_taylor(rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1],
                                           pk.MU_SUN, veff, N=10, units=units, color=(alphas[-i - 1], 0, 1 - alphas[-i - 1]), ax=axes)
            xbwd[-i - 2] = rbwd[-i - 2][0] / units
            ybwd[-i - 2] = rbwd[-i - 2][1] / units
            zbwd[-i - 2] = rbwd[-i - 2][2] / units
        if plot_segments:
            axes.scatter(xbwd[1:], ybwd[1:], zbwd[1:], marker='o', s=5, c='k')

        # Plotting the thrust vectors
        if plot_thrusts:
            throttles = np.array(throttles)
            xlim = axes.get_xlim()
            xrange = xlim[1] - xlim[0]
            ylim = axes.get_ylim()
            yrange = ylim[1] - ylim[0]
            zlim = axes.get_zlim()
            zrange = zlim[1] - zlim[0]

            scale = 0.1

            throttles[:, 0] *= xrange
            throttles[:, 1] *= yrange
            throttles[:, 2] *= zrange

            throttles *= scale

            for (x, y, z, t) in zip(xfwd[:-1] + xbwd[:-1], yfwd[:-1] + ybwd[:-1], zfwd[:-1] + zbwd[:-1], throttles):
                axes.plot([x, x + t[0]], [y, y + t[1]], [z, z + t[2]], c='g')

        return axes

    def plot_dists_thrust(self, x, axes=None):
        """
        axes = prob.plot_dists_thrust(self, x, axes=None)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Returns:
            matplotlib.axes: axes where to plot

        Plots the distance of the spacecraft from the Earth/Sun and the thrust profile
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.add_subplot(111)

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = x[1]
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [np.linalg.norm(x[3 + 3 * i: 6 + 3 * i])
                     for i in range(n_seg)]

        dist_earth = [0.0] * (n_seg + 2)  # distances spacecraft - Earth
        dist_sun = [0.0] * (n_seg + 2)  # distances spacecraft - Sun
        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, _, _ = self._propagate(
            x)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / pk.AU
        yfwd[0] = rfwd[0][1] / pk.AU
        zfwd[0] = rfwd[0][2] / pk.AU
        r_E = [ri / pk.AU for ri in self.__earth.eph(pk.epoch(fwd_grid[0]))[0]]
        dist_earth[0] = np.linalg.norm(
            [r_E[0] - xfwd[0], r_E[1] - yfwd[0], r_E[2] - zfwd[0]])
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])

        for i in range(fwd_seg):
            xfwd[i + 1] = rfwd[i + 1][0] / pk.AU
            yfwd[i + 1] = rfwd[i + 1][1] / pk.AU
            zfwd[i + 1] = rfwd[i + 1][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.__earth.eph(pk.epoch(fwd_grid[i + 1]))[0]]
            dist_earth[
                i + 1] = np.linalg.norm([r_E[0] - xfwd[i + 1], r_E[1] - yfwd[i + 1], r_E[2] - zfwd[i + 1]])
            dist_sun[
                i + 1] = np.linalg.norm([xfwd[i + 1], yfwd[i + 1], zfwd[i + 1]])

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / pk.AU
        ybwd[-1] = rbwd[-1][1] / pk.AU
        zbwd[-1] = rbwd[-1][2] / pk.AU
        r_E = [
            ri / pk.AU for ri in self.__earth.eph(pk.epoch(bwd_grid[-1]))[0]]
        dist_earth[-1] = np.linalg.norm([r_E[0] - xbwd[-1],
                                         r_E[1] - ybwd[-1], r_E[2] - zbwd[-1]])
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])

        for i in range(bwd_seg):
            xbwd[-i - 2] = rbwd[-i - 2][0] / pk.AU
            ybwd[-i - 2] = rbwd[-i - 2][1] / pk.AU
            zbwd[-i - 2] = rbwd[-i - 2][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.__earth.eph(pk.epoch(bwd_grid[-i - 2]))[0]]
            dist_earth[-i - 2] = np.linalg.norm(
                [r_E[0] - xbwd[-i - 2], r_E[1] - ybwd[-i - 2], r_E[2] - zbwd[-i - 2]])
            dist_sun[-i -
                     2] = np.linalg.norm([xbwd[-i - 2], ybwd[-i - 2], zbwd[-i - 2]])

        axes.set_xlabel("t [mjd2000]")
        # Plot Earth distance
        axes.plot(times, dist_earth, c='b', label="sc-Earth")
        # Plot Sun distance
        axes.plot(times, dist_sun, c='y', label="sc-Sun")
        axes.set_ylabel("distance [AU]", color='k')
        axes.set_ylim(bottom=0.)
        axes.tick_params('y', colors='k')
        axes.legend(loc=2)
        # draw threshold where Earth gravity equals 0.1*Tmax
        if self.__earth_gravity:
            axes.axhline(y=np.sqrt(pk.MU_EARTH * self.__sc.mass / (self.__sc.thrust * 0.1)
                                   ) / pk.AU, c='g', ls=":", lw=1, label="earth_g = 0.1*Tmax")
        # Plot thrust profile
        axes = axes.twinx()
        if self.__sep:
            max_thrust = self.__sc.thrust
            thrusts = np.linalg.norm(
                np.array(ufwd + ubwd), axis=1) / max_thrust
            # plot maximum thrust achievable at that distance from the Sun
            distsSun = dist_sun[:fwd_seg] + \
                dist_sun[-bwd_seg:] + [dist_sun[-1]]
            Tmaxs = [self._sep_model(d)[0] / max_thrust for d in distsSun]
            axes.step(np.concatenate(
                (fwd_grid, bwd_grid[1:])), Tmaxs, where="post", c='lightgray', linestyle=':')
        else:
            thrusts = throttles.copy()
        # duplicate the last for plotting
        thrusts = np.append(thrusts, thrusts[-1])
        axes.step(np.concatenate(
            (fwd_grid, bwd_grid[1:])), thrusts, where="post", c='r', linestyle='--')
        axes.set_ylabel("T/Tmax$_{1AU}$", color='r')
        axes.tick_params('y', colors='r')
        axes.set_xlim([times[0], times[-1]])
        axes.set_ylim([0, max(thrusts) + 0.2])

        return axes

    def double_segments(self, x):
        """
        new_prob, new_x = prob.double_segments(self,x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Returns:
            - the new udp having twice the segments
            - list: the new chromosome to be used as initial guess

        Returns the decision vector encoding a low trust trajectory having double the number of segments with respect to x
        and a 'similar' throttle history. In case high fidelity is True, and x is a feasible trajectory, the returned decision vector
        also encodes a feasible trajectory that can be further optimized
        Returns the problem and the decision vector encoding a low-thrust trajectory having double the number of
        segments with respect to x and the same throttle history. If x is a feasible trajectory, the new chromosome is "slightly
        unfeasible", due to the new refined Earth's gravity and/or SEP thrust that are now computed in the 2 halves of each segment.
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        new_x = np.append(x[:3], np.repeat(x[3:].reshape((-1, 3)), 2, axis=0))

        new_prob = lt_margo(
            target=self.target,
            n_seg=2 * self.__n_seg,
            grid_type=self.__grid_type,
            t0=[pk.epoch(self.__lb[0]), pk.epoch(self.__ub[0])],
            tof=[self.__lb[1], self.__ub[1]],
            m0=self.__sc.mass,
            Tmax=self.__sc.thrust,
            Isp=self.__sc.isp,
            earth_gravity=self.__earth_gravity,
            sep=self.__sep,
            start=self.__start
        )

        return new_prob, new_x

    def pretty(self, x):
        """
        prob.pretty(x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Prints human readable information on the trajectory represented by the decision vector x
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        n_seg = self.__n_seg
        m_i = self.__sc.mass
        t0 = x[0]
        T = x[1]
        m_f = x[2]
        thrusts = [np.linalg.norm(x[3 + 3 * i: 6 + 3 * i])
                   for i in range(n_seg)]

        tf = t0 + T
        mP = m_i - m_f
        deltaV = self.__sc.isp * pk.G0 * np.log(m_i / m_f)

        dt = np.append(self.__fwd_dt, self.__bwd_dt) * T / pk.DAY2SEC
        time_thrusts_on = sum(dt[i] for i in range(
            len(thrusts)) if thrusts[i] > 0.1)

        print("Departure:", pk.epoch(t0), "(", t0, "mjd2000 )")
        print("Time of flight:", T, "days")
        print("Arrival:", pk.epoch(tf), "(", tf, "mjd2000 )")
        print("Delta-v:", deltaV, "m/s")
        print("Propellant consumption:", mP, "kg")
        print("Thrust-on time:", time_thrusts_on, "days")
