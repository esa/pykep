import pykep as pk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class _direct_base(object):
    """Base class for direct trajectory optimisation problems with one only leg.

    All inheriting classes will adopt, ``plot_traj``, ``plot_control``, and ``get_traj``.
    """

    def __init__(self, mass=1000., thrust=0.3, isp=3000., nseg=10, mu=pk.MU_SUN, hf=False):

        # segements
        if isinstance(nseg, int):
            self.nseg = nseg
        else:
            raise TypeError("nseg must be supplied as int.")

        # spacecraft
        self.sc = pk.sims_flanagan.spacecraft(mass, thrust, isp)

        # leg
        self.leg = pk.sims_flanagan.leg()
        self.leg.set_spacecraft(self.sc)
        self.leg.set_mu(mu)
        self.leg.high_fidelity = hf

    def get_nobj(self):
        return 1

    def get_nec(self):
        return 7

    def get_nic(self):
        return self.nseg

    def plot_traj(self, z, units=pk.AU, N=20, axes=None):
        """This function plots the 3 dimensional spacecraft trajectory, given a solution chromosome.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - units (``float``, ``int``): units by which to scale the trajectory dimensions.
            - N (``int``): Number of points to be plotted along one arc.
        """

        # a call to the fitness on the chromosome z will change the class data member leg and set it
        # to represent the data in the chromosome z
        self.fitness(z)

        # creates a figure if needed
        if axes is None:
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        # plots a small Sun
        axes.scatter([0], [0], [0], color='y')

        # plots the leg
        pk.orbit_plots.plot_sf_leg(self.leg, units=units, N=20, ax=axes)

        # plots problem specifics
        self._plot_traj(z, axes, units)

        return axes

    def plot_control(self, z, mark="k.-", time=True, axes=None):
        """Plots the control profile of the trajectory, as a function of time.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - time (``bool``): If ``True``, x-axis is time in ``mjd2000``. If ``False``, x-axis is node index.

        """

        # data
        traj = self.get_traj(z)

        # time
        t = traj[:, 0] - traj[0,0]

        # throttle
        u = traj[:, 8]

        # figure
        if axes is None:
            plt.figure()
            axes = plt.gca()

        # with time
        if time:
            axes.plot(t, u, mark)
            plt.xlabel("Time [days]")

        # without time
        elif not time:
            axes.plot(u, mark)
            plt.xlabel("Segment number")
        # no se
        else:
            raise RuntimeError("Something is wrong!")

        # label
        plt.ylabel("Throttle [ND]")

        return axes

    def get_traj(self, z):
        """Retrieves the trajectory information.
        ::

            traj = [[t0, x0, y0, z0, vx0, vy0, vz0, m0, u0, ux0, uy0, uz0]
                    ...
                    [tf, xf, yf, zf, vxf, vyf, vzf, mf, uf, uxf, uyf, uzf]]

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
        """

        # set leg
        self.fitness(z)

        # get states
        x = list(self.leg.get_states())
        # clean states
        for i in range(len(x)):
            # remove matchpoint duplicate
            x[i].pop(self.nseg)
            # convert to numpy.ndarray
            x[i] = np.asarray(x[i], np.float64)
            # time and mass
            if i in [0, 3]:
                x[i].reshape((self.nseg * 2 + 1, 1))
            # position and velocity
            elif i in [1, 2]:
                x[i].reshape((self.nseg * 2 + 1, 3))
            else:
                raise RuntimeError("Something is wrong!")

        # unpack states
        t, r, v, m = x

        # control
        u = self._get_controls(z)
        # since controls are only defined at midpoints we need to add the values at the non midpoint nodes
        tmp = [0] * len(u)*2
        for i in range(self.nseg):
            tmp[i*6] = u[i*3]
            tmp[i*6+1] = u[i*3+1]
            tmp[i*6+2] = u[i*3+2]
            tmp[i*6+3] = u[i*3]
            tmp[i*6+4] = u[i*3+1]
            tmp[i*6+5] = u[i*3+2]
        tmp.append(u[-3])
        tmp.append(u[-2])
        tmp.append(u[-1])
        u = np.asarray(tmp, np.float64).reshape((self.nseg*2+1, 3))
        # throttle
        umag = np.linalg.norm(u, axis=1).reshape((self.nseg*2+1, 1))

        # full dataset [t, x, y, z, vx, vy, vz, m, u, ux, uy, uz]
        return np.hstack((t.reshape((self.nseg*2+1, 1)), r, v, m.reshape((self.nseg*2+1, 1)), umag, u))

    def pretty(self, z):
        data = self.get_traj(z) 
        self._pretty(z)

        print("\nSpacecraft Initial Position (m)  : [{!r}, {!r}, {!r}]".format(data[0,1], data[0,2], data[0,3]))
        print("Spacecraft Initial Velocity (m/s): [{!r}, {!r}, {!r}]".format(data[0,4], data[0,5], data[0,6]))
        print("Spacecraft Initial Mass  (kg)    : {!r}".format(data[0,7]))

        print("Spacecraft Final Position (m)  : [{!r}, {!r}, {!r}]".format(data[-1,1], data[-1,2], data[-1,3]))
        print("Spacecraft Final Velocity (m/s): [{!r}, {!r}, {!r}]".format(data[-1,4], data[-1,5], data[-1,6]))
        print("Spacecraft Final Mass  (kg)    : {!r}".format(data[-1,7]))
           


class direct_pl2pl(_direct_base):
    """Represents a direct transcription transfer between solar system planets.

    This problem works by manipulating the starting epoch t0, the transfer time T the final mass mf and the controls 
    ::

        z = [t0, T, mf, Vxi, Vyi, Vzi, Vxf, Vyf, Vzf, controls]
    """

    def __init__(self,
                 p0="earth",
                 pf="mars",
                 mass=1000,
                 thrust=0.3,
                 isp=3000,
                 nseg=20,
                 t0=[500, 1000],
                 tof=[200, 500],
                 vinf_dep=1e-3,
                 vinf_arr=1e-3,
                 hf=False):
        """Initialises a direct transcription orbit to orbit problem.

        Args:
            - p0 (``str``): Departure planet name. (will be used to construct a planet.jpl_lp object)
            - pf (``str``): Arrival planet name. (will be used to construct a planet.jpl_lp object)
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - nseg (``int``): Number of colocation nodes.
            - t0 (``list``): Launch epochs bounds [mjd2000].
            - tof (``list``): Transfer time bounds [days].
            - vinf_dep (``float``): allowed launch DV [km/s] 
            - vinf_arr (``float``): allowed arrival DV [km/s]
            - hf (``bool``): High-fidelity. Activates a continuous representation for the thrust.
        """

        # initialise base
        _direct_base.__init__(self, mass, thrust, isp, nseg, pk.MU_SUN, hf)

        # planets
        if all([isinstance(pl, str) for pl in [p0, pf]]):
            self.p0 = pk.planet.jpl_lp(p0)
            self.pf = pk.planet.jpl_lp(pf)
        else:
            raise TypeError("Planet names must be supplied as str.")

        # bounds TODO check
        self.t0 = t0
        self.tof = tof

        # boundary conditions on velocity
        self.vinf_dep = vinf_dep * 1000  # (in m)
        self.vinf_arr = vinf_arr * 1000  # (in m)

        # The class is built around solar system plaents hence mu is always the
        # SUN
        self.mu = pk.MU_SUN

    def fitness(self, z):

        # epochs (mjd2000)
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # final mass
        mf = z[2]

        # controls
        u = z[9:]

        # compute Cartesian states of planets
        r0, v0 = self.p0.eph(t0)
        rf, vf = self.pf.eph(tf)

        # add the vinfs from the chromosome
        v0 = [a + b for a, b in zip(v0, z[3:6])]
        vf = [a + b for a, b in zip(vf, z[6:9])]

        # spacecraft states
        x0 = pk.sims_flanagan.sc_state(r0, v0, self.sc.mass)
        xf = pk.sims_flanagan.sc_state(rf, vf, mf)

        # set leg
        self.leg.set(t0, x0, u, tf, xf)

        # compute equality constraints
        ceq = np.asarray(self.leg.mismatch_constraints(), np.float64)

        # nondimensionalise equality constraints
        ceq[0:3] /= pk.AU
        ceq[3:6] /= pk.EARTH_VELOCITY
        ceq[6] /= self.sc.mass

        # compute inequality constraints
        cineq = np.asarray(self.leg.throttles_constraints(), np.float64)

        # compute inequality constraints on departure and arrival velocities
        v_dep_con = (z[3] ** 2 + z[4] ** 2 + z[5] ** 2 - self.vinf_dep ** 2)
        v_arr_con = (z[6] ** 2 + z[7] ** 2 + z[8] ** 2 - self.vinf_arr ** 2)

        # nondimensionalize inequality constraints
        v_dep_con /= pk.EARTH_VELOCITY ** 2
        v_arr_con /= pk.EARTH_VELOCITY ** 2

        return np.hstack(([-mf], ceq, cineq, [v_dep_con, v_arr_con]))

    def get_nic(self):
        return super(direct_pl2pl, self).get_nic() + 2

    def get_bounds(self):
        lb = [self.t0[0], self.tof[0], self.sc.mass * 0.1] + \
            [-self.vinf_dep] * 3 + [-self.vinf_arr] * 3 + \
            [-1, -1, -1] * self.nseg
        ub = [self.t0[1], self.tof[1], self.sc.mass] + \
            [self.vinf_dep] * 3 + [self.vinf_arr] * 3 + \
            [1, 1, 1] * self.nseg
        return (lb, ub)

    def _plot_traj(self, z, axis, units):

        # times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # plot Keplerian
        pk.orbit_plots.plot_planet(
            self.p0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(
            self.pf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)

    def _pretty(self, z):
        print("\nLow-thrust NEP transfer from " + self.p0.name + " to " + self.pf.name)
        print("\nLaunch epoch: {!r} MJD2000, a.k.a. {!r}".format(z[0], pk.epoch(z[0])))
        print("Arrival epoch: {!r} MJD2000, a.k.a. {!r}".format(z[0]+z[1], pk.epoch(z[0]+z[1])))
        print("Time of flight (days): {!r} ".format(z[1]))
        print("\nLaunch DV (km/s) {!r} - [{!r},{!r},{!r}]".format(np.sqrt(z[3]**2+z[4]**2+z[5]**2) / 1000, z[3] / 1000, z[4] / 1000, z[5] / 1000))
        print("Arrival DV (km/s) {!r} - [{!r},{!r},{!r}]".format(np.sqrt(z[6]**2+z[7]**2+z[8]**2) / 1000, z[6] / 1000, z[7] / 1000, z[8] / 1000))

    @staticmethod
    def _get_controls(z):
        return z[9:]


class direct_or2or(_direct_base):
    """Represents a direct transcription transfer between orbits.

    This problem works by manipulating the time of flight ``T``, final mass ``mf`` and mean anomolies ``M0, Mf``.
    ::

        z = [T, mf, M0, Mf, controls]
    """

    def __init__(self, elem0, elemf, mass, thrust, isp, nseg, Tlb, Tub, E0lb, E0ub, Eflb, Efub, mu=pk.MU_SUN, hf=True):
        """Initialises a direct transcription orbit to orbit problem.

        Args:
            - elem0 (``list``, ``tuple``, ``numpy.ndarray``): Departure Keplerian elements. The eccentric anomlay will be manipulated.
            - elemf (``list``, ``tuple``, ``numpy.ndarray``): Arrival Keplerian elements. The eccentric anomlay will be manipulated.
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - nseg (``int``): Number of colocation nodes.
            - Tlb (``float``, ``int``): Minimimum time of flight [mjd2000].
            - Tub (``float``, ``int``): Maximum time of flight [mjd2000].
            - E0lb (``float``, ``int``): Minimum departure eccentric anomoly [rad].
            - E0ub (``float``, ``int``): Maximum departure eccentric anomoly [rad].
            - Eflb (``float``, ``int``): Minimum arrival eccentric anomoly [rad].
            - E0fb (``float``, ``int``): Maximum arrival eccentric anomoly [rad].
            - mu (``float``): Gravitational parameter of primary body [m^3/s^2].
            - hf (``bool``): ``True`` for continous thrust, ``False`` for impulsive thrust.
        """

        if all([(isinstance(elem, list) or isinstance(elem, tuple) or isinstance(elem, np.ndarray)) for elem in [elem0, elemf]]):
            elem0 = np.asarray(elem0, np.float64)
            elemf = np.asarray(elemf, np.float64)
        else:
            raise ValueError(
                "Both elem0 and elemf must be supplied as instances of list, tuple, or numpy.ndarray.")

        if all([elem.size == 6 for elem in [elem0, elemf]]):
            self.elem0 = elem0
            self.elemf = elemf
        else:
            raise TypeError("Both elem0 and elemf must be 6-dimensional.")

        if not all([(isinstance(T, float) or isinstance(T, int)) for T in [Tlb, Tub]]):
            raise TypeError(
                "Both Tlb and Tub must be supplied as instances of either float or int.")
        elif not Tlb < Tub:
            raise ValueError("Tlb must be less than Tub.")
        else:
            self.Tlb = float(Tlb)
            self.Tub = float(Tub)

        # initialise base
        _direct_base.__init__(self, mass, thrust, isp, nseg, mu, hf)

        self.mu = mu
        self.E0lb = E0lb
        self.E0ub = E0ub
        self.Eflb = Eflb
        self.Efub = Efub

    def fitness(self, z):

        # epochs (mjd2000)
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # final mass
        mf = z[1]

        # eccentric anomolies
        E0 = z[2]
        Ef = z[3]

        # controls
        u = z[4:]

        # set Keplerian elements
        self.elem0[5] = E0
        self.elemf[5] = Ef

        # compute Cartesian states
        r0, v0 = pk.par2ic(self.elem0, self.mu)
        rf, vf = pk.par2ic(self.elemf, self.mu)

        # spacecraft states
        x0 = pk.sims_flanagan.sc_state(r0, v0, self.sc.mass)
        xf = pk.sims_flanagan.sc_state(rf, vf, mf)

        # set leg
        self.leg.set(t0, x0, u, tf, xf)

        # compute equality constraints
        ceq = np.asarray(self.leg.mismatch_constraints(), np.float64)

        # nondimensionalise equality constraints
        ceq[0:3] /= pk.AU
        ceq[3:6] /= pk.EARTH_VELOCITY
        ceq[6] /= self.sc.mass

        # compute inequality constraints
        cineq = np.asarray(self.leg.throttles_constraints(), np.float64)

        return np.hstack(([-mf], ceq, cineq))

    def get_bounds(self):
        pi = 3.14159265359
        lb = [self.Tlb, self.sc.mass / 10, self.E0lb,
              self.Eflb, *(-1, -1, -1) * self.nseg]
        ub = [self.Tub, self.sc.mass, self.E0ub,
              self.Efub, *(1, 1, 1) * self.nseg]
        return (lb, ub)

    def _plot_traj(self, z, axis, units):

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # Keplerian
        kep0 = pk.planet.keplerian(t0, self.elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # plot Keplerian
        pk.orbit_plots.plot_planet(
            kep0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(
            kepf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)

    @staticmethod
    def _get_controls(z):
        return z[4:]
