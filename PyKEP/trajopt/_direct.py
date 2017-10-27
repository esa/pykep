import PyKEP as pk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class _direct_base(object):
    """Base class for direct trajectory optimisation problems.

    All inheriting classes will adopt, ``plot_traj``, ``plot_control``, and ``get_traj``.
    """

    def __init__(self, mass, thrust, isp, nseg, mu=pk.MU_SUN, hf=True):

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

    def plot_traj(self, z, units=pk.AU, N=20):
        """This function plots the 3 dimensional spacecraft trajectory, given a solution chromosome.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - units (``float``, ``int``): units by which to scale the trajectory dimensions.
            - N (``int``): Number of points to be plotted along one arc.
        """

        # set problem
        self.fitness(z)

        # figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        pk.orbit_plots.plot_sf_leg(self.leg, units=units, N=20, ax=axis)

        # plot problem specifics
        self._plot_traj(z, axis, units)

        # show plot
        plt.show()

    def plot_control(self, z, mark="k.-", time=True):
        """Plots the control profile of the trajectory, as a function of time.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - time (``bool``): If ``True``, x-axis is time in ``mjd2000``. If ``False``, x-axis is node index.

        """

        # data
        traj = self.get_traj(z)

        # time
        t = traj[:, 0]

        # throttle
        u = traj[:, 8]

        # fiure
        plt.figure()

        # with time
        if time:
            plt.plot(t, u, mark)
            plt.xlabel("Time [days]")

        # without time
        elif not time:
            plt.plot(u, mark)
            plt.xlabel("Segement number")

        # no se
        else:
            raise RuntimeError("Något är fel!")

        # label
        plt.ylabel("Throttle [ND]")

        # show
        plt.show()

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
            # keep only midpoints
            x[i] = x[i][1:self.nseg*2:2]
            # convert to numpy.ndarray
            x[i] = np.asarray(x[i], np.float64)
            # time and mass
            if i in [0, 3]:
                x[i].reshape((self.nseg, 1))
            # position and velocity
            elif i in [1, 2]:
                x[i].reshape((self.nseg, 3))
            else:
                raise RuntimeError("Something is worng!")

        # unpack states
        t, r, v, m = x

        # control
        u = np.asarray(self._get_controls(z), np.float64).reshape((self.nseg, 3))
        # throttle
        umag = np.linalg.norm(u, axis=1).reshape((self.nseg, 1))

        # full dataset [t, x, y, z, vx, vy, vz, m, u, ux, uy, uz]
        return np.hstack((t.reshape((self.nseg, 1)), r, v, m.reshape((self.nseg, 1)), umag, u))

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
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].
            - hf (``bool``): ``True`` for continous thrust, ``False`` for impulsive thrust.
        """

        if all([(isinstance(elem, list) or isinstance(elem, tuple) or isinstance(elem, np.ndarray)) for elem in [elem0, elemf]]):
            elem0 = np.asarray(elem0, np.float64)
            elemf = np.asarray(elemf, np.float64)
        else:
            raise ValueError("Both elem0 and elemf must be supplied as instances of list, tuple, or numpy.ndarray.")

        if all([elem.size == 6 for elem in [elem0, elemf]]):
            self.elem0 = elem0
            self.elemf = elemf
        else:
            raise TypeError("Both elem0 and elemf must be 6-dimensional.")

        if not all([(isinstance(T, float) or isinstance(T, int)) for T in [Tlb, Tub]]):
            raise TypeError("Both Tlb and Tub must be supplied as instances of either float or int.")
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
        lb = [self.Tlb, self.sc.mass/10, self.E0lb, self.Eflb, *(-1, -1, -1)*self.nseg]
        ub = [self.Tub, self.sc.mass, self.E0ub, self.Efub, *(1, 1, 1)*self.nseg]
        return (lb, ub)

    def _plot_traj(self, z, axis, units):

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # Keplerian
        kep0 = pk.planet.keplerian(t0, self.elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # plot Keplerian
        pk.orbit_plots.plot_planet(kep0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(kepf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)

    @staticmethod
    def _get_controls(z):
        return z[4:]
