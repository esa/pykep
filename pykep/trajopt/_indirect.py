import pykep as pk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class _indirect_base(object):
    """Base class for indirect trajectory optimisation problems.

    All inheriting classes will adopt ``plot_traj`` and ``plot_control`` and ``pretty`` and ``get_traj`` and have
    to implement the methods ``_plot_traj`` and ``_pretty``: 
    """

    def __init__(
        self, mass, thrust, isp,
        mu, freemass, freetime, alpha, bound,
        atol, rtol
    ):

        # spacecraft
        self.sc = pk.sims_flanagan.spacecraft(mass, thrust, isp)

        # indirect leg
        self.leg = pk.pontryagin.leg(
            sc=self.sc, mu=mu, freemass=freemass, freetime=freetime, alpha=alpha, bound=bound
        )

        # integration parameters
        if all([(isinstance(par, float) or isinstance(par, int)) for par in [atol, rtol]]):
            self.atol = float(atol)
            self.rtol = float(rtol)
        else:
            raise TypeError(
                "Both atol and rtol must be an instance of either float or int.")

    def fitness(self):
        """This function will be redefined in the inheriting classes
        """
        pass

    def _plot_traj(self):
        """This function will be redefined in the inheriting classes
        """
        pass

    def _pretty(self):
        """This function will be redefined in the inheriting classes
        """
        pass

    def get_nobj(self):
        return 1

    def get_nec(self):
        return self.leg.nec

    def plot_traj(self, z, mark="k", atol=1e-12, rtol=1e-12, units=pk.AU, axes=None, quiver=False, length=1):
        """This function plots the 3 dimensional spacecraft trajectory, given a solution chromosome.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - atol (``float``, ``int``): absolute integration tolerance.
            - rtol (``float``, ``int``): relative integration tolerance.
            - units (``float``, ``int``): units by which to scale the trajectory dimensions.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - quiver (``bool``): when True the thrust is visualized with a quiver plot
            - length (``float``): Length of thrust arrow if quiver is True

        """

        # set problem
        self.fitness(z)

        # figure
        if axes is None:
            fig = plt.figure()
            axes = fig.gca(projection='3d')
        elif not isinstance(axes, Axes3D):
            raise TypeError(
                "Axis must be instance of matplotlib.axes._subplots.Axes3DSubplot.")

        # problem specifics
        self._plot_traj(z, axes, units)

        # Sun
        axes.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot_traj(axes, mark, atol, rtol, units,
                           quiver=quiver, length=length)

        return axes

    def plot_control(self, z, mark="k.-", atol=1e-12, rtol=1e-12, axes=None):
        """Plots the control profile of the trajectory, as a function of time.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - atol (``float``, ``int``): absolute integration tolerance.
            - rtol (``float``, ``int``): relative integration tolerance.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
        """

        # set problem
        self.fitness(z)

        # create figure
        if axes is None:
            fig = plt.figure()
            axes = fig.gca()

        # leg
        self.leg.plot('tof', 'u', mark=mark, atol=atol, rtol=rtol,
                      xlabel="Time [mjd2000]", ylabel="Throttle [ND]", axes=axes)

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

        # get states [t, x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm,
        # umag, ux, uy, uz, H]
        x = self.leg.get_states()
        # we make sure the throttles have the correct magnitude (get_states
        # returns a direction as defined by the primer)
        for t in x:
            t[-4] *= t[-5]
            t[-3] *= t[-5]
            t[-2] *= t[-5]
        return x

    def pretty(self, z):
        data = self.get_traj(z)
        self._pretty(z)

        print("\nSpacecraft Initial Position (m)  : [{!r}, {!r}, {!r}]".format(
            data[0, 1], data[0, 2], data[0, 3]))
        print("Spacecraft Initial Velocity (m/s)  : [{!r}, {!r}, {!r}]".format(
            data[0, 4], data[0, 5], data[0, 6]))
        print("Spacecraft Initial Mass  (kg)      : {!r}".format(data[0, 7]))

        print("Spacecraft Final Position (m)  : [{!r}, {!r}, {!r}]".format(
            data[-1, 1], data[-1, 2], data[-1, 3]))
        print("Spacecraft Final Velocity (m/s): [{!r}, {!r}, {!r}]".format(
            data[-1, 4], data[-1, 5], data[-1, 6]))
        print("Spacecraft Final Mass  (kg)    : {!r}".format(data[-1, 7]))
        print("Used propellant  (kg)          : {!r}".format(
            data[0, 7] - data[-1, 7]))


class indirect_pt2pt(_indirect_base):
    """
    Represents an indirect trajectory optimisation problem between two Cartesian states with heliocentric dynamics.
    The class can be used as UDP in pagmo.

    The decision chromosome is
    ::

        z = [T, l0]

    """

    def __init__(self,
                 x0=[-51051524893.335152, -142842795180.97464, 1139935.2553601924,
                     30488.847061907356, -10612.482697050367, -204.23284335657095, 1000],
                 xf=[24753885674.871033, 231247560000.17883, 4236305010.4256544, -
                     23171.900670190855, 4635.6817290400222, 666.44019588506023, 910.48383959441833],
                 thrust=0.3,
                 isp=3000,
                 mu=pk.MU_SUN,
                 tof=[276.15166075931495, 276.15166075931495],
                 freetime=False,
                 alpha=0,    # quadratic control
                 bound=False,
                 atol=1e-12,
                 rtol=1e-12):
        """
        Constructs an instance of the ``pykep.trajopt.indirect_pt2pt`` problem.

        Args:
            - x0 (``list``, ``tuple``, ``numpy.ndarray``): Departure state [m, m, m, m/s, m/s, m/s, kg].
            - xf (``list``, ``tuple``, ``numpy.ndarray``): Arrival state [m, m, m, m/s, m/s, m/s, kg].
            - tof (``list``): Transfer time bounds [days].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - mu (``float``): Gravitational parameter of primary body [m^3/s^2].
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parameter (0 -quadratic control, 1 - mass optimal)
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
        """

        # Cartesian states
        if not all([(isinstance(x, list) or isinstance(x, tuple) or isinstance(x, np.ndarray)) for x in [x0, xf]]):
            raise TypeError(
                "Both x0 and xf must be supplied as an instance of either list, tuple, or numpy.ndarray.")
        elif not all([len(x) == 7 for x in [x0, xf]]):
            raise TypeError(
                "Both x0 and xf must be supplied with 7 dimensions.")
        else:
            self.x0 = pk.sims_flanagan.sc_state()
            self.x0.set(x0)
            self.xf = pk.sims_flanagan.sc_state()
            self.xf.set(xf)

        self.tof = tof

        # initialise base
        _indirect_base.__init__(
            self, x0[-1], thrust, isp, mu, True, freetime, alpha, bound,
            atol, rtol
        )

    def fitness(self, z):

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # costates
        l0 = np.asarray(z[1:])

        # set leg
        self.leg.set(t0, self.x0, l0, tf, self.xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        return np.hstack(([1, ceq]))

    def get_bounds(self):
        lb = [self.tof[0]] + [-100] * 7
        ub = [self.tof[1]] + [100] * 7
        return (lb, ub)

    def _plot_traj(self, z, axes, units):

        # states
        x0 = self.leg.x0
        xf = self.leg.xf

        # times
        t0 = pk.epoch(self.leg.t0)
        tf = pk.epoch(self.leg.tf)

        # Computes the osculating Keplerian elements at start and arrival
        elem0 = list(pk.ic2par(x0[0:3], x0[3:6], self.leg.mu))
        elemf = list(pk.ic2par(xf[0:3], xf[3:6], self.leg.mu))

        # Converts the eccentric anomaly into eccentric anomaly
        elem0[5] = elem0[5] - elem0[1] * np.sin(elem0[5])
        elemf[5] = elemf[5] - elemf[1] * np.sin(elemf[5])

        # Creates two virtual keplerian planets with the said elements
        kep0 = pk.planet.keplerian(t0, elem0)
        kepf = pk.planet.keplerian(tf, elemf)

        # Plots the departure and arrival osculating orbits
        pk.orbit_plots.plot_planet(
            kep0, t0, units=units, color=(0.8, 0.8, 0.8), axes=axes)
        pk.orbit_plots.plot_planet(
            kepf, tf, units=units, color=(0.8, 0.8, 0.8), axes=axes)

    def _pretty(self, z):
        print("\nPoint to point transfer: ")
        print("\nFrom: " + str(self.x0))
        print("To: " + str(self.xf))
        print("Time of flight (days): {!r} ".format(z[0]))


class indirect_or2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between two orbits.

    Decision chromosome is
    ::

        z = [T, M0, Mf, l0]

    """

    def __init__(self,
                 elem0=[149598261129.93335, 0.016711230601231957,
                        2.640492490927786e-07, 3.141592653589793, 4.938194050401601, 0],
                 elemf=[227943822376.03537, 0.09339409892101332,
                        0.032283207367640024, 0.8649771996521327, 5.000312830124232, 0],
                 mass=1000,
                 thrust=0.3,
                 isp=2500,
                 atol=1e-12,
                 rtol=1e-12,
                 tof=[100, 700],
                 freetime=True,
                 alpha=0,
                 bound=False,
                 mu=pk.MU_SUN):
        """Initialises ``pykep.trajopt.indirect_or2or`` problem.

        Args:
            - elem0 (``list``, ``tuple``, ``numpy.ndarray``): Departure Keplerian elements (mutable eccentric anomaly).
            - elemf (``list``, ``tuple``, ``numpy.ndarray``): Arrival Keplerian elements (mutable eccentric anomaly).
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - tof (``list``): Transfer time bounds [days].
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parameter, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parameter of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, mu, True, freetime, alpha, bound,
            atol, rtol
        )

        # Keplerian elements
        self.elem0 = np.asarray(elem0)
        self.elemf = np.asarray(elemf)

        # Time of flight bounds
        self.tof = tof

    def fitness(self, z):

        # departure and arrival times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # departure and arrival eccentric anomolies
        M0 = z[1]
        Mf = z[2]

        # departure costates
        l0 = np.asarray(z[3:])

        # set Keplerian elements
        elem0 = np.hstack([self.elem0[:5], [M0]])
        elemf = np.hstack([self.elemf[:5], [Mf]])

        # compute Cartesian states
        r0, v0 = pk.par2ic(elem0, self.leg.mu)
        rf, vf = pk.par2ic(elemf, self.leg.mu)

        # departure and arrival states (xf[6] is unused)
        x0 = pk.sims_flanagan.sc_state(r0, v0, self.sc.mass)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        obj = self.leg.trajectory[-1, -1]

        # Transversality conditions
        # At start
        lambdas0 = np.array(self.leg.trajectory[0, 7:13])
        r0norm = np.sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2])
        tmp = - pk.MU_SUN / r0norm**3
        tangent = np.array([v0[0], v0[1], v0[2], tmp *
                            r0[0], tmp * r0[1], tmp * r0[2]])
        tangent_norm = np.linalg.norm(tangent)
        tangent = tangent / tangent_norm
        T0 = np.dot(lambdas0, tangent)

        # At end
        lambdasf = np.array(self.leg.trajectory[-1, 7:13])
        rfnorm = np.sqrt(rf[0] * rf[0] + rf[1] * rf[1] + rf[2] * rf[2])
        tmp = - pk.MU_SUN / rfnorm**3
        tangent = np.array([vf[0], vf[1], vf[2], tmp *
                            rf[0], tmp * rf[1], tmp * rf[2]])
        tangent_norm = np.linalg.norm(tangent)
        tangent = tangent / tangent_norm
        Tf = np.dot(lambdasf, tangent)

        return np.hstack(([1], ceq, [T0, Tf]))

    def get_nec(self):
        return self.leg.nec + 2

    def get_bounds(self):
        lb = [self.tof[0], -4 * np.pi, -4 * np.pi] + [-1e2] * 7
        ub = [self.tof[1], 4 * np.pi, 4 * np.pi] + [1e2] * 7
        return (lb, ub)

    def _plot_traj(self, z, axes, units):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(pykep.trajopt.indirect_or2or).plot_traj(pop.champion_x)
        """

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # Mean Anomalies
        M0 = z[1] - self.elem0[1] * np.sin(z[1])
        Mf = z[2] - self.elemf[1] * np.sin(z[2])

        elem0 = np.hstack([self.elem0[:5], [M0]])
        elemf = np.hstack([self.elemf[:5], [Mf]])

        # Keplerian points
        kep0 = pk.planet.keplerian(t0, elem0)
        kepf = pk.planet.keplerian(tf, elemf)

        # planets
        pk.orbit_plots.plot_planet(
            kep0, t0=t0, units=units, axes=axes, color=(0.8, 0.8, 0.8))
        pk.orbit_plots.plot_planet(
            kepf, t0=tf, units=units, axes=axes, color=(0.8, 0.8, 0.8))

    def _pretty(self, z):
        print("\nOrbit to orbit transfer: ")
        print("\nFrom: " + str(list(self.elem0)))
        print("To: " + str(list(self.elemf)))
        print("Time of flight (days): {!r} ".format(z[0]))
        print("Starting mean anomaly (rad): {!r} ".format(z[1]))
        print("Arrival mean anomaly (rad): {!r} ".format(z[2]))


class indirect_pt2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between a Cartesian state and an orbit.

    Decision chromosome is
    ::

        z = [T, Mf, l0]

    """

    def __init__(self, x0, elemf, mass, thrust, isp, atol, rtol, tof, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
        """Initialises ``pykep.trajopt.indirect_pt2or`` problem.

        Args:
            - x0 (``list``, ``tuple``, ``numpy.ndarray``): Departure state [m, m, m, m/s, m/s, m/s, kg].
            - elemf (``list``, ``tuple``, ``numpy.ndarray``): Arrival Keplerian elements SI units. (mean anomaly will be changed).
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - tof (``list``): Transfer time bounds [days].
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parameter, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parameter of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, mu, True, freetime, alpha, bound,
            atol, rtol
        )

        # departure state and arrival Keplerian elements
        self.x0 = np.asarray(x0, np.float64)
        self.elemf = np.asarray(elemf, np.float64)
        self.tof = tof

    def fitness(self, z):

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # final eccentric anomaly
        Mf = z[1]

        # initial costates
        l0 = np.asarray(z[2:])

        # set arrival Keplerian elements
        self.elemf[5] = Mf

        # departure state
        x0 = pk.sims_flanagan.sc_state(self.x0[0:3], self.x0[3:6], self.x0[6])

        # compute Cartesian arrival state
        rf, vf = pk.par2ic(self.elemf, self.leg.mu)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        # mf = self.leg.trajectory[-1, 6]

        # Transversality condition at the end
        lambdasf = np.array(self.leg.trajectory[-1, 7:13])
        rfnorm = np.sqrt(rf[0] * rf[0] + rf[1] * rf[1] + rf[2] * rf[2])
        tmp = - pk.MU_SUN / rfnorm**3
        tangent = np.array([vf[0], vf[1], vf[2], tmp *
                            rf[0], tmp * rf[1], tmp * rf[2]])
        tangent_norm = np.linalg.norm(tangent)
        tangent = tangent / tangent_norm
        Tf = np.dot(lambdasf, tangent)

        return np.hstack(([1], ceq, [Tf]))

    def get_nec(self):
        return self.leg.nec + 1

    def get_bounds(self):
        lb = [self.tof[0], -4 * np.pi] + [-1e2] * 7
        ub = [self.tof[1], 4 * np.pi] + [1e2] * 7
        return (lb, ub)

    def _plot_traj(self, z, axes, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(pykep.trajopt.indirect_pt2or).plot_traj(pop.champion_x)
        """

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[1])

        # Keplerian elements of the osculating orbit at start
        elem0 = list(pk.ic2par(self.x0[0:3], self.x0[3:6], self.leg.mu))
        # Eccentric to Mean Anomaly
        elem0[5] = elem0[5] - elem0[1] * np.sin(elem0[5])

        # Mean Anomaly at the target orbit
        Mf = z[1] - self.elemf[1] * np.sin(z[1])

        elemf = np.hstack([self.elemf[:5], [Mf]])

        # Keplerian elements
        kep0 = pk.planet.keplerian(t0, elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # plot departure and arrival
        pk.orbit_plots.plot_planet(
            kep0, t0, units=units, color=(0.8, 0.8, 0.8), axes=axes)
        pk.orbit_plots.plot_planet(
            kepf, tf, units=units, color=(0.8, 0.8, 0.8), axes=axes)

    def _pretty(self, z):
        """
        prob.pretty(x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Prints human readable information on the trajectory represented by the decision vector x
        """
        print("\nPoint to orbit transfer: ")
        print("\nFrom (cartesian): ", list(self.x0))
        print("To (osculating elements): ", list(self.elemf))
        print("Time of flight (days): {!r} ".format(z[0]))
        print("Arrival mean anomaly (rad): {!r} ".format(z[1]))


class indirect_pt2pl(_indirect_base):
    """
    Represents an indirect trajectory optimisation problem between a Cartesian state and a planet (rendezvous).
    Since the terminal conditions on the planet are not fixed, the transversality condition H=0 is deactivated
    and optimization of T happens via an explicit minimization of the objective (hybrid direct-indirect method)

    Decision chromosome is
    ::

        z = [T, l0]

    """

    def __init__(self,
                 x0=[-24482087316.947845, -150000284705.77328, -196089391.29376224,
                     31677.87649549203, -5859.747563624047, -351.75278222719828, 1000],
                 pf="mars",
                 mass=1000,
                 thrust=0.3,
                 isp=3000,
                 tof=[230, 280],
                 t0=1251.0286746844447,
                 mu=pk.MU_SUN,
                 alpha=0,
                 bound=False,
                 atol=1e-12,
                 rtol=1e-12
                 ):
        """Initialises ``pykep.trajopt.indirect_pt2or`` problem.

        Args:
            - x0 (``list``, ``tuple``, ``numpy.ndarray``): Departure state [m, m, m, m/s, m/s, m/s, kg].
            - pf (``str``): Arrival planet name. (will be used to construct a planet.jpl_lp object)            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - tof (``list``): Transfer time bounds [days].
            - t0 (``float``): launch epoch [MJD2000].
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parameter, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parameter of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, mu, True, False, alpha, bound,
            atol, rtol
        )

        # departure epoch
        self.t0 = pk.epoch(t0)
        # departure state
        self.x0 = np.asarray(x0, np.float64)
        # arrival planet
        self.pf = pk.planet.jpl_lp(pf)
        # bounds on the time of flight
        self.tof = tof
        # store the alfa value (immutable)
        self._alpha = alpha

    def fitness(self, z):

        # times
        t0 = self.t0
        tf = pk.epoch(t0.mjd2000 + z[0])

        # intial costates
        l0 = np.asarray(z[1:])

        # arrival conditions
        rf, vf = self.pf.eph(tf)

        # departure state
        x0 = pk.sims_flanagan.sc_state(self.x0[0:3], self.x0[3:6], self.x0[6])

        # arrival state (mass will be ignored)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        obj = self.leg.trajectory[-1, -1] * self.leg._dynamics.c2 * 1000

        return np.hstack(([obj], ceq))

    def get_bounds(self):
        lb = [self.tof[0]] + [-1e2] * 7
        ub = [self.tof[1]] + [1e2] * 7
        return (lb, ub)

    def _plot_traj(self, z, axes, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(pykep.trajopt.indirect_pt2or).plot_traj(pop.champion_x)
        """

        # states
        x0 = self.x0

        # times
        t0 = self.t0
        tf = pk.epoch(t0.mjd2000 + z[0])

        # Computes the osculating Keplerian elements at start
        elem0 = list(pk.ic2par(x0[0:3], x0[3:6], self.leg.mu))

        # Converts the eccentric anomaly into eccentric anomaly
        elem0[5] = elem0[5] - elem0[1] * np.sin(elem0[5])

        # Creates a virtual keplerian planet with the said elements
        kep0 = pk.planet.keplerian(t0, elem0)

        # Plots the departure and arrival osculating orbits
        pk.orbit_plots.plot_planet(
            kep0, t0, units=units, color=(0.8, 0.8, 0.8), axes=axes)
        pk.orbit_plots.plot_planet(
            self.pf, tf, units=units, color=(0.8, 0.8, 0.8), axes=axes)

    def _pretty(self, z):
        print("\nPlanet to orbit transfer, alpha is: ",  self._alpha)
        print("\nFrom (cartesian): " + str(list(self.x0)))
        print("Launch epoch: {!r} MJD2000, a.k.a. {!r}".format(
            self.t0.mjd2000, self.t0))
        print("\nTo (planet): " + self.pf.name)
        print("Time of flight (days): {!r} ".format(z[0]))
