import PyKEP as pk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class _indirect_base(object):

    def __init__(self, mass, thrust, isp, atol, rtol, freemass, freetime, alpha, bound, mu=pk.MU_SUN):

        # spacecraft
        self.sc = pk.sims_flanagan.spacecraft(mass, thrust, isp)

        # indirect leg
        self.leg = pk.pontryagin.leg(sc=self.sc, mu=mu, freetime=freetime, alpha=alpha, bound=bound)

        # integration parametres
        if all([(isinstance(par, float), isinstance(par, int)) for par in [atol, rtol]]):
            self.atol = float(atol)
            self.rtol = float(rtol)
        else:
            raise TypeError("Both atol and rtol must be an instance of either float or int.")

    def get_nobj(self):
        return 1

    def get_nec(self):
        return self.leg.nec

class indirect_pt2pt(_indirect_base):
    """Represents an indirect trajectory optimisation problem between two Cartesian states with heliocentric dynamics.

    This class represents the
    `NLP <https://esa.github.io/pagmo2/docs/python/tutorials/coding_udp_constrained.html>`_
    transcription of the indirect optimal control problem of optimising a
    spacecraft trajectory (alla moda di `Pontryagin's maximum principle <https://en.wikipedia.org/wiki/Pontryagin%27s_maximum_principle>`_)
    with respect to its fuel consumption between two Cartesian states, subject
    to heliocentric dynamics (gravitational influence of the sun) and propulsive
    control (firing of thrusters).

    The decision chromosome is
    ::

        z = [t0, T, l0]
    """

    def __init__(
        self, x0, xf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub,
        freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN
    ):
        """Initialises ``PyKEP.trajopt.indirect_pt2pt`` problem.

        Args:
            - x0 (``list``, ``tuple``, ``numpy.ndarray``): Departure state [m, m, m, m/s, m/s, m/s, kg].
            - xf (``list``, ``tuple``, ``numpy.ndarray``): Arrival state [m, m, m, m/s, m/s, m/s, kg].
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - t0lb (``float``, ``int``): Minimimum departure time [mjd2000].
            - t0ub (``float``, ``int``): Maximum departure time [mjd2000].
            - Tlb (``float``, ``int``): Minimimum time of flight [mjd2000].
            - Tub (``float``, ``int``): Maximum time of flight [mjd2000].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].
        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, atol, rtol, freemass, freetime, alpha, bound
        )

        # bounds
        self.t0lb = float(t0lb)
        self.t0ub = float(t0ub)
        self.Tlb = float(Tlb)
        self.Tub = float(Tub)

        # Cartesian states
        self.x0 = np.asarray(x0)
        self.xf = np.asarray(xf)

    def fitness(self, z):

        # times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # costates
        l0 = np.asarray(z[2:])

        # states
        x0 = pk.sims_flanagan.sc_state()
        x0.set(self.x0)
        xf = pk.sims_flanagan.sc_state()
        xf.set(self.xf)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf, ceq]))

    def get_bounds(self):
        lb = [self.t0lb, self.Tlb, *[-1e-2]*7]
        ub = [self.t0ub, self.Tub, *[1e-2]*7]
        return (lb, ub)

    def plot_traj(self, z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_pt2pt).plot_traj(pop.champion_x)
        """

        # set up figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # departure and arrival times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # departure costates
        l0 = z[2:]

        # set and plot leg
        self.leg.set(t0, self.x0, l0, tf, self.xf)
        self.leg.plot(axis, atol=atol, rtol=rtol, units=units)

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot(axis, mark, atol, rtol, units)

        # states
        x0 = self.leg.x0
        xf = self.leg.xf

        # Keplerian elements
        kep0 = pk.planet.keplerian(t0, pk.ic2par(x0[0:3], x0[3:6], self.mu))
        kepf = pk.planet.keplerian(tf, pk.ic2par(xf[0:3], xf[3:6], self.mu))

        # plot departure and arrival
        for kep, t in zip([kep0, kepf], [t0, tf]):
            pk.orbit_plots.plot_planet(kep, t, units=units, color=(0.8, 0.8, 0.8), ax=axis)

        # show plot
        plt.show()

class indirect_pl2pl(_indirect_base):
    """Represents an indirect trajectory optimisation problem between two planets.

    Decision chromosome is
    ::

        z = [t0, T, l0]

    """

    def __init__(self, p0, pf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
        """Initialises ``PyKEP.trajopt.indirect_pl2pl`` problem.

        Args:
            - p0 (``str``): Departure planet name.
            - pf (``str``): Arrival planet name.
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - t0lb (``float``, ``int``): Minimimum departure time [mjd2000].
            - t0ub (``float``, ``int``): Maximum departure time [mjd2000].
            - Tlb (``float``, ``int``): Minimimum time of flight [mjd2000].
            - Tub (``float``, ``int``): Maximum time of flight [mjd2000].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].
        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, atol, rtol, freemass, freetime, alpha, bound
        )

        # bounds
        self.t0lb = float(t0lb)
        self.t0ub = float(t0ub)
        self.Tlb = float(Tlb)
        self.Tub = float(Tub)

        # planets
        if all([isinstance(pl, str) for pl in [p0, pf]]):
            self.p0 = pk.planet.jpl_lp(p0)
            self.pf = pk.planet.jpl_lp(pf)
        else:
            raise TypeError("Planet names must be supplied as str.")

    def fitness(self, z):

        # times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # initial costates
        l0 = np.asarray(z[2:])

        # Cartesian states
        r0, v0 = self.p0.eph(t0)
        rf, vf = self.pf.eph(tf)

        # states
        x0 = pk.sims_flanagan.sc_state(r0, v0, self.sc.mass)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass/10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # propagate leg
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf, ceq]))

    def get_bounds(self):
        lb = [self.t0lb, self.Tlb, *[-1e2]*7]
        ub = [self.t0ub, self.Tub, *[1e2]*7]
        return (lb, ub)

    def plot_traj(self, z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_pl2pl).plot_traj(pop.champion_x)
        """

        # set up figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # departure and arrival times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # departure costates
        l0 = z[2:]

        # set and plot leg
        self.leg.set(t0, self.x0, l0, tf, self.xf)
        self.leg.plot(axis, atol=atol, rtol=rtol, units=units)

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot(axis, mark, atol, rtol, units)

        # planets
        pk.orbit_plots.plot_planet(self.p0, t0=t0, units=units, ax=axis, color=(0.8, 0.8, 0.8))
        pk.orbit_plots.plot_planet(self.pf, t0=tf, units=units, ax=axis, color=(0.8, 0.8, 0.8))

        # show plot
        plt.show()

class indirect_or2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between two orbits.

    Decision chromosome is
    ::

        z = [T, M0, Mf, l0]

    """

    def __init__(self, elem0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
        """Initialises ``PyKEP.trajopt.indirect_or2or`` problem.

        Args:
            - elem0 (``list``, ``tuple``, ``numpy.ndarray``): Departure Keplerian elements (mutable mean anomoly).
            - elemf (``list``, ``tuple``, ``numpy.ndarray``): Arrival Keplerian elements (mutable mean anomoly).
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - Tlb (``float``, ``int``): Minimimum time of flight [mjd2000].
            - Tub (``float``, ``int``): Maximum time of flight [mjd2000].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, atol, rtol, freemass, freetime, alpha, bound
        )

        # bounds
        self.Tlb = float(Tlb)
        self.Tub = float(Tub)

        # Keplerian elements
        self.elem0 = np.asarray(elem0)
        self.elemf = np.asarray(elemf)

    def fitness(self, z):

        # departure and arrival times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # departure and arrival mean anomolies
        M0 = z[1]
        Mf = z[2]

        # departure costates
        l0 = np.asarray(z[3:])

        # set Keplerian elements
        self.elem0[5] = M0
        self.elemf[5] = Mf

        # compute Cartesian states
        r0, v0 = pk.par2ic(self.elem0, self.leg.mu)
        rf, vf = pk.par2ic(self.elemf, self.leg.mu)

        # departure and arrival states (xf[6] is unused)
        x0 = pk.sims_flanagan.sc_state(r0, v0, self.sc.mass)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass/10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf], ceq))

    def get_bounds(self):
        pi = 3.14159265359
        lb = [self.Tlb, -4*pi, -4*pi, *[-1e2]*7]
        ub = [self.Tub, 4*pi, 4*pi, *[1e2]*7]
        return (lb, ub)

    def plot_traj(self, z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_or2or).plot_traj(pop.champion_x)
        """

        # set up figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # set problem
        self.fitness(z)

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot(axis, mark, atol, rtol, units)

        # Keplerian points
        kep0 = pk.planet.keplerian(t0, self.elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # planets
        pk.orbit_plots.plot_planet(kep0, t0=t0, units=units, ax=axis, color=(0.8, 0.8, 0.8))
        pk.orbit_plots.plot_planet(kepf, t0=tf, units=units, ax=axis, color=(0.8, 0.8, 0.8))

        # show plot
        plt.show()

class indirect_pt2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between a Cartesian state and an orbit.

    Decision chromosome is
    ::

        z = [T, Mf, l0]

    """

    def __init__(self, x0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
        """Initialises ``PyKEP.trajopt.indirect_pt2or`` problem.

        Args:
            - x0 (``list``, ``tuple``, ``numpy.ndarray``): Departure state [m, m, m, m/s, m/s, m/s, kg].
            - elemf (``list``, ``tuple``, ``numpy.ndarray``): Arrival Keplerian elements (mutable mean anomoly).
            - mass (``float``, ``int``): Spacecraft wet mass [kg].
            - thrust (``float``, ``int``): Spacecraft maximum thrust [N].
            - isp (``float``, ``int``): Spacecraft specific impulse [s].
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - Tlb (``float``, ``int``): Minimimum time of flight [mjd2000].
            - Tub (``float``, ``int``): Maximum time of flight [mjd2000].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, atol, rtol, freemass, freetime, alpha, bound, mu
        )

        # bounds
        self.Tlb = float(Tlb)
        self.Tub = float(Tub)

        # departure state and arrival Keplerian elements
        self.x0 = np.asarray(x0)
        self.elemf = np.asarray(elemf)

    def fitness(self, z):

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # final mean anomoly
        Mf = z[1]

        # intial costates
        l0 = np.asarray(z[2:])

        # set arrival Keplerian elements
        self.elemf[5] = Mf

        # departure state
        x0 = pk.sims_flanagan.sc_state(self.x0[0:3], self.x0[3:6], self.x0[6])

        # compute Cartesian arrival state
        rf, vf = pk.par2ic(self.elemf, self.leg.mu)
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass/10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf], ceq))

    def get_bounds(self):
        pi = 3.14159265359
        lb = [self.Tlb, -4*pi, *[-1e2]*7]
        ub = [self.Tub, 4*pi, *[1e2]*7]
        return (lb, ub)

    def plot_traj(self, z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_pt2or).plot_traj(pop.champion_x)
        """

        # set up figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # set problem
        self.fitness(z)

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot(axis, mark, atol, rtol, units)

        # Keplerian points
        kep0 = pk.planet.keplerian(t0, pk.ic2par(self.x0[0:3], self.x0[3:6], mu=self.mu))
        kepf = pk.planet.keplerian(tf, self.elemf)

        # plot planets
        pk.orbit_plots.plot_planet(kep0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(kepf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)

        plt.show()
