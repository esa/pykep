import PyKEP as pk
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class _indirect_base(object):
    """Base class for indirect trajectory optimisation problems.

    All inheriting classes will adopt ``plot_traj`` and ``plot_control``.
    """

    def __init__(
        self, mass, thrust, isp,
        mu, freemass, freetime, alpha, bound,
        atol, rtol,
        t0lb=None, t0ub=None, Tlb=None, Tub=None,
        M0lb=None, M0ub=None, Mflb=None, Mfub=None
    ):

        # spacecraft
        self.sc = pk.sims_flanagan.spacecraft(mass, thrust, isp)

        # indirect leg
        self.leg = pk.pontryagin.leg(
            sc=self.sc, mu=mu, freemass=freemass, freetime=freetime, alpha=alpha, bound=bound
        )

        # integration parametres
        if all([(isinstance(par, float) or isinstance(par, int)) for par in [atol, rtol]]):
            self.atol = float(atol)
            self.rtol = float(rtol)
        else:
            raise TypeError(
                "Both atol and rtol must be an instance of either float or int.")

        # bounds
        lbs = ["t0lb", "Tlb", "M0lb", "Mflb"]
        ubs = ["t0ub", "Tub", "M0ub", "Mfub"]

        for lb, ub in zip(lbs, ubs):

            # values
            lbval = eval(lb)
            ubval = eval(ub)
            bvals = [lbval, ubval]

            # if the pair isn't supplied, ignore
            if any([b is None for b in bvals]):
                continue

            # check types
            elif not all([(isinstance(b, float) or isinstance(b, int)) for b in bvals]):
                raise TypeError("Lower and uper bounds must be supplied as an instance of float or int.")

            # check validity
            elif not lbval <= ubval:
                raise TypeError("Lower bound must be less than upper bound.")

            # set attribute
            else:
                [setattr(self, b, bval) for b, bval in zip([lb, ub], bvals)]

    def get_nobj(self):
        return 1

    def get_nec(self):
        return self.leg.nec

    def plot_traj(self, z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU):
        """This function plots the 3 dimensional spacecraft trajectory, given a solution chromosome.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - atol (``float``, ``int``): absolute integration tolerance.
            - rtol (``float``, ``int``): relative integration tolerance.
            - units (``float``, ``int``): units by which to scale the trajectory dimensions.

        """

        # set problem
        self.fitness(z)

        # figure
        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # sun
        axis.scatter([0], [0], [0], color='y')

        # leg
        self.leg.plot_traj(axis, mark, atol, rtol, units)

        # problem specifics
        self._plot_traj(z, axis, units)

        # show
        plt.show()

    def plot_control(self, z, mark="k.-", atol=1e-12, rtol=1e-12):
        """Plots the control profile of the trajectory, as a function of time.

        Args:
            - z (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - mark (``string``): matplotlib marker.
            - atol (``float``, ``int``): absolute integration tolerance.
            - rtol (``float``, ``int``): relative integration tolerance.

        """

        # set problem
        self.fitness(z)

        # leg
        self.leg.plot('t', 'u', mark=mark, atol=atol, rtol=rtol, xlabel="Time [mjd2000]", ylabel="Throttle [ND]")

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

    def __init__(self, x0, xf, mass, thrust, isp, mu, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, atol=1e-10, rtol=1e-10):
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
            self, mass, thrust, isp, mu, freemass, freetime, alpha, bound,
            atol, rtol, t0lb=t0lb, t0ub=t0ub, Tlb=Tlb, Tub=Tub
        )

        # Cartesian states
        if not all([(isinstance(x, list) or isinstance(x, tuple) or isinstance(x, np.ndarray)) for x in [x0, xf]]):
            raise TypeError("Both x0 and xf must be supplied as an instance of either list, tuple, or numpy.ndarray.")
        elif not all([len(x) == 7 for x in [x0, xf]]):
            raise TypeError("Both x0 and xf must be supplied with 7 dimensions.")
        else:
            self.x0 = pk.sims_flanagan.sc_state()
            self.x0.set(x0)
            self.xf = pk.sims_flanagan.sc_state()
            self.xf.set(xf)

    def fitness(self, z):

        # times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # costates
        l0 = np.asarray(z[2:])

        # set leg
        self.leg.set(t0, self.x0, l0, tf, self.xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf, ceq]))

    def get_bounds(self):
        lb = [self.t0lb, self.Tlb, *[-1e2] * 7]
        ub = [self.t0ub, self.Tub, *[1e2] * 7]
        return (lb, ub)

    def _plot_traj(self, z, axis, units):

        # states
        x0 = self.leg.x0
        xf = self.leg.xf

        # times
        t0 = pk.epoch(self.leg.t0)
        tf = pk.epoch(self.leg.tf)

        # Keplerian elements
        elem0 = pk.ic2par(x0[0:3], x0[3:6], self.leg.mu)
        elemf = pk.ic2par(xf[0:3], xf[3:6], self.leg.mu)

        # Keplerian elements
        kep0 = pk.planet.keplerian(t0, elem0)
        kepf = pk.planet.keplerian(tf, elemf)

        # plot departure and arrival
        pk.orbit_plots.plot_planet(kep0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(kepf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)

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
            self, mass, thrust, isp, mu, freemass, freetime, alpha, bound,
            atol, rtol, t0lb=t0lb, t0ub=t0ub, Tlb=Tlb, Tub=Tub
        )

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
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # propagate leg
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf, ceq]))

    def get_bounds(self):
        lb = [self.t0lb, self.Tlb, *[-1e2] * 7]
        ub = [self.t0ub, self.Tub, *[1e2] * 7]
        return (lb, ub)

    def _plot_traj(self, z, axis, units):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_pl2pl).plot_traj(pop.champion_x)
        """

        # departure and arrival times
        t0 = pk.epoch(z[0])
        tf = pk.epoch(z[0] + z[1])

        # planets
        pk.orbit_plots.plot_planet(
            self.p0, t0=t0, units=units, ax=axis, color=(0.8, 0.8, 0.8))
        pk.orbit_plots.plot_planet(
            self.pf, t0=tf, units=units, ax=axis, color=(0.8, 0.8, 0.8))

class indirect_or2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between two orbits.

    Decision chromosome is
    ::

        z = [T, M0, Mf, l0]

    """

    def __init__(self, elem0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
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
            - M0lb (``float``, ``int``): Minimum departure mean anomoly [rad].
            - M0ub (``float``, ``int``): Maximum departure mean anomoly [rad].
            - Mflb (``float``, ``int``): Minimum arrival mean anomoly [rad].
            - M0fb (``float``, ``int``): Maximum arrival mean anomoly [rad].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.
            - mu (``float``): Gravitational parametre of primary body [m^3/s^2].

        """

        # initialise base
        _indirect_base.__init__(
            self, mass, thrust, isp, mu, freemass, freetime, alpha, bound,
            atol, rtol, Tlb=Tlb, Tub=Tub, M0lb=M0lb, M0ub=M0ub, Mflb=Mflb, Mfub=Mfub
        )

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
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf], ceq))

    def get_bounds(self):
        lb = [self.Tlb, self.M0lb, self.Mflb, *[-1e2] * 7]
        ub = [self.Tub, self.M0ub, self.Mfub, *[1e2] * 7]
        return (lb, ub)

    def _plot_traj(self, z, axis, units):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_or2or).plot_traj(pop.champion_x)
        """

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[0])

        # Keplerian points
        kep0 = pk.planet.keplerian(t0, self.elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # planets
        pk.orbit_plots.plot_planet(
            kep0, t0=t0, units=units, ax=axis, color=(0.8, 0.8, 0.8))
        pk.orbit_plots.plot_planet(
            kepf, t0=tf, units=units, ax=axis, color=(0.8, 0.8, 0.8))

class indirect_pt2or(_indirect_base):
    """Represents an indirect trajectory optimisation problem between a Cartesian state and an orbit.

    Decision chromosome is
    ::

        z = [T, Mf, l0]

    """

    def __init__(self, x0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN):
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
            self, mass, thrust, isp, mu, freemass, freetime, alpha, bound,
            atol, rtol, Tlb=Tlb, Tub=Tub, Mflb=Mflb, Mfub=Mfub
        )


        # departure state and arrival Keplerian elements
        self.x0 = np.asarray(x0, np.float64)
        self.elemf = np.asarray(elemf, np.float64)

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
        xf = pk.sims_flanagan.sc_state(rf, vf, self.sc.mass / 10)

        # set leg
        self.leg.set(t0, x0, l0, tf, xf)

        # equality constraints
        ceq = self.leg.mismatch_constraints(atol=self.atol, rtol=self.rtol)

        # final mass
        mf = self.leg.trajectory[-1, 6]

        return np.hstack(([-mf], ceq))

    def get_bounds(self):
        pi = 3.14159265359
        lb = [self.Tlb, self.Mflb, *[-1e2] * 7]
        ub = [self.Tub, self.Mfub, *[1e2] * 7]
        return (lb, ub)

    def _plot_traj(self, z, axis, units=pk.AU):
        """Plots spacecraft trajectory.

        Args:
            - z (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Examples:
            >>> prob.extract(PyKEP.trajopt.indirect_pt2or).plot_traj(pop.champion_x)
        """

        # times
        t0 = pk.epoch(0)
        tf = pk.epoch(z[1])

        # Keplerian elements
        elem0 = pk.ic2par(self.x0[0:3], self.x0[3:6], self.leg.mu)

        # Keplerian elements
        kep0 = pk.planet.keplerian(t0, elem0)
        kepf = pk.planet.keplerian(tf, self.elemf)

        # plot departure and arrival
        pk.orbit_plots.plot_planet(kep0, t0, units=units, color=(0.8, 0.8, 0.8), ax=axis)
        pk.orbit_plots.plot_planet(kepf, tf, units=units, color=(0.8, 0.8, 0.8), ax=axis)
