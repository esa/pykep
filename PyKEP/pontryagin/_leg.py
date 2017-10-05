from _dynamics import _dynamics
import PyKEP as pk
from scipy.integrate import ode
import numpy as np


class leg(object):
    """Indirect optimal control transcription trajectory leg.

    This class represents an indirect optimal control transcription
    (alla moda di Pontryagin's maximum principle) of a generic two-point
    boundary trajectory.

    Attributes:
        t0 (`float`): Departure time [mjd2000].
        x0 (`numpy.ndarray`): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
        l0 (`numpy.ndarray`): Costate variables [ND, ND, ND, ND, ND, ND, ND].
        tf (`float`): Arrival time [mjd2000].
        xf (`numpy.ndarray`): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].
        sc (`pk.sims_flanagan.spacecraft`): Generic spacecraft with propulsive properties.
        mu (`float`): Gravitational parametre of primary body [m^3/s^2]
        freemass (`bool`): Activates final mass transversality condition.
            Allows final mass to vary.
        freetime (`bool`): Activates final time transversality condition.
            Allows final time to vary.
        alpha (`float`): Homotopy parametre, governing the degree
            to which the theoretical control law is intended to reduce
            propellant expenditure or energy.
            Setting the parametre to 1 enforces a mass-optimal control law,
            with a characteristic bang-bang control profile (either full throttle or off).
            Setting the parametre to 0 enforces a pure quadratic control law.
        bound (`bool`): Activates bounded control, in which the control
            throttle is bounded between 0 and 1, otherwise the control
            throttle is allowed to unbounded.
        nec (`int`): Number of equality constraints as determined by
            `freemass` and `freetime`.

    """

    def __init__(
        self, t0=None, x0=None, l0=None, tf=None, xf=None,
        sc=pk.sims_flanagan.spacecraft(1000, 0.3, 2500), mu=pk.MU_SUN,
        freemass=True, freetime=True, alpha=1, bound=True
    ):
        """Initialises an indirect optimal control transcription trajectory leg.

        Args:
            t0 (`PyKEP.epoch`, `None`): Departure time [mjd2000].
            x0 (`PyKEP.sims_flanagan.sc_state`, `None`): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
            l0 (`numpy.ndarray`, `list`, `tuple`, `None`): Costate variables [ND, ND, ND, ND, ND, ND, ND].
            tf (`PyKEP.epoch`, `None`): Arrival time [mjd2000].
            xf (`PyKEP.epoch`, `None`): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].
            sc (`PyKEP.sims_flanagan.spacecraft`): Generic spacecraft with propulsive properties.
            mu (`float`, `int`): Gravitational parametre of primary body [m^3/s^2].
            freemass (`bool`): Activates final mass transversality condition.
                Allows final mass to vary.
            freetime (`bool`): Activates final time transversality condition.
                Allows final time to vary.
            alpha (`float`, `int`): Homotopy parametre, governing the degree
                to which the theoretical control law is intended to reduce
                propellant expenditure or energy.
                Setting the parametre to 1 enforces a mass-optimal control law,
                with a characteristic bang-bang control profile (either full throttle or off).
                Setting the parametre to 0 enforces a pure quadratic control law.
            bound (`bool`): Activates bounded control, in which the control
                throttle is bounded between 0 and 1, otherwise the control
                throttle is allowed to unbounded.

        Raises:
            TypeError: If `sc` is not supplied as an instance of `PyKEP.sims_flanagan.spacecraft`.
            TypeError: If `mu` is neither supplied as an instance of `float` or `int`.
            ValueError: If `mu` is not supplied as a positive number.
            TypeError: If either `freemass`, `freetime`, or `bound` is not supplied as an instance of `bool`.
            AttributeError: If equality constraint dimension cannot be determined form `freemass` and `freetime`.
            TypeError: If `alpha` is neither supplied as an instance of `float` or `int`.
            ValueError: If `alpha` is not in between 0 and 1.
            ValueError: If `alpha == 1` and `bound == False`.
                Control cannot be unbounded with mass-optimal control.
            TypeError: If either `t0` or `tf` is not supplied as an instance of `PyKEP.epoch`.
            ValueError: If departure time, `t0.mjd2000` does not occur before arrival time, `tf.mjd2000`.
            TypeError: If either `x0` or xf` is not supplied as an instance of `PyKEP.sims_flanagan.sc_state`.
            TypeError: If costate, `l0`, is neither supplied as an instance of `numpy.ndarray`, `list`, or `tuple`.
            TypeError: If costate, `l0`, is does not have 7 elements.
                Each element of `l0` corresponds to the respective elements of `x0`.

        Examples:

            >>> l = PyKEP.pontryagin.leg()
            >>> l.mismatch_constraints()

            >>> l = PyKEP.pontryagin.leg(t0, x0, l0, tf, xf)

        Note:
            - Typically spacecraft trajectories are optimised to reduce final mass, which theoretically results in a bang-bang control profile. However, a bang-bang control profile (either 0 or 1) is problematic to an optimiser due to its discontinous nature. The trajectory optimisation problem becomes easier when one first optimises with unbounded quadratic control (`alpha == 0 and bound == False`), stores the solution, then uses the solution to optimise with bounded quadratic control (`alpha == 0 and bound == True`). Then using the solution from bounded quadratic control, one can incrementally optimise for increasing homotopy parametres (e.g. `alphas = numpy.linspace(0, 1, 200)`) until a mass-optimal control solution converges.
            - Boundary conditions `t0`, `x0`, `l0`, `tf`, and `xf` are optional in the constructor. If some, but not all, boundary conditions are supplied, they are not set and disregarded. If the boundary conditions are not set upon construction (`__init__`), they must be set with `set(t0, x0, l0, tf, xf)`.
            - `propagate` will not work unless the boundary conditions have been set either through `__init__` or `set`.
        """

        # check spacecraft
        if isinstance(sc, pk.sims_flanagan.spacecraft):
            self.spacecraft = sc
        else:
            raise TypeError(
                "Spacecraft, sc, must be of pontryagin.spacecraft type.")

        # check gravitational parametre
        if not (isinstance(mu, float) or not isinstance(mu, int)):
            raise TypeError(
                "Gravitational parametre, mu, must be supplied as float or int.")
        elif not mu > 0:
            raise ValueError("Gravitational parametre, mu, must be positive.")
        else:
            self.mu = float(mu)

        # check freetime, freemass, and bound
        if not all([isinstance(param, bool) for param in [freemass, freetime, bound]]):
            raise TypeError(
                "Freemass, freetime, bound parametres must supplied as booleans.")
        else:
            self.freemass = bool(freemass)
            self.freetime = bool(freetime)
            self.bound    = bool(bound)
            print(self.freetime)

        # equality constraint dimensionality
        if self.freemass and self.freetime:
            self.nec = 8
        elif self.freemass and not self.freetime:
            self.nec = 7
        elif not self.freemass and self.freetime:
            self.nec = 8
        elif not self.freemass and not self.freetime:
            self.nec = 7
        else:
            raise AttributeError("Could not determine equality constraint dimensionality.")

        # check alpha
        if not (isinstance(alpha, float) or isinstance(alpha, int)):
            raise TypeError(
                "Homotopy parametre, alpha, must be supplied as float or int.")
        elif not (alpha >= 0 and alpha <= 1):
            raise ValueError(
                "Homotopy parametre, alpha, must be between 0 and 1.")
        elif (alpha == 1 and self.bound == False):
            raise ValueError(
                "If homotopy parametre, alpha, is 1, control must be bounded.")
        else:
            self.alpha = float(alpha)

        # if any of the necessary boundaries are not supplied
        if any(elem is None for elem in [t0, x0, l0, tf, xf]):
            pass

        # check validity if they're all supplied
        else:

            # check departure and arrival times
            if not all([isinstance(t, pk.epoch) for t in [t0, tf]]):
                raise TypeError(
                    "Departure and arrival times, t0 & tf, must be supplied as PyKEP.epoch.")
            elif not t0.mjd2000 < tf.mjd2000:
                raise ValueError(
                    "Departure time must occur before arrival time.")
            else:
                self.t0 = float(t0.mjd2000)
                self.tf = float(tf.mjd2000)

            # check departure and arrival states
            states = [x0, xf]
            if not all([isinstance(state, pk.sims_flanagan.sc_state) for state in states]):
                raise TypeError(
                    "x0 and xf must be instances of PyKEP.sims_flanagan.sc_state.")
            else:
                self.x0 = np.asarray(x0.get(), np.float64)
                self.xf = np.asarray(xf.get(), np.float64)

            # check departure costate
            if not (isinstance(l0, list) or isinstance(l0, np.ndarray) or isinstance(l0, tuple)):
                raise TypeError(
                    "Departure costate, l0, must be supplied as either a list, numpy.ndarray, or tuple.")
            elif len(l0) != 7:
                raise TypeError("Costate vector, l0, must be 7-dimensional.")
            else:
                self.l0 = np.asarray(l0, np.float64)

        # dynamics
        self._dynamics = _dynamics(sc=self.spacecraft, mu=mu)

        # integrator
        self._integrator = ode(
            lambda t, fs: self._dynamics._eom_fullstate(fs),
            lambda t, fs: self._dynamics._eom_fullstate_jac(fs)
        )

    def _recorder(self, t, fs):

        # append time
        self._times = np.append(self._times, t)

        # append fullstate
        self._trajectory = np.vstack((self._trajectory, fs))

    def set(self, t0, x0, l0, tf, xf):
        """Sets the departure and arrival boundary conditions of the trajectory.

        This instance method sets the departure time `t0`, departure state `x0`,
        departure costate `l0`, arrival time `tf`, arrival state `xf`.

        Args:
            t0 (`PyKEP.epoch`): Departure time [mjd2000].
            x0 (`PyKEP.sims_flanagan.sc_state`): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
            l0 (`numpy.ndarray`, `list`, `tuple`): Costate variables [ND, ND, ND, ND, ND, ND, ND].
            tf (`PyKEP.epoch`): Arrival time [mjd2000].
            xf (`PyKEP.epoch`): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].

        Raises:
            TypeError: If either `t0` or `tf` is not supplied as an instance of `PyKEP.epoch`.
            TypeError: If either `x0` or `xf` is not supllied as an instance of `PyKEP.sims_flanagan.sc_state`.
            TypeError: If `l0` is neither supllied as a `numpy.ndarray`, `list`, or `tuple`.

        Note:
            If ``Args`` were not supplied in `__init__`, `set` must be called before calling `mismatch_constraints`.

        Examples:

            >>> l = leg() # t0, x0, l0, tf, xf not set yet
            >>> l.set(t0, x0, l0, tf, xf) # now they're set

        """

        # check inputs
        if not all([isinstance(t, pk.epoch) for t in [t0, tf]]):
            raise TypeError("t0 and tf must be instances of PyKEP.epoch.")
        elif not all([isinstance(state, pk.sims_flanagan.sc_state) for state in [x0, xf]]):
            raise TypeError(
                "x0 and xf must be instances of PyKEP.sims_flanagan.sc_state.")
        elif not (isinstance(l0, tuple) or isinstance(l0, list) or isinstance(l0, np.ndarray)):
            raise TypeError(
                "l0 must be either a tuple, list, or numpy.ndarray.")
        elif len(l0) != 7:
            raise TypeError("Costate vector must be 7-dimensional.")
        else:

            # departure
            self.t0 = float(t0.mjd2000)
            self.x0 = np.asarray(x0.get(), np.float64)
            self.l0 = np.asarray(l0, np.float64)

            # arrival
            self.tf = float(tf.mjd2000)
            self.xf = np.asarray(xf.get(), np.float64)

    def _propagate(self, atol, rtol):

        # nondimensionalise departure state
        x0       = self.x0
        x0[0:3] /= self._dynamics.L
        x0[3:6] /= self._dynamics.V
        x0[6]   /= self._dynamics.M

        # nondimensional fullstate
        fs0 = np.hstack((x0, self.l0))

        # convert mjd2000 to mjs2000
        t0 = self.t0*24*60*60
        tf = self.tf*24*60*60

        # nondimensionalise times
        t0 /= self._dynamics.T
        tf /= self._dynamics.T

        # clear trajectory history
        self._times = np.empty((1, 0), dtype=np.float64)
        self._trajectory = np.empty((0, 14), dtype=np.float64)

        # set integration method
        self._integrator.set_integrator("dopri5", atol=atol, rtol=rtol)

        # set recorder
        self._integrator.set_solout(self._recorder)

        # set initial conditions
        self._integrator.set_initial_value(fs0, t0)

        # numerically integrate
        self._integrator.integrate(tf)

    def mismatch_constraints(self, atol=1e-5, rtol=1e-5):
        """Returns the nondimensional mismatch equality constraints of the arrival boundary conditions.

        This method propagates the nondimensional dynamics of the spacecraft
        from the departure time `t0` to the arrival time `tf`, then evaluates
        the nondimensional mismatch between the desired arrival state `xf`,
        arrival mass costate `lmf == 0` (`if freemass == True`), and arrival
        Hamiltonian `H == 0` (`if freetime == True`) and the supplied desired
        arrival boundary conditions.

        Args:
            atol (`float`): Absolute integration solution tolerance.
            rtol (`float`): Relative integration solution tolerance.

        Returns:
            `numpy.ndarray`: Equality constraints vector composed of the
            arrival mismatch in position & velocity, arrival mass costate if
            `freemass == True`, arrival mismatch in mass if `freemass == False`,
            and arrival Hamiltonian if `freetime == True`.::

            ceq = [drf, dvf, dmf]    # freemass == False; freetime == False
            ceq = [drf, dvf, lmf]    # freemass == True; freetime == False
            ceq = [drf, dvf, lmf, H] # freemass == True; freetime == True

        Raises:
            TypeError: If either `atol` or `rtol` is supplied as neither an instance of `float` or `int`.
            AttributeError: If boundary conditions `t0`, `x0`, `l0`, `tf`, and `x0` have not been set through either `__init__` or `set`.

        Note:
            - This method uses an explicit Runge-Kutta method of order (4)5
            due to Dormand & Prince with adaptive stepsize control.
            Smaller values of `atol` and `rtol` will increase the accuracy of
            a converged trajectory optimisation solution. However, smaller
            values result in slower integration executions and thus slower
            optimiser iterations. It is preferable to first solve a trajectory
            optimisation problem with higher values of `atol` and `rtol` to
            quickly converge to an approximate solution, subsequently resolving
            the problem with smaller tolerance values using the previous solution.
            This method is akin to mesh refinement in direct trajectory optimisation.
            - State arrival boundary conditions in position and velocity are
            always returned. If the final mass transversality condition is
            activated (`freemass == True`), than the final mass costate is
            returned as part of the equality constraint vector. Otherwise
            (`freemass == False`), the final mass mismatch is returned as part
            of the equality constraint vector. If the final time transversality
            condition is activated `freetime == True`, then than the final
            Hamiltonian will be returned as part of the equality constraint
            vector. Otherwise (`freetime == False`), the Hamiltonian will not
            be supplied as part of the equality constraint vector.

        Examples:

            >>> l = PyKEP.pontryagin.leg(freemass=True, freetime=True)
            >>> l.set(t0, x0, l0, tf, xf)
            >>> l.mismatch_constraints(atol=1e-10, rtol=1e-10)
            [ -1.15402379e+01  -2.23345886e+00  -9.30022917e-01  -2.53440778e+00
              -3.44246359e+00  -3.96669697e-01  -2.82967399e+03  -8.58526037e-01]
            >>> l.nec
            8

            >>> l = PyKEP.pontryagin.leg(freemass=False, freetime=True)
            >>> l.set(t0, x0, l0, tf, xf)
            >>> l.mismatch_constraints(atol=1e-10, rtol=1e-10)
            [-0.82435194  0.4375439   0.04615264  0.69192818 -0.18327442 -0.00930848
              0.46335841  1.66752446]
            >>> l.nec
            8

            >>> l = PyKEP.pontryagin.leg(freemass=True, freetime=False)
            >>> l.set(t0, x0, l0, tf, xf)
            >>> l.mismatch_constraints(atol=1e-10, rtol=1e-10)
            [  4.95782514e+00  -7.94403974e+00  -1.10158930e-02   7.27278513e+00
              -1.70998019e+00   9.13925064e-02  -1.45194673e+06]
            >>> l.nec
            7

            >>> l = PyKEP.pontryagin.leg(freemass=False, freetime=False)
            >>> l.set(t0, x0, l0, tf, xf)
            >>> l.mismatch_constraints(atol=1e-10, rtol=1e-10)
            [-0.90868565  0.23238016  0.04596579  0.61543688 -0.50023407 -0.0058185
              0.62170522]
            >>> l.nec
            7

            >>> l = PyKEP.pontryagin.leg()
            >>> l.mismatch_constraints(atol=1e-10, rtol=1e-10)
            AttributeError: Cannot propagate dynamics, as boundary conditions t0, x0, l0, tf, and xf have not been set. Use set(t0, x0, l0, tf, xf) to set boundary conditions.
        """

        if not all([(isinstance(tol, float) or isinstance(tol, int)) for tol in [atol, rtol]]):
            raise TypeError("Both atol and rtol must be supplied as instances of either float or int.")
        if any([hasattr(self, atr) == False for atr in ["t0", "x0", "l0", "tf", "xf"]]):
            raise AttributeError("Cannot propagate dynamics, as boundary conditions t0, x0, l0, tf, and xf have not been set. Use set(t0, x0, l0, tf, xf) to set boundary conditions.")
        else:
            atol = float(atol)
            rtol = float(atol)

        # propagate trajectory
        self._propagate(atol=atol, rtol=rtol)

        # nondimensional position and velocity arrival mismatch
        drf = self._trajectory[-1, 0:3] - (self.xf[0:3] / self._dynamics.L)
        dvf = self._trajectory[-1, 3:6] - (self.xf[3:6] / self._dynamics.V)

        # free arrival mass
        if self.freemass:
            lmf = self._trajectory[-1, 13]
        # fixed arrival mass
        else:
            dmf = self._trajectory[-1, 6] - self.xf[6] / self._dynamics.M

        # free arrival time
        if self.freetime:
            H = self._dynamics._hamiltonian(self._trajectory[-1])
        # fixed arrival time
        else:
            pass

        # create equality constraints
        if (self.freemass and self.freetime):
            ceq = np.hstack((drf, dvf, [lmf], [H]))
        elif (self.freemass and not self.freetime):
            ceq = np.hstack((drf, dvf, [lmf]))
        elif (self.freetime and not self.freemass):
            ceq = np.hstack((drf, dvf, [dmf], [H]))
        elif (not self.freemass and not self.freetime):
            ceq = np.hstack((drf, dvf, [dmf]))
        else:
            raise AttributeError("Could not determine equality constraint vector.")

        return ceq

    def get_states(self, atol=1e-12, rtol=1e-12):
        """
        """

        # traj = [t, x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz, lm, u, ux, uy, uz]
        # traj.shape = (npts, 19)

        # propagate trajectory
        self._propagate(atol=atol, rtol=rtol)

        # get times
        t = self._times.reshape(self._times.size, 1)
        # redimensionalise times
        t *= self._dynamics.T

        # get controls
        u = np.asarray([self._dynamics._pontryagin(fs)
                        for fs in self._trajectory])

        # get Hamiltonian
        H = np.asarray([self._dynamics._hamiltonian(fs)
                        for fs in self._trajectory])
        H = H.reshape(H.size, 1)

        # get trajectory
        traj = self._trajectory
        # redimensionalise trajectory
        traj[:, 0:3] *= self._dynamics.L
        traj[:, 3:6] *= self._dynamics.V
        traj[:, 6]   *= self._dynamics.M

        # assemble full trajectory history
        traj = np.hstack((t, traj, u, H))

        return traj


if __name__ == "__main__":

    # create spacecraft
    sc = pk.sims_flanagan.spacecraft(1000, 0.3, 2500)

    # create boundaries
    t0 = pk.epoch(0)
    tf = pk.epoch(2000)

    # create planets
    p0 = pk.planet.jpl_lp("earth")
    pf = pk.planet.jpl_lp("mars")

    # get states
    x0, v0 = p0.eph(t0)
    xf, vf = pf.eph(tf)

    # create spacecraft state
    x0 = pk.sims_flanagan.sc_state(x0, v0, sc.mass)
    xf = pk.sims_flanagan.sc_state(xf, vf, sc.mass / 10)

    # departure costate
    l0 = np.random.randn(7)

    # case 1: default constructor
    l = leg()
    l.set(t0, x0, l0, tf, xf)
    l.mismatch_constraints()
    print(l.get_states())
