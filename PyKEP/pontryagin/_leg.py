from PyKEP.pontryagin._dynamics import _dynamics
from PyKEP.core import MU_SUN, epoch, AU
from PyKEP.sims_flanagan import spacecraft, sc_state
from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


class leg(object):
    """Indirect optimal control transcription trajectory leg.

    This class represents an indirect optimal control transcription
    (alla moda di `Pontryagin's maximum principle <https://en.wikipedia.org/wiki/Pontryagin%27s_maximum_principle>`_) of a generic
    two-point boundary trajectory.

    Attributes:
        - t0 (``float``): Departure time [mjd2000].
        - x0 (``numpy.ndarray``): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
        - l0 (``numpy.ndarray``): Costate variables [ND, ND, ND, ND, ND, ND, ND].
        - tf (``float``): Arrival time [mjd2000].
        - xf (``numpy.ndarray``): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].
        - sc (``PyKEP.sims_flanagan.spacecraft``): Generic spacecraft with propulsive properties.
        - mu (``float``): Gravitational parametre of primary body [m^3/s^2].
        - freemass (``bool``): Activates final mass transversality condition.
        - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
        - alpha (``float``): Homotopy parametre.
        - bound (``bool``): Activates bounded control.
        - nec (``int``): Number of equality constraints.
        - trajectory (``numpy.ndarray``): Array of nondimensional fullstates.
        - times (``numpy.ndarray``): Vector of nondimensional integration times.
    """

    def __init__(self, t0=None, x0=None, l0=None, tf=None, xf=None, sc=spacecraft(1000, 0.3, 2500), mu=MU_SUN, freemass=True, freetime=True, alpha=1, bound=True):
        """Initialises an indirect optimal control transcription trajectory leg.

        Args:
            - t0 (``PyKEP.epoch``, ``None``): Departure time [mjd2000].
            - x0 (``PyKEP.sims_flanagan.sc_state``, ``None``): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
            - l0 (``numpy.ndarray``, ``list``, ``tuple``, ``None``): Costate variables [ND, ND, ND, ND, ND, ND, ND].
            - tf (``PyKEP.epoch``, ``None``): Arrival time [mjd2000].
            - xf (``PyKEP.sims_flanagan.sc_state``, ``None``): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].
            - sc (``PyKEP.sims_flanagan.spacecraft``): Generic spacecraft with propulsive properties.
            - mu (``float``, ``int``): Gravitational parametre of primary body [m^3/s^2].
            - freemass (``bool``): Activates final mass transversality condition.
            - freetime (``bool``): Activates final time transversality condition. Allows final time to vary.
            - alpha (``float``, ``int``): Homotopy parametre, governing the degree to which the theoretical control law is intended to reduce propellant expenditure or energy. Setting the parametre to 1 enforces a mass-optimal control law, with a characteristic bang-bang control profile (either full throttle or off). Setting the parametre to 0 enforces a pure quadratic control law.
            - bound (``bool``): Activates bounded control, in which the control throttle is bounded between 0 and 1, otherwise the control throttle is allowed to unbounded.

        Raises:
            - TypeError: If ``sc`` is not supplied as an instance of `PyKEP.sims_flanagan.spacecraft`.
            - TypeError: If ``mu`` is neither supplied as an instance of ``float`` nor ``int``.
            - ValueError: If ``mu`` is not supplied as a positive number.
            - TypeError: If either ``freemass``, ``freetime``, or ``bound`` is not supplied as an instance of ``bool``.
            - AttributeError: If equality constraint dimension cannot be determined form ``freemass`` and ``freetime``.
            - TypeError: If ``alpha`` is neither supplied as an instance of ``float`` nor ``int``.
            - ValueError: If ``alpha`` is not in between 0 and 1.
            - ValueError: If ``alpha == 1`` and ``bound == False``. Control cannot be unbounded with mass-optimal control.
            - TypeError: If either ``t0`` or ``tf`` is not supplied as an instance of ``PyKEP.epoch``.
            - ValueError: If departure time, ``t0.mjd2000`` does not occur before arrival time, `t`f.mjd2000``.
            - TypeError: If either ``x0`` or ``xf`` is not supplied as an instance of ``PyKEP.sims_flanagan.sc_state``.
            - TypeError: If costate, ``l0``, is neither supplied as an instance of ``numpy.ndarray``, ``list``, nor ``tuple``.
            - TypeError: If costate, ``l0``, is does not have 7 elements. Each element of ``l0`` corresponds to the respective elements of ``x0``.

        Examples:
            >>> l = PyKEP.pontryagin.leg()
            >>> l = PyKEP.pontryagin.leg(t0, x0, l0, tf, xf)

        .. note::
            Typically spacecraft trajectories are optimised to reduce final
            mass, which theoretically results in a bang-bang control profile.
            However, a bang-bang control profile (either 0 or 1) is problematic
            to an optimiser due to its discontinous nature. The trajectory
            optimisation problem becomes easier when one first optimises with
            unbounded quadratic control (``alpha == 0 and bound == False``),
            stores the solution, then uses the solution to optimise with
            bounded quadratic control (``alpha == 0 and bound == True``). Then
            using the solution from bounded quadratic control, one can
            incrementally optimise for increasing homotopy parametres
            (e.g. ``alphas = numpy.linspace(0, 1, 200)``) until a mass-optimal
            control solution converges.

        .. note::
            Boundary conditions ``t0``, ``x0``, ``l0``, ``tf``, and ``xf``
            are optional in the constructor. If some, but not all, boundary
            conditions are supplied, they are not set and disregarded. If the
            boundary conditions are not set upon construction (``__init__``),
            they must be set with ``set(t0, x0, l0, tf, xf)``.

        .. note::
            ``mismatch_constraints`` will not work unless the boundary
            conditions have been set either through ``__init__`` or ``set``.
        """

        # check spacecraft
        if isinstance(sc, spacecraft):
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

        # dynamics
        self._dynamics = _dynamics(sc=self.spacecraft, mu=mu, alpha=alpha, bound=bound)

        # check freetime, freemass, and bound
        if not all([isinstance(param, bool) for param in [freemass, freetime, bound]]):
            raise TypeError(
                "Freemass, freetime, bound parametres must supplied as booleans.")
        else:
            self.freemass = bool(freemass)
            self.freetime = bool(freetime)
            self.bound = bool(bound)

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
            raise AttributeError(
                "Could not determine equality constraint dimensionality.")

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
            if not all([isinstance(t, epoch) for t in [t0, tf]]):
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
            if not all([isinstance(state, sc_state) for state in states]):
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

        # integrator
        self._integrator = ode(
            lambda t, fs: self._dynamics._eom_fullstate(fs),
            #lambda t, fs: self._dynamics._eom_fullstate_jac(fs)
        )

    def _recorder(self, t, fs):

        # append time
        self.times = np.append(self.times, t)

        # append fullstate
        self.trajectory = np.vstack((self.trajectory, fs))

    def set(self, t0, x0, l0, tf, xf):
        """Sets the departure and arrival boundary conditions of the trajectory.

        This instance method sets the departure time `t0`, departure state `x0`,
        departure costate `l0`, arrival time `tf`, arrival state `xf`.

        Args:
            - t0 (``PyKEP.epoch``): Departure time [mjd2000].
            - x0 (``PyKEP.sims_flanagan.sc_state``): Cartesian departure state [m, m, m, m/s, m/s, m/s, kg].
            - l0 (``numpy.ndarray``, ``list``, ``tuple``): Costate variables [ND, ND, ND, ND, ND, ND, ND].
            - tf (``PyKEP.epoch``): Arrival time [mjd2000].
            - xf (``PyKEP.epoch``): Cartesian arrival state [m, m, m, m/s, m/s, m/s, kg].

        Raises:
            - TypeError: If either ``t0`` or ``tf`` is not supplied as an instance of ``PyKEP.epoch``.
            - TypeError: If either ``x0`` or ``xf`` is not supllied as an instance of ``PyKEP.sims_flanagan.sc_state``.
            - TypeError: If ``l0`` is neither supllied as a ``numpy.ndarray``, ``list``, nor ``tuple``.

        .. note::
            If ``args`` were not supplied in ``__init__``, ``set`` must be called before calling ``mismatch_constraints``.

        Examples:
            >>> l = leg() # t0, x0, l0, tf, xf not set yet
            >>> l.set(t0, x0, l0, tf, xf) # now they're set
        """

        # check inputs
        if not all([isinstance(t, epoch) for t in [t0, tf]]):
            raise TypeError("t0 and tf must be instances of PyKEP.epoch.")
        elif not t0.mjd2000 < tf.mjd2000:
            raise ValueError("t0 must be less than tf.")
        elif not all([isinstance(state, sc_state) for state in [x0, xf]]):
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
        x0       = np.copy(self.x0) #NOTE: very important
        x0[0:3] /= self._dynamics.L
        x0[3:6] /= self._dynamics.V
        x0[6]   /= self._dynamics.M

        # nondimensional fullstate
        fs0 = np.hstack((x0, self.l0, [0]))

        # convert mjd2000 to mjs2000
        t0 = self.t0 * 24 * 60 * 60
        tf = self.tf * 24 * 60 * 60

        # nondimensionalise times
        t0 /= self._dynamics.T
        tf /= self._dynamics.T

        # clear trajectory history
        self.times      = np.empty((1, 0), dtype=np.float64)
        self.trajectory = np.empty((0, 15), dtype=np.float64)

        # set integration method
        self._integrator.set_integrator("dop853", atol=atol, rtol=rtol)

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
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.

        Returns:
            - ``numpy.ndarray``: Equality constraints vector ``ceq`` composed of the arrival mismatch in position & velocity, arrival mass costate if ``freemass == True``, arrival mismatch in mass if ``freemass == False``, and arrival Hamiltonian if ``freetime == True``.
            ::

                ceq = [drf, dvf, dmf]    # freemass == False; freetime == False
                ceq = [drf, dvf, lmf]    # freemass == True; freetime == False
                ceq = [drf, dvf, lmf, H] # freemass == True; freetime == True

        Raises:
            - TypeError: If either ``atol`` or ``rtol`` is supplied as neither an instance of ``float`` nor ``int``.
            - AttributeError: If boundary conditions ``t0``, ``x0``, ``l0``, ``tf``, and ``x0`` have not been set through either ``__init__`` or ``set``.

        .. note::
            This method uses an explicit Runge-Kutta method of order (4)5
            due to Dormand & Prince with adaptive stepsize control.
            Smaller values of ``atol`` and ``rtol`` will increase the accuracy of
            a converged trajectory optimisation solution. However, smaller
            values result in slower integration executions and thus slower
            optimiser iterations. It is preferable to first solve a trajectory
            optimisation problem with higher values of ``atol`` and ``rtol`` to
            quickly converge to an approximate solution, subsequently resolving
            the problem with smaller tolerance values using the previous solution.
            This method is akin to mesh refinement in direct trajectory optimisation.

        .. note::
            State arrival boundary conditions in position and velocity are
            always returned. If the final mass transversality condition is
            activated (``freemass == True``), than the final mass costate is
            returned as part of the equality constraint vector. Otherwise
            (``freemass == False``), the final mass mismatch is returned as part
            of the equality constraint vector. If the final time transversality
            condition is activated ``freetime == True``, then than the final
            Hamiltonian will be returned as part of the equality constraint
            vector. Otherwise (``freetime == False``), the Hamiltonian will not
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
            raise TypeError(
                "Both atol and rtol must be supplied as instances of either float or int.")
        if any([hasattr(self, atr) == False for atr in ["t0", "x0", "l0", "tf", "xf"]]):
            raise AttributeError(
                "Cannot propagate dynamics, as boundary conditions t0, x0, l0, tf, and xf have not been set. Use set(t0, x0, l0, tf, xf) to set boundary conditions.")
        else:
            atol = float(atol)
            rtol = float(atol)

        # propagate trajectory
        self._propagate(atol, rtol)

        # desired nondimensional arrival states
        brf = self.xf[0:3] / self._dynamics.L
        bvf = self.xf[3:6] / self._dynamics.V

        # propagated nondimensional arrival states
        rf = self.trajectory[-1, 0:3]
        vf = self.trajectory[-1, 3:6]

        # nondimensional position and velocity arrival mismatch
        drf = rf - brf
        dvf = vf - bvf

        # free arrival mass
        if self.freemass:
            lmf = self.trajectory[-1, 13]
        # fixed arrival mass
        else:
            dmf = self.trajectory[-1, 6] - self.xf[6] / self._dynamics.M

        # free arrival time
        if self.freetime:
            Hf = self._dynamics._hamiltonian(self.trajectory[-1])
        # fixed arrival time
        else:
            pass

        # create equality constraints
        if (self.freemass and self.freetime):
            ceq = np.hstack((drf, dvf, [lmf], [Hf]))
        elif (self.freemass and not self.freetime):
            ceq = np.hstack((drf, dvf, [lmf]))
        elif (self.freetime and not self.freemass):
            ceq = np.hstack((drf, dvf, [dmf], [Hf]))
        elif (not self.freemass and not self.freetime):
            ceq = np.hstack((drf, dvf, [dmf]))
        else:
            raise AttributeError(
                "Could not determine equality constraint vector.")

        return ceq

    def get_states(self, atol=1e-12, rtol=1e-12):
        """Returns the trajectory data of time, states, costates, and controls.

        This method propagates the spacecrafts dynamics with a chosen ``atol``
        ``rtol`` from ``t0``, ``x0``, ``l0`` to ``t0``, then returns a `numpy.ndarray`
        with the number of rows corresponding to the number of data points
        along the trajectory (as determinded by ``atol`` and ``rtol``), and 19
        columns.

        Args:
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.

        Returns:
            - ``numpy.ndarray``: Trajectory data array of shape ``(npts, 20)``, where
            ``npts`` is the number of data points along the trajectory, as
            determind by the integration error tolerances ``atol`` and ``rtol``.

            The returned array is characterised by:
            ::

                [[t0, x0, y0, z0, vx0, vy0, vz0, m0, lx0, ly0, lz0, lvx0, lvy0, lvz0, lm0, u0, ux0, uy0, uz0, H0],
                ...
                 [tf, xf, yf, zf, vxf, vyf, vzf, mf, lxf, lyf, lzf, lvxf, lvyf, lvzf, lmf, uf, uxf, uyf, uzf, Hf]]

        Raises:
            - TypeError: If either ``atol`` or ``rtol`` is supplied as neither an instance of ``float`` nor ``int``.
            - AttributeError: If boundary conditions ``t0``, ``x0``, ``l0``, ``tf``, and ``x0`` have not been set through either ``__init__`` or ``set``.

        .. note::
            The returned array has the units of
            [mjd2000, m, m, m, m/s, m/s, m/s, kg, ND, ND, ND, ND, ND, ND, ND, ND, ND, ND, ND, ND].

        .. note::
            Setting either ``atol`` or ``rtol`` smaller than ``1e-12`` may result in numerical integration difficulties.

        Examples:
            >>> l = leg(t0, x0, l0, tf, xf)
            >>> traj = l.get_states(atol=1e-12, rtol=1e-12)
            >>> t = traj[:, 0] # times
            >>> r = traj[:, 1:4] # positions
            >>> v = traj[:, 4:7] # velocities
            >>> m = traj[:, 7] # masses
            >>> lm = traj[:, 14] # mass costates
            >>> u = traj[:, 15] # control throttles
            >>> H = traj[:, 19] # Hamiltonians
            >>> plt.plot(t, u) # plot the control throttle history
            >>> plt.show()
        """

        if not all([(isinstance(tol, float) or isinstance(tol, int)) for tol in [atol, rtol]]):
            raise TypeError(
                "Both atol and rtol must be supplied as instances of either float or int.")
        if any([hasattr(self, atr) == False for atr in ["t0", "x0", "l0", "tf", "xf"]]):
            raise AttributeError(
                "Cannot propagate dynamics, as boundary conditions t0, x0, l0, tf, and xf have not been set. Use set(t0, x0, l0, tf, xf) to set boundary conditions.")
        else:
            atol = float(atol)
            rtol = float(atol)

        # propagate trajectory
        self._propagate(atol, rtol)

        # nondimensional times
        t = np.copy(self.times)
        # mjs2000 times
        t *= self._dynamics.T
        # mjd2000 times
        t /= 24 * 60 * 60
        # reshape
        t = self.times.reshape(t.size, 1)

        # controls
        u = np.asarray([self._dynamics._pontryagin(fs)
                        for fs in self.trajectory])

        # get Hamiltonian
        H = np.asarray([self._dynamics._hamiltonian(fs)
                        for fs in self.trajectory])
        H = H.reshape(H.size, 1)

        # get trajectory
        traj = np.copy(self.trajectory)
        # redimensionalise trajectory
        traj[:, 0:3] *= self._dynamics.L
        traj[:, 3:6] *= self._dynamics.V
        traj[:,   6] *= self._dynamics.M

        # assemble full trajectory history
        traj = np.hstack((t, traj, u, H))

        return traj

    def plot_traj(self, axes, mark="k.-", atol=1e-11, rtol=1e-11, units=AU):
        """Plots trajectory onto a 3D axis.

        Args:
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axis onto which to plot the trajectory.
            - mark (``str``): Marker style.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - units (``float``, ``int``): Length unit by which to normalise data.

        Raises:
            - TypeError: If ``axes`` is not an instance of ``mpl_toolkits.mplot3d.Axes3D``.
            - TypeError: If ``mark`` is not an instance of ``str``.
            - TypeError: If either ``atol``, ``rtol``, or ``units`` is neither an instance of ``float`` nor ``int``.

        Examples:
            >>> sc = pk.sims_flanagan.spacecraft(1000, 0.3, 2500) # spacecraft
            >>> p0 = pk.planet.jpl_lp("earth")
            >>> pf = pk.planet.jpl_lp("mars")
            >>> t0 = pk.epoch(0)
            >>> tf = pk.epoch(1000)
            >>> r0, v0 = p0.eph(t0)
            >>> rf, vf = pf.eph(tf)
            >>> x0 = pk.sims_flanagan.sc_state(r0, v0, sc.mass)
            >>> xf = pk.sims_flanagan.sc_state(rf, vf, sc.mass/10)
            >>> l0 = np.random.randn(7)
            >>> l = leg(t0, x0, l0, tf, xf)
            >>> fig = plt.figure()
            >>> axes = fig.gca(projection='3d')
            >>> l.plot(axes)
            >>> plt.show()
        """

        if not isinstance(axes, Axes3D):
            raise TypeError(
                "Axis must be instance of matplotlib.axes._subplots.Axes3DSubplot.")
        elif not isinstance(mark, str):
            raise TypeError("Mark must be instance of string.")
        elif not all((isinstance(par, float) or isinstance(par, int)) for par in [atol, rtol, units]):
            raise TypeError(
                "atol, rtol, and units must be either instances of float or int.")
        else:
            pass

        # get trajectory
        traj = self.get_states(atol=atol, rtol=rtol)[:, 1:4]

        # normalise
        traj /= units

        # plot trajectory
        axes.plot(traj[:, 0], traj[:, 1], traj[:, 2], mark)

        return axes

    def plot(self, x, y, mark="k.-", atol=1e-12, rtol=1e-12, unitsx=1, unitsy=1, xlabel=False, ylabel=False, axes = None):
        """Plots in two dimensions of the leg's trajectory data.
        ::

            keys = ['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'lx', 'ly', 'lz', 'lvx', 'lvy', 'lvz', 'lm', 'u', 'ux', 'uy', 'uz', 'H']

        Args:
            - x (``str``): x-axis dimension in ``keys``.
            - y (``str``): y-axis dimension in ``keys``.
            - mark (``str``): Marker style.
            - atol (``float``, ``int``): Absolute integration solution tolerance.
            - rtol (``float``, ``int``): Relative integration solution tolerance.
            - unitsx (``float``, ``int``): Unit by which to normalise x-axis data.
            - unitsy (``float``, ``int``): Unit by which to normalise y-axis data.
            - xlabel (``str``, ``bool``): x-axis label. If label is ``False``, no label is placed; if ``True``, dimension name is placed.
            - ylabel (``str``, ``bool``): y-axis label. If label is ``False``, no label is placed; if ``True``, dimension name is placed.

        """

        keys = ['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'lx', 'ly', 'lz', 'lvx', 'lvy', 'lvz', 'lm', 'obj', 'u', 'ux', 'uy', 'uz', 'H']

        if not all([isinstance(dim, str) for dim in [x, y, mark]]):
            raise TypeError("x, y, and mark must be supplied as instances of str.")
        elif not all([dim in keys for dim in [x, y]]):
            raise ValueError("Both x and y must be in " + str(keys) + ".")
        elif not all([(isinstance(par, float) or isinstance(par, int)) for par in [atol, rtol, unitsx, unitsy]]):
            raise TypeError("atol, rtol, unitsx, and unitsy must be supplied as an instance of either int or float.")
        elif not all([(isinstance(label, str) or isinstance(label, bool)) for label in [xlabel, ylabel]]):
            raise TypeError("xlabel and ylabel must be supplied an instance of either str or bool.")
        else:

            # get trajectory
            traj = self.get_states(atol=atol, rtol=rtol)

            # get components and normalise
            xi, yi = keys.index(x), keys.index(y)

            # create figure
            if axes is None:
                fig = plt.figure()
                axes = fig.gca()

            # plot
            axes.plot(traj[:, xi]/unitsx, traj[:, yi]/unitsy, mark)

            # labels
            if isinstance(xlabel, str):
                plt.xlabel(xlabel)
            elif xlabel is True:
                plt.xlabel(x)

            if isinstance(ylabel, str):
                plt.ylabel(ylabel)
            elif ylabel is True:
                plt.ylabel(y)

            return axes
