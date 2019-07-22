from numpy.linalg import norm
from math import sqrt, asin, acos

from pykep.core import epoch, fb_con, EARTH_VELOCITY, AU, MU_SUN
from pykep.planet import jpl_lp
from pykep.sims_flanagan import leg, spacecraft, sc_state

class mga_lt_nep:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) problem representing a low-thrust
    interplanetary trajectory modelled as a Multiple Gravity Assist trajectory using multiple :class:`pykep.sims_flanagan_leg` 

    - Yam, C.H., di Lorenzo, D., and Izzo, D., Low-Thrust Trajectory Design as a Constrained Global Optimization Problem, Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering, 225(11), pp.1243-1251, 2011.

    The decision vector (chromosome) is::

      [t0] + 
      [T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1] + 
      [T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2] + 
      ... + 
      [throttles1] + 
      [throttles2] + 
      ...

    .. note::

      The resulting problem is non linearly constrained. 
    """

    def __init__(self,
                 seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('mercury')],
                 n_seg = [5, 20],
                 t0 = [3000, 4000],
                 tof = [[100, 1000], [200, 2000]],
                 vinf_dep = 3,
                 vinf_arr = 2,
                 mass = 2000.0,
                 Tmax = 0.5,
                 Isp = 3500.0,
                 fb_rel_vel = 6,
                 multi_objective = False,
                 high_fidelity = False):
        """
        Args::

            - seq (```list of pykep.planet```): defines the encounter sequence for the trajectoty (including the initial planet).
            - n_seg (```list``` of ```int```): the number of segments to be used for each leg.
            - t0 (```list``` of ```floats```): the launch window (in mjd2000).
            - tof (```list``` of ```2D-list```): minimum and maximum time of each leg (days).
            - vinf_dep (```float```): maximum launch hyperbolic velocity allowed (in km/sec).
            - vinf_arr (```float```): maximum arrival hyperbolic velocity allowed (in km/sec).
            - mass (```float```): spacecraft starting mass. (in kg).
            - Tmax (```float```): maximum thrust (in N).
            - Isp (```float```):: engine specific impulse (in s).
.           - fb_rel_vel (```float```): determines the bounds on the maximum allowed relative velocity at all fly-bys (in km/sec).
            - multi-objective (```bool```): when True defines the problem as a multi-objective problem, returning total DV and time of flight.
            - high_fidelity (```bool```):makes the trajectory computations slower, but actually dynamically feasible.
        """
        # We define some data members (we use the double underscore to
        # indicate they are private)
        self._seq = seq
        self._n_seg = n_seg
        self._t0 = t0
        self._tof = tof
        self._vinf_dep = vinf_dep * 1000
        self._vinf_arr = vinf_arr * 1000
        self._mass = mass
        self._Tmax = Tmax
        self._fb_rel_vel = fb_rel_vel
        self._multiobjective = multi_objective
        self._high_fidelity = high_fidelity

        self._n_legs = len(seq) - 1
        self._sc = spacecraft(mass, Tmax, Isp)
        self._leg = leg()
        self._leg.set_mu(MU_SUN)
        self._leg.set_spacecraft(self._sc)

    def get_bounds(self):
        # Convenience aliases
        t0 = self._t0
        mass = self._mass
        fb_rel_vel = self._fb_rel_vel
        n_legs = self._n_legs
        n_seg = self._n_seg
        vinf_dep = self._vinf_dep
        vinf_arr = self._vinf_arr
        tof = self._tof

        # Basic bounds
        lb = [t0[0]] + [0., mass / 2., -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -
                                fb_rel_vel, -fb_rel_vel, -fb_rel_vel] * n_legs + [-1, -1, -1] * sum(n_seg)
        ub = [t0[1]] + [1, mass, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel,
                                fb_rel_vel, fb_rel_vel] * n_legs + [1, 1, 1] * sum(n_seg)
        # bounds on the vinfs......
        lb[3:6] = [-vinf_dep] * 3
        ub[3:6] = [vinf_dep] * 3
        lb[-sum(n_seg) * 3 - 3:-sum(n_seg) * 3] = [-vinf_arr] * 3
        ub[-sum(n_seg) * 3 - 3:-sum(n_seg) * 3] = [vinf_arr] * 3
        # and for the time of flight
        lb[1:1 + 8 * n_legs:8] = [el[0] for el in tof]
        ub[1:1 + 8 * n_legs:8] = [el[1] for el in tof]

        return (lb, ub)

    def fitness(self, x):
        # The final mass is the fitness
        objfun = [x[2 + 8 * (self._n_legs - 1)]]
        eq_c = []
        ineq_c = []
        if self._multiobjective:
            objfun = objfun + [sum(x[1:8*self._n_legs:8])]
        # We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self._n_legs + 1))
        r_P = list([None] * (self._n_legs + 1))
        v_P = list([None] * (self._n_legs + 1))
        for i in range(len(self._seq)):
            t_P[i] = epoch(x[0] + sum(x[1:i*8:8]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])
        # We assemble the constraints. 
        # 1 - Mismatch Constraints
        for i in range(self._n_legs):
            # Departure velocity of the spacecraft in the heliocentric frame
            v0 = [a + b for a, b in zip(v_P[i], x[3 + 8 * i:6 + 8 * i])]
            if i==0:
                m0 = self._mass
            else:
                m0 = x[2 + 8 * (i-1)]
            x0 = sc_state(r-P[i], v0, m0)
            vf = [a + b for a, b in zip(v_P[i+1], x[6 + 8 * i:9 + 8 * i])]
            xf = sc_state(r_P[i+1], vf, x[2 + 8 * i])
            self._leg.set(t_P[i], x0, x[-3 * sum(self._n_seg[i:]):-sum(self._n_seg[i+1:]) * 3], tf, xf)
            mismatch = list(self._leg.mismatch_constraints())
            # Making the mismatch non dimensional (assumes an heliocentric interplanetary trajectory)
            mismatch[0] /= AU
            mismatch[1] /= AU
            mismatch[2] /= AU
            mismatch[3] /= EARTH_VELOCITY
            mismatch[4] /= EARTH_VELOCITY
            mismatch[5] /= EARTH_VELOCITY
            mismatch[6] /= self._mass
            eq_c = eq_c + mismatch

        # 2 - Fly-by constraints
        for i in range(self._n_legs - 1):
            DV_eq, alpha_ineq = fb_con(x[6 + i * 8:9 + i * 8], x[11 + i * 8:14 + i * 8], self._seq[i+1])
            eq_c = eq_c + [DV_eq / (EARTH_VELOCITY * EARTH_VELOCITY)]
            ineq_c = ineq_c + [alpha_ineq]

        # 3 - Departure and arrival Vinf
        # departure
        v_dep_con = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] -
                     self._vinf_dep * self._vinf_dep) / (EARTH_VELOCITY * EARTH_VELOCITY)
        # arrival
        n_fb = self._n_legs - 1
        v_arr_con = (x[6 + n_fb * 8] * x[6 + n_fb * 8] + x[7 + n_fb * 8] * x[7 + n_fb * 8] + x[8 + n_fb * 8] * x[8 + n_fb * 8] -
                     self._vinf_arr * self._vinf_arr) / (EARTH_VELOCITY * EARTH_VELOCITY)
        ineq_c = ineq_c + [v_dep_con, v_arr_con]

        return objfun + eq_c + ineq_c

    # And this helps visualizing the trajectory
    def plot(self, x, ax=None):
        """
        ax = prob.plot(x, ax=None)

        - x: encoded trajectory
        - ax: matplotlib axis where to plot. If None figure and axis will be created
        - [out] ax: matplotlib axis where to plot

        Plots the trajectory represented by a decision vector x on the 3d axis ax

        Example::

          ax = prob.plot(x)
        """
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from pykep import epoch, AU
        from pykep.sims_flanagan import sc_state
        from pykep.orbit_plots import plot_planet, plot_sf_leg

        # Creating the axis if necessary
        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.gca(projection='3d')
        else:
            axis = ax

        # Plotting the Sun ........
        axis.scatter([0], [0], [0], color='y')

        # Plotting the legs .......
        # 1 - We decode the chromosome extracting the time of flights
        T = list([0] * (self.__n_legs))
        for i in range(self.__n_legs):
            T[i] = x[1 + i * 8]

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 1))
        r_P = list([None] * (self.__n_legs + 1))
        v_P = list([None] * (self.__n_legs + 1))

        for i, planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq[i].eph(t_P[i])

        # 3 - We iterate through legs to compute mismatches and throttles
        # constraints
        ceq = list()
        cineq = list()
        m0 = self.__sc.mass
        for i in range(self.__n_legs):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(3 + i * 8):(6 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a,
                 b in zip(v_P[i + 1], x[(6 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[2 + i * 8])
            throttles = x[(1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(
                1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            m0 = x[2 + 8 * i]
            plot_sf_leg(self.__leg, units=AU, N=10, axes=axis)

        # Plotting planets
        for i, planet in enumerate(self.seq):
            plot_planet(planet, t_P[i], units=AU,
                        legend=True, color=(0.7, 0.7, 1), axes=axis)

        plt.show()
        return axis

    def high_fidelity(self, boolean):
        """
        prob.high_fidelity(status)

        - status: either True or False (True sets high fidelity on)

        Sets the trajectory high fidelity mode

        Example::

          prob.high_fidelity(True)
        """
        # We avoid here that objfun and constraint are kept that have been
        # evaluated wrt a different fidelity
        self.reset_caches()
        # We set the propagation fidelity
        self.__leg.high_fidelity = boolean

    def ic_from_mga_1dsm(self, x):
        """
        x_lt = prob.ic_from_mga_1dsm(x_mga)

        - x_mga: compatible trajectory as encoded by an mga_1dsm problem

        Returns an initial guess for the low-thrust trajectory, converting the mga_1dsm solution x_dsm. The user
        is responsible that x_mga makes sense (i.e. it is a viable mga_1dsm representation). The conversion is done by importing in the
        low-thrust encoding a) the launch date b) all the legs durations, c) the in and out relative velocities at each planet.
        All throttles are put to zero.

        Example::

          x_lt= prob.ic_from_mga_1dsm(x_mga)
        """
        from math import pi, cos, sin, acos
        from scipy.linalg import norm
        from pykep import propagate_lagrangian, lambert_problem, DAY2SEC, fb_prop

        retval = list([0.0] * self.dimension)
        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T
        T = list([0] * (self.__n_legs))

        for i in range(len(T)):
            T[i] = log(x[2 + 4 * i])
        total = sum(T)
        T = [x[1] * time / total for time in T]

        retval[0] = x[0]
        for i in range(self.__n_legs):
            retval[1 + 8 * i] = T[i]
            retval[2 + 8 * i] = self.__sc.mass

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 1))
        r_P = list([None] * (self.__n_legs + 1))
        v_P = list([None] * (self.__n_legs + 1))
        DV = list([None] * (self.__n_legs + 1))

        for i, planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        theta = 2 * pi * x[1]
        phi = acos(2 * x[2] - 1) - pi / 2

        Vinfx = x[3] * cos(phi) * cos(theta)
        Vinfy = x[3] * cos(phi) * sin(theta)
        Vinfz = x[3] * sin(phi)

        retval[3:6] = [Vinfx, Vinfy, Vinfz]

        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(r_P[0], v0, x[4] * T[0] * DAY2SEC, MU_SUN)

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * DAY2SEC
        l = lambert_problem(r, r_P[1], dt, MU_SUN)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        retval[6:9] = [a - b for a, b in zip(v_end_l, v_P[1])]

        # 4 - And we proceed with each successive leg
        for i in range(1, self.__n_legs):
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            7 + (i - 1) * 4] * self.seq[i].radius, x[6 + (i - 1) * 4], self.seq[i].mu_self)
            retval[3 + i * 8:6 + i * 8] = [a -
                                           b for a, b in zip(v_out, v_P[i])]
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, MU_SUN)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem(r, r_P[i + 1], dt, MU_SUN)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
            retval[6 + i * 8:9 + i * 8] = [a -
                                           b for a, b in zip(v_end_l, v_P[i + 1])]
        return retval

    def double_segments(self, x):
        """
        x_doubled = prob.double_segments(x)

        - x: compatible trajectory as encoded by an mga_1dsm mga_lt_nep

        Returns the decision vector encoding a low trust trajectory having double the number of segments with respect to x
        and a 'similar' throttle history. In case high fidelity is True, and x is a feasible trajectory, the returned decision vector
        also encodes a feasible trajectory that can be further optimized

        Example::

          prob = traj.mga_lt_nep(nseg=[[10],[20]])
          pop = population(prob,1)
          .......OPTIMIZE.......
          x = prob.double_segments(pop.champion.x)
          prob = traj.mga_lt_nep(nseg=[[20],[40]])
          pop = population(prob)
          pop.push_back(x)
          .......OPTIMIZE AGAIN......
        """
        y = list()
        y.extend(x[:-sum(self.__n_seg) * 3])
        for i in range(sum(self.__n_seg)):
            y.extend(x[-(sum(self.__n_seg) - i) * 3:-
                       (sum(self.__n_seg) - 1 - i) * 3] * 2)
        y.extend(x[-3:] * 2)
        return y
