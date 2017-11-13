from PyGMO.problem import base as base_problem
from pykep.core import epoch, fb_con, EARTH_VELOCITY, AU, MU_SUN
from pykep.planet import jpl_lp
from pykep.sims_flanagan._sims_flanagan import leg, spacecraft, sc_state


class mga_lt_nep(base_problem):

    """
    This class is a PyGMO (http://esa.github.io/pygmo/) problem representing a low-thrust
    interplanetary trajectory modelled as a Multiple Gravity Assist trajectory with sims_flanagan legs

    - Yam, C.H., di Lorenzo, D., and Izzo, D.,	 Low-Thrust Trajectory Design as a Constrained Global Optimization Problem,  Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering, 225(11), pp.1243-1251, 2011.

    The decision vector (chromosome) is::

      [t0] + [T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1] + [T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2] + ... + [throttles1] + [throttles2] + ...

    .. note::

      The resulting problem is non linearly constrained. The resulting trajectory is not time-bounded.
    """

    def __init__(self,
                 seq=[jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')],
                 n_seg=[10] * 2,
                 t0=[epoch(0), epoch(1000)],
                 tof=[[200, 500], [200, 500]],
                 vinf_dep=2.5,
                 vinf_arr=2.0,
                 mass=4000.0,
                 Tmax=1.0,
                 Isp=2000.0,
                 fb_rel_vel=6,
                 multi_objective=False,
                 high_fidelity=False):
        """
        prob = mga_lt_nep(seq = [jpl_lp('earth'),jpl_lp('venus'),jpl_lp('earth')], n_seg = [10]*2,
        t0 = [epoch(0),epoch(1000)], tof = [[200,500],[200,500]], vinf_dep=2.5, vinf_arr=2.0, mass=4000.0, Tmax=1.0, Isp=2000.0,
        multi_objective = False, fb_rel_vel = 6, high_fidelity=False)

        - seq: list of pykep.planet defining the encounter sequence for the trajectoty (including the initial planet)
        - n_seg: list of integers containing the number of segments to be used for each leg (len(n_seg) = len(seq)-1)
        - t0: list of pykep epochs defining the launch window
        - tof: minimum and maximum time of each leg (days)
        - vinf_dep: maximum launch hyperbolic velocity allowed (in km/sec)
        - vinf_arr: maximum arrival hyperbolic velocity allowed (in km/sec)
        - mass: spacecraft starting mass
        - Tmax: maximum thrust
        - Isp: engine specific impulse
        - fb_rel_vel = determines the bounds on the maximum allowed relative velocity at all fly-bys (in km/sec)
        - multi-objective: when True defines the problem as a multi-objective problem, returning total DV and time of flight
        - high_fidelity = makes the trajectory computations slower, but actually dynamically feasible.
        """

        # 1) We compute the problem dimensions .... and call the base problem constructor
        self.__n_legs = len(seq) - 1
        n_fb = self.__n_legs - 1
        # 1a) The decision vector length
        dim = 1 + self.__n_legs * 8 + sum(n_seg) * 3
        # 1b) The total number of constraints (mismatch + fly-by + boundary + throttles
        c_dim = self.__n_legs * 7 + n_fb * 2 + 2 + sum(n_seg)
        # 1c) The number of inequality constraints (boundary + fly-by angle + throttles)
        c_ineq_dim = 2 + n_fb + sum(n_seg)
        # 1d) the number of objectives
        f_dim = multi_objective + 1
        # First we call the constructor for the base PyGMO problem
        # As our problem is n dimensional, box-bounded (may be multi-objective), we write
        # (dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
        super(mga_lt_nep, self).__init__(dim, 0, f_dim, c_dim, c_ineq_dim, 1e-4)

        # 2) We then define some class data members
        # public:
        self.seq = seq
        # private:
        self.__n_seg = n_seg
        self.__vinf_dep = vinf_dep * 1000
        self.__vinf_arr = vinf_arr * 1000
        self.__sc = spacecraft(mass, Tmax, Isp)
        self.__leg = leg()
        self.__leg.set_mu(MU_SUN)
        self.__leg.set_spacecraft(self.__sc)
        self.__leg.high_fidelity = high_fidelity
        fb_rel_vel *= 1000
        # 3) We compute the bounds
        lb = [t0[0].mjd2000] + [0, mass / 2, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel] * self.__n_legs + [-1, -1, -1] * sum(self.__n_seg)
        ub = [t0[1].mjd2000] + [1, mass, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel] * self.__n_legs + [1, 1, 1] * sum(self.__n_seg)
        # 3a ... and account for the bounds on the vinfs......
        lb[3:6] = [-self.__vinf_dep] * 3
        ub[3:6] = [self.__vinf_dep] * 3
        lb[-sum(self.__n_seg) * 3 - 3:-sum(self.__n_seg) * 3] = [-self.__vinf_arr] * 3
        ub[-sum(self.__n_seg) * 3 - 3:-sum(self.__n_seg) * 3] = [self.__vinf_arr] * 3
        # 3b... and for the time of flight
        lb[1:1 + 8 * self.__n_legs:8] = [el[0] for el in tof]
        ub[1:1 + 8 * self.__n_legs:8] = [el[1] for el in tof]

        # 4) And we set the bounds
        self.set_bounds(lb, ub)

    # Objective function
    def _objfun_impl(self, x):
        if self.f_dimension == 1:
            return (-x[2 + (self.__n_legs - 1) * 8],)
        else:
            return (-x[2 + (self.__n_legs - 1) * 8], sum(x[1:1 + 8 * self.__n_legs:8]))

    # Constraints function
    def _compute_constraints_impl(self, x):
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

        # 3 - We iterate through legs to compute mismatches and throttles constraints
        ceq = list()
        cineq = list()
        m0 = self.__sc.mass
        for i in range(self.__n_legs):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(3 + i * 8):(6 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a, b in zip(v_P[i + 1], x[(6 + i * 8):(9 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[2 + i * 8])
            throttles = x[(1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            m0 = x[2 + 8 * i]
            ceq.extend(self.__leg.mismatch_constraints())
            cineq.extend(self.__leg.throttles_constraints())

        # Adding the boundary constraints
        # departure
        v_dep_con = (x[3] ** 2 + x[4] ** 2 + x[5] ** 2 - self.__vinf_dep ** 2) / (EARTH_VELOCITY ** 2)
        # arrival
        v_arr_con = (x[6 + (self.__n_legs - 1) * 8] ** 2 + x[7 + (self.__n_legs - 1) * 8] ** 2 + x[8 + (self.__n_legs - 1) * 8] ** 2 - self.__vinf_arr ** 2) / (EARTH_VELOCITY ** 2)
        cineq.append(v_dep_con * 100)
        cineq.append(v_arr_con * 100)

        # We add the fly-by constraints
        for i in range(self.__n_legs - 1):
            DV_eq, alpha_ineq = fb_con(x[6 + i * 8:9 + i * 8], x[11 + i * 8:14 + i * 8], self.seq[i + 1])
            ceq.append(DV_eq / (EARTH_VELOCITY ** 2))
            cineq.append(alpha_ineq)

        # Making the mismatches non dimensional
        for i in range(self.__n_legs):
            ceq[0 + i * 7] /= AU
            ceq[1 + i * 7] /= AU
            ceq[2 + i * 7] /= AU
            ceq[3 + i * 7] /= EARTH_VELOCITY
            ceq[4 + i * 7] /= EARTH_VELOCITY
            ceq[5 + i * 7] /= EARTH_VELOCITY
            ceq[6 + i * 7] /= self.__sc.mass

        # We assemble the constraint vector
        retval = list()
        retval.extend(ceq)
        retval.extend(cineq)

        return retval

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

        # 3 - We iterate through legs to compute mismatches and throttles constraints
        ceq = list()
        cineq = list()
        m0 = self.__sc.mass
        for i in range(self.__n_legs):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(3 + i * 8):(6 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a, b in zip(v_P[i + 1], x[(6 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[2 + i * 8])
            throttles = x[(1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(1 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            m0 = x[2 + 8 * i]
            plot_sf_leg(self.__leg, units=AU, N=10, ax=axis)

        # Plotting planets
        for i, planet in enumerate(self.seq):
            plot_planet(planet, t_P[i], units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)

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
        # We avoid here that objfun and constraint are kept that have been evaluated wrt a different fidelity
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
        # 1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
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
            v_out = fb_prop(v_end_l, v_P[i], x[7 + (i - 1) * 4] * self.seq[i].radius, x[6 + (i - 1) * 4], self.seq[i].mu_self)
            retval[3 + i * 8:6 + i * 8] = [a - b for a, b in zip(v_out, v_P[i])]
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, MU_SUN)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem(r, r_P[i + 1], dt, MU_SUN)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
            retval[6 + i * 8:9 + i * 8] = [a - b for a, b in zip(v_end_l, v_P[i + 1])]
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
            y.extend(x[-(sum(self.__n_seg) - i) * 3:-(sum(self.__n_seg) - 1 - i) * 3] * 2)
        y.extend(x[-3:] * 2)
        return y
