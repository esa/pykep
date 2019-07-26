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

      The resulting optimization problem is non linearly constrained and unfeasible a.e. (almost everywhere)
    """

    def __init__(self,
                 seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('mercury')],
                 n_seg = [5, 20],
                 t0 = [3000, 4000],
                 tof = [[100, 1000], [200, 2000]],
                 vinf_dep = 3,
                 vinf_arr = 2,
                 mass = [1200., 2000.0],
                 Tmax = 0.5,
                 Isp = 3500.0,
                 fb_rel_vel = 6,
                 multi_objective = False,
                 high_fidelity = False):
        """
        Args:
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
        self._sc = spacecraft(mass[1], Tmax, Isp)
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
        lb = [t0[0]] + [0., mass[0], -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -
                                fb_rel_vel, -fb_rel_vel, -fb_rel_vel] * n_legs + [-1, -1, -1] * sum(n_seg)
        ub = [t0[1]] + [1, mass[1], fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel,
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
        objfun = [-x[2 + 8 * (self._n_legs - 1)]]
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
                m0 = self._mass[1]
            else:
                m0 = x[2 + 8 * (i-1)]
            x0 = sc_state(r_P[i], v0, m0)
            vf = [a + b for a, b in zip(v_P[i+1], x[6 + 8 * i:9 + 8 * i])]
            xf = sc_state(r_P[i+1], vf, x[2 + 8 * i])
            idx_start = 1 + 8 * self._n_legs + sum(self._n_seg[:i]) * 3
            idx_end   = 1 + 8 * self._n_legs + sum(self._n_seg[:i+1]) * 3
            self._leg.set(t_P[i], x0, x[idx_start:idx_end], t_P[i+1], xf)
            mismatch = list(self._leg.mismatch_constraints())
            # Making the mismatch non dimensional (assumes an heliocentric interplanetary trajectory)
            mismatch[0] /= AU
            mismatch[1] /= AU
            mismatch[2] /= AU
            mismatch[3] /= EARTH_VELOCITY
            mismatch[4] /= EARTH_VELOCITY
            mismatch[5] /= EARTH_VELOCITY
            mismatch[6] /= self._mass[1]
            eq_c = eq_c + mismatch
            ineq_c = ineq_c + list(self._leg.throttles_constraints())

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

    def get_nec(self):
        return self._n_legs * 7 + (self._n_legs - 1)

    def get_nic(self):
        return sum(self._n_seg) + (self._n_legs - 1) + 2


    # And this helps visualizing the trajectory
    def plot(self, x, axes=None):
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
        from pykep.orbit_plots import plot_sf_leg, plot_planet

        # Creating the axis if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            ax = fig.gca(projection='3d')
        else:
            ax = axes

        # Plotting the Sun ........
        ax.scatter([0], [0], [0], color=['y'])

         # We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self._n_legs + 1))
        r_P = list([None] * (self._n_legs + 1))
        v_P = list([None] * (self._n_legs + 1))
        for i in range(len(self._seq)):
            t_P[i] = epoch(x[0] + sum(x[1:i*8:8]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])
            plot_planet(self._seq[i], t0 = t_P[i],
                        units=AU, legend=False, color=(0.7, 0.7, 0.7), s=30, axes=ax)

        # We assemble the constraints. 
        # 1 - Mismatch Constraints
        for i in range(self._n_legs):
            # Departure velocity of the spacecraft in the heliocentric frame
            v0 = [a + b for a, b in zip(v_P[i], x[3 + 8 * i:6 + 8 * i])]
            if i==0:
                m0 = self._mass[1]
            else:
                m0 = x[2 + 8 * (i-1)]
            x0 = sc_state(r_P[i], v0, m0)
            vf = [a + b for a, b in zip(v_P[i+1], x[6 + 8 * i:9 + 8 * i])]
            xf = sc_state(r_P[i+1], vf, x[2 + 8 * i])
            idx_start = 1 + 8 * self._n_legs + sum(self._n_seg[:i]) * 3
            idx_end   = 1 + 8 * self._n_legs + sum(self._n_seg[:i+1]) * 3
            self._leg.set(t_P[i], x0, x[idx_start:idx_end], t_P[i+1], xf)
            plot_sf_leg(self._leg, units=AU, N=10, axes=ax, legend=False)

        return axes

    def high_fidelity(self, status = False):
        """
        prob.high_fidelity(status = True)

        - status: either True or False (True sets high fidelity on)

        Sets the trajectory high fidelity mode

        Example::

          prob.high_fidelity(True)
        """
        # We set the propagation fidelity
        self._leg.high_fidelity = status