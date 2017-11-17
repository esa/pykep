from pykep.core import epoch, DAY2SEC, MU_SUN, lambert_problem, propagate_lagrangian, fb_prop, AU
from pykep.planet import jpl_lp
from math import pi, cos, sin, acos, log
from scipy.linalg import norm


class mga_1dsm:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) problem representing a Multiple Gravity Assist
    trajectory allowing one only impulsive Deep Space Manouvre between each leg.

    - Izzo, Dario. "Global optimization and space pruning for spacecraft trajectory design." Spacecraft Trajectory Optimization 1 (2010): 178-200.

    The decision vector is::

      [t0, T] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ...

    ... in the units: [mjd2000, days] + [nd,nd,km/s,nd,years,nd] + [rad,nd,nd,nd] + ....

    Each leg time-of-flight can be decoded as follows, T_n = T log(alpha_n) / \sum_i(log(alpha_i)).

    .. note::

       The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
    """

    def __init__(self,
                 seq=[jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')],
                 t0=[epoch(0), epoch(1000)],
                 tof=[1.0, 5.0],
                 vinf=[0.5, 2.5],
                 add_vinf_dep=False,
                 add_vinf_arr=True,
                 multi_objective=False):
        """
        pykep.trajopt.mga_1dsm(seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr=True)

        - seq: list of pykep planets defining the encounter sequence (including the starting launch)
        - t0: list of two epochs defining the launch window
        - tof: list of two floats defining the minimum and maximum allowed mission lenght (days)
        - vinf: list of two floats defining the minimum and maximum allowed initial hyperbolic velocity (at launch), in km/sec
        - multi_objective: when True constructs a multiobjective problem (dv, T)
        - add_vinf_dep: when True the computed Dv includes the initial hyperbolic velocity (at launch)
        - add_vinf_arr: when True the computed Dv includes the final hyperbolic velocity (at the last planet)
        """

        # Sanity checks ...... all planets need to have the same
        # mu_central_body
        if ([r.mu_central_body for r in seq].count(seq[0].mu_central_body) != len(seq)):
            raise ValueError(
                'All planets in the sequence need to have exactly the same mu_central_body')
        self.__add_vinf_dep = add_vinf_dep
        self.__add_vinf_arr = add_vinf_arr
        self.__n_legs = len(seq) - 1
        self._t0 = t0
        self._tof = tof
        self._vinf = vinf
        self._obj_dim = multi_objective + 1

        # We then define all planets in the sequence  and the common central
        # body gravity as data members
        self.seq = seq
        self.common_mu = seq[0].mu_central_body

    def get_nobj(self):
        return self._obj_dim

    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        vinf = self._vinf
        seq = self.seq
        lb = [t0[0].mjd2000, tof[0]] + [0.0, 0.0, vinf[0] * 1000, 1e-5,
                                        1e-5] + [-2 * pi, 1.1, 1e-5, 1e-5] * (self.__n_legs - 1)
        ub = [t0[1].mjd2000, tof[1]] + [1.0, 1.0, vinf[1] * 1000, 1.0 - 1e-5,
                                        1.0 - 1e-5] + [2 * pi, 30.0, 1.0 - 1e-5, 1.0 - 1e-5] * (self.__n_legs - 1)
        # Accounting that each planet has a different safe radius......
        for i, pl in enumerate(seq[1:-1]):
            lb[8 + 4 * i] = pl.safe_radius / pl.radius
        return (lb, ub)

    def _decode_times_and_vinf(self, x):
        T = list([0] * (self.__n_legs))
        for i in range(len(T)):
            T[i] = -log(x[6 + 4 * i])
        alpha_sum = sum(T)

        theta = 2 * pi * x[2]
        phi = acos(2 * x[3] - 1) - pi / 2

        Vinfx = x[4] * cos(phi) * cos(theta)
        Vinfy = x[4] * cos(phi) * sin(theta)
        Vinfz = x[4] * sin(phi)

        return ([x[1] * time / alpha_sum for time in T], Vinfx, Vinfy, Vinfz)

    # Objective function
    def fitness(self, x):
        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T and the cartesian components of vinf
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 1))
        r_P = list([None] * (self.__n_legs + 1))
        v_P = list([None] * (self.__n_legs + 1))
        DV = list([0.0] * (self.__n_legs + 1))
        for i, planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(
            r_P[0], v0, x[5] * T[0] * DAY2SEC, self.common_mu)

        # Lambert arc to reach seq[1]
        dt = (1 - x[5]) * T[0] * DAY2SEC
        l = lambert_problem(r, r_P[1], dt, self.common_mu, False, False)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.__n_legs):
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            8 + (i - 1) * 4] * self.seq[i].radius, x[7 + (i - 1) * 4], self.seq[i].mu_self)
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[9 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[9 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem(r, r_P[i + 1], dt,
                                self.common_mu, False, False)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])

        # Last Delta-v
        if self.__add_vinf_arr:
            DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])

        if self.__add_vinf_dep:
            DV[0] += x[3]

        if self._obj_dim == 1:
            return (sum(DV),)
        else:
            return (sum(DV), sum(T))

    def pretty(self, x):
        """
        prob.plot(x)

        - x: encoded trajectory

        Prints human readable information on the trajectory represented by the decision vector x

        Example::

          print(prob.pretty(x))
        """
        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T and the cartesian components of vinf
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 1))
        r_P = list([None] * (self.__n_legs + 1))
        v_P = list([None] * (self.__n_legs + 1))
        DV = list([None] * (self.__n_legs + 1))

        for i, planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        print("First Leg: " + self.seq[0].name + " to " + self.seq[1].name)
        print("Departure: " + str(t_P[0]) +
              " (" + str(t_P[0].mjd2000) + " mjd2000) ")
        print("Duration: " + str(T[0]) + "days")
        print("VINF: " + str(x[4] / 1000) + " km/sec")

        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(
            r_P[0], v0, x[5] * T[0] * DAY2SEC, self.common_mu)

        print("DSM after " + str(x[5] * T[0]) + " days")

        # Lambert arc to reach seq[1]
        dt = (1 - x[5]) * T[0] * DAY2SEC
        l = lambert_problem(r, r_P[1], dt, self.common_mu, False, False)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])
        print("DSM magnitude: " + str(DV[0]) + "m/s")

        # 4 - And we proceed with each successive leg
        for i in range(1, self.__n_legs):
            print("\nleg no. " + str(i + 1) + ": " +
                  self.seq[i].name + " to " + self.seq[i + 1].name)
            print("Duration: " + str(T[i]) + "days")
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            8 + (i - 1) * 4] * self.seq[i].radius, x[7 + (i - 1) * 4], self.seq[i].mu_self)
            print(
                "Fly-by epoch: " + str(t_P[i]) + " (" + str(t_P[i].mjd2000) + " mjd2000) ")
            print(
                "Fly-by radius: " + str(x[8 + (i - 1) * 4]) + " planetary radii")
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[9 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            print("DSM after " + str(x[9 + (i - 1) * 4] * T[i]) + " days")
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[9 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem(r, r_P[i + 1], dt,
                                self.common_mu, False, False)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
            print("DSM magnitude: " + str(DV[i]) + "m/s")

        # Last Delta-v
        print("\nArrival at " + self.seq[-1].name)
        DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
        print(
            "Arrival epoch: " + str(t_P[-1]) + " (" + str(t_P[-1].mjd2000) + " mjd2000) ")
        print("Arrival Vinf: " + str(DV[-1]) + "m/s")
        print("Total mission time: " + str(sum(T) / 365.25) + " years")

    # Plot of the trajectory
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
        from pykep.orbit_plots import plot_planet, plot_lambert, plot_kepler

        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.gca(projection='3d')
        else:
            axis = ax

        axis.scatter(0, 0, 0, color='y')

        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T and the cartesian components of vinf
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 1))
        r_P = list([None] * (self.__n_legs + 1))
        v_P = list([None] * (self.__n_legs + 1))
        DV = list([None] * (self.__n_legs + 1))

        for i, planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = planet.eph(t_P[i])
            plot_planet(planet, t0=t_P[i], color=(
                0.8, 0.6, 0.8), legend=True, units=AU, ax=axis)

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(
            r_P[0], v0, x[5] * T[0] * DAY2SEC, self.common_mu)

        plot_kepler(r_P[0], v0, x[5] * T[0] * DAY2SEC, self.common_mu,
                    N=100, color='b', legend=False, units=AU, ax=axis)

        # Lambert arc to reach seq[1]
        dt = (1 - x[5]) * T[0] * DAY2SEC
        l = lambert_problem(r, r_P[1], dt, self.common_mu, False, False)
        plot_lambert(l, sol=0, color='r', legend=False, units=AU, ax=axis)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.__n_legs):
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            8 + (i - 1) * 4] * self.seq[i].radius, x[7 + (i - 1) * 4], self.seq[i].mu_self)
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[9 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            plot_kepler(r_P[i], v_out, x[9 + (i - 1) * 4] * T[i] * DAY2SEC,
                        self.common_mu, N=100, color='b', legend=False, units=AU, ax=axis)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[9 + (i - 1) * 4]) * T[i] * DAY2SEC

            l = lambert_problem(r, r_P[i + 1], dt,
                                self.common_mu, False, False)
            plot_lambert(l, sol=0, color='r', legend=False,
                         units=AU, N=1000, ax=axis)

            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        plt.show()
        return axis

    def get_extra_info(self):
        return ("\n\t Sequence: " + [pl.name for pl in self.seq].__repr__() +
                "\n\t Add launcher vinf to the objective?: " + self.__add_vinf_dep.__repr__() +
                "\n\t Add final vinf to the objective?: " + self.__add_vinf_arr.__repr__())
