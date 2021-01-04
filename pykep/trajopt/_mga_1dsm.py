from pykep.core import epoch, DAY2SEC, MU_SUN, lambert_problem, propagate_lagrangian, fb_prop, AU, epoch
from pykep.planet import jpl_lp
from pykep.trajopt._lambert import lambert_problem_multirev
from math import pi, cos, sin, acos, log, sqrt
import numpy as np
from typing import Any, List, Tuple

# Avoiding scipy dependency
def norm(x):
    return sqrt(sum([it * it for it in x]))

class mga_1dsm:
    r"""
    This class transcribes a Multiple Gravity Assist trajectory with one deep space maneuver per leg (MGA-1DSM) into an optimisation problem.
    It may be used as a User Defined Problem (UDP) for the pygmo (http://esa.github.io/pygmo/) optimisation suite.

    - Izzo, Dario. "Global optimization and space pruning for spacecraft trajectory design." Spacecraft Trajectory Optimization 1 (2010): 178-200.

    The decision vector (chromosome) is::

      direct encoding: [t0] + [u, v, Vinf, eta1, T1] + [beta, rp/rV, eta2, T2] + ... 
      alpha encoding:  [t0] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ... + [T]
      eta encoding:    [t0] + [u, v, Vinf, eta1, n1] + [beta, rp/rV, eta2, n2] + ...

      where t0 is a mjd2000, Vinf is in km/s, T in days, beta in radians and the rest non dimensional.

    .. note::

       The time of flights of a MGA-1DSM trajectory (and in general) can be encoded in different ways.
       When they are directly present in the decision vector, we have a *direct* encoding. This is the most 'evolvable' encoding
       but also the one that requires the most problem knowledge (e.g. to define the bounds on each leg) and is not 
       very flexible in dealing with constraints on the total time of flight. The *alpha* and *eta* encodings, instead, allow
       to only specify bounds on the time of flight of the entire trajectory, and not on the single legs: a property that is attractive
       for multi-objective optimization, for example.

       In the *alpha* encoding each leg time-of-flight is decoded as follows, T_i = T log(alpha_i) / \sum_n(log(alpha_n)).
       In the *eta* encoding  each leg time-of-flight is decoded as follows, T_i = (tof_max - \sum_0^(i-1)(T_j)) * eta_i

       The chromosome dimension for the direct and eta encoding is the same, while the alpha encoding requires one more gene.

    .. note::

       The resulting problem is box-bounded (unconstrained).
    """

    def __init__(self,
                 seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')],
                 t0 = [epoch(0), epoch(1000)],
                 tof = [[10, 300], [10, 300]],
                 vinf = [0.5, 2.5],
                 add_vinf_dep = False,
                 add_vinf_arr = True,
                 tof_encoding = 'direct',
                 multi_objective = False,
                 orbit_insertion = False,
                 e_target = None,
                 rp_target = None,
                 eta_lb = 0.1,
                 eta_ub = 0.9,
                 rp_ub = 30,
                 max_revs = 0
                 ):
        """
        pykep.trajopt.mga_1dsm(seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = [0.5, 2.5], multi_objective = False, add_vinf_dep = False, add_vinf_arr=True)

        - seq (``list`` of ``pykep.planet``): the encounter sequence (including the starting launch)
        - t0 (``list`` of ``pykep.epoch`` or ``floats``): the launch window (in mjd2000 if floats)
        - tof (``list`` or ``float``): bounds on the time of flight (days). If *tof_encoding* is 'direct', this contains a list
            of 2D lists defining the upper and lower bounds on each leg. If *tof_encoding* is 'alpha',
            this contains a list of two floats containing the lower and upper bounds on the time-of-flight. If *tof_encoding*
            is 'eta' tof is a float defining the upper bound on the time-of-flight
        - vinf (``list``): the minimum and maximum allowed initial hyperbolic velocity (at launch), in km/sec
        - add_vinf_dep (``bool``): when True the computed Dv includes the initial hyperbolic velocity (at launch)
        - add_vinf_arr (``bool``): when True the computed Dv includes the final hyperbolic velocity (at the last planet)
        - tof_encoding (``str``): one of 'direct', 'alpha' or 'eta'. Selects the encoding for the time of flights
        - multi_objective (``bool``): when True constructs a multiobjective problem (dv, T)
        - orbit_insertion (``bool``): when True the arrival dv is computed as that required to acquire a target orbit defined by e_target and rp_target
        - e_target (``float``): if orbit_insertion is True this defines the target orbit eccentricity around the final planet
        - rp_target (``float``): if orbit_insertion is True this defines the target orbit pericenter around the final planet (in m)
        - max_revs (``int``): maximal number of revolutions for lambert transfer
        """

        # Sanity checks
        # 1 - Planets need to have the same mu_central_body
        if ([r.mu_central_body for r in seq].count(seq[0].mu_central_body) != len(seq)):
            raise ValueError(
                'All planets in the sequence need to have identical mu_central_body')
        # 2 - tof encoding needs to be one of 'alpha', 'eta', 'direct'
        if tof_encoding not in ['alpha', 'eta', 'direct']:
            raise TypeError(
                'tof encoding must be one of \'alpha\', \'eta\', \'direct\'')
        # 3 - tof is expected to have different content depending on the tof_encoding
        if tof_encoding == 'direct':
            if np.shape(np.array(tof)) != (len(seq) - 1, 2):
                raise TypeError(
                    'tof_encoding is ' + tof_encoding + ' and tof must be a list of two dimensional lists and with length equal to the number of legs')
        if tof_encoding == 'alpha':
            if np.shape(np.array(tof)) != (2,):
                raise TypeError(
                    'tof_encoding is ' + tof_encoding + ' and tof must be a list of two floats')
        if tof_encoding == 'eta':
            if np.shape(np.array(tof)) != ():
                raise TypeError(
                    'tof_encoding is ' + tof_encoding + ' and tof must be a float')
        # 4 - Check launch window t0. If defined in terms of floats transform into epochs
        if len(t0) != 2:
            raise TypeError(
                    't0 is ' + t0 + ' while should be a list of two floats or epochs')
        if type(t0[0]) is not epoch:
            t0[0] = epoch(t0[0])
        if type(t0[1]) is not epoch:
            t0[1] = epoch(t0[1])
        # 5 - Check that if orbit insertion is selected e_target and r_p are
        # defined
        if orbit_insertion:
            if rp_target is None:
                raise ValueError(
                    'The rp_target needs to be specified when orbit insertion is selected')
            if e_target is None:
                raise ValueError(
                    'The e_target needs to be specified when orbit insertion is selected')
            if add_vinf_arr is False:
                raise ValueError(
                    'When orbit insertion is selected, the add_vinf_arr must be True')

        self._seq = seq
        self._t0 = t0
        self._tof = tof
        self._vinf = vinf
        self._add_vinf_dep = add_vinf_dep
        self._add_vinf_arr = add_vinf_arr
        self._tof_encoding = tof_encoding
        self._multi_objective = multi_objective
        self._orbit_insertion = orbit_insertion
        self._e_target = e_target
        self._rp_target = rp_target
        self._eta_lb = eta_lb
        self._eta_ub = eta_ub
        self._rp_ub = rp_ub

        self.n_legs = len(seq) - 1
        self.common_mu = seq[0].mu_central_body
        self.max_revs = max_revs

    def get_nobj(self):
        return self._multi_objective + 1

    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        vinf = self._vinf
        seq = self._seq
        # Base for all possiblities (eta encoding)
        lb = [t0[0].mjd2000] + [0.0, 0.0, vinf[0] * 1000, self._eta_lb,
                                1e-3] + [-2 * pi, np.nan, self._eta_lb, 1e-3] * (self.n_legs - 1)
        ub = [t0[1].mjd2000] + [1.0, 1.0, vinf[1] * 1000, self._eta_ub,
                                1.0 - 1e-3] + [2 * pi, self._rp_ub, self._eta_ub, 1.0 - 1e-3] * (self.n_legs - 1)
        # Distinguishing among cases (only direct and alpha)
        if self._tof_encoding == 'alpha':
            lb = lb + [tof[0]]
            ub = ub + [tof[1]]
        elif self._tof_encoding == 'direct':
            for i in range(self.n_legs):
                lb[5 + 4 * i] = tof[i][0]
                ub[5 + 4 * i] = tof[i][1]

        # Setting the minimum rp/rP using the planet safe radius
        for i, pl in enumerate(seq[1:-1]):
            lb[7 + 4 * i] = pl.safe_radius / pl.radius
        return (lb, ub)

    def _decode_times_and_vinf(self, x):
        # 1 - we decode the times of flight
        if self._tof_encoding == 'alpha':
            # decision vector is  [t0] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ... + [T]
            T = list([0] * (self.n_legs))
            for i in range(len(T)):
                T[i] = -log(x[5 + 4 * i])
            alpha_sum = sum(T)
            retval_T = [x[-1] * time / alpha_sum for time in T]
        elif self._tof_encoding == 'direct':
            # decision vector is  [t0] + [u, v, Vinf, eta1, T1] + [beta, rp/rV, eta2, T2] + ...
            retval_T = x[5::4]
        elif self._tof_encoding == 'eta':
            # decision vector is  [t0] + [u, v, Vinf, eta1, n1] + [beta, rp/rV, eta2, n2] + ...
            dt = self._tof
            T = [0] * self.n_legs
            for i in range(self.n_legs):
                T[i] = (dt - sum(T[:i])) * x[5 + 4*i]
            retval_T = T

        # 2 - we decode the hyperbolic velocity at departure
        theta = 2 * pi * x[1]
        phi = acos(2 * x[2] - 1) - pi / 2

        Vinfx = x[3] * cos(phi) * cos(theta)
        Vinfy = x[3] * cos(phi) * sin(theta)
        Vinfz = x[3] * sin(phi)

        return (retval_T, Vinfx, Vinfy, Vinfz)

    def _compute_dvs(self, x: List[float]) -> Tuple[
        List[float], # DVs
        List[Any], # Lambert legs
        List[float], # T
        List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]], # ballistic legs
        List[float], # episodes of ballistic legs
    ]:
        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T and the cartesian components of vinf
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.n_legs + 1))
        r_P = list([None] * (self.n_legs + 1))
        v_P = list([None] * (self.n_legs + 1))
        DV = list([0.0] * (self.n_legs + 1))
        for i in range(len(self._seq)):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])
        ballistic_legs = []
        ballistic_ep = []
        lamberts = []

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        ballistic_legs.append((r_P[0], v0))
        ballistic_ep.append(t_P[0].mjd2000)
        r, v = propagate_lagrangian(
            r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * DAY2SEC
        l = lambert_problem_multirev(v, lambert_problem(
                    r, r_P[1], dt, self.common_mu, cw=False, max_revs=self.max_revs))
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
        lamberts.append(l)

        ballistic_legs.append((r, v_beg_l))
        ballistic_ep.append(t_P[0].mjd2000 + x[4] * T[0])

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
            ballistic_legs.append((r_P[i], v_out))
            ballistic_ep.append(t_P[i].mjd2000)
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem_multirev(v, lambert_problem(r, r_P[i + 1], dt,
                                  self.common_mu, cw=False, max_revs=self.max_revs))
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            lamberts.append(l)
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])

            ballistic_legs.append((r, v_beg_l))
            ballistic_ep.append(t_P[i].mjd2000 + x[8 + (i - 1) * 4] * T[i])

        # Last Delta-v
        if self._add_vinf_arr:
            DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
            if self._orbit_insertion:
                # In this case we compute the insertion DV as a single pericenter
                # burn
                DVper = np.sqrt(DV[-1] * DV[-1] + 2 *
                                self._seq[-1].mu_self / self._rp_target)
                DVper2 = np.sqrt(2 * self._seq[-1].mu_self / self._rp_target -
                                self._seq[-1].mu_self / self._rp_target * (1. - self._e_target))
                DV[-1] = np.abs(DVper - DVper2)

        if self._add_vinf_dep:
            DV[0] += x[3]

        return (DV, lamberts, T, ballistic_legs, ballistic_ep)

    # Objective function
    def fitness(self, x):
        DV, _, T, _, _ = self._compute_dvs(x)
        if not self._multi_objective:
            return [sum(DV),]
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
        t_P = list([None] * (self.n_legs + 1))
        r_P = list([None] * (self.n_legs + 1))
        v_P = list([None] * (self.n_legs + 1))
        DV = list([0.0] * (self.n_legs + 1))
        for i in range(len(self._seq)):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        print("First Leg: " + self._seq[0].name + " to " + self._seq[1].name)
        print("Departure: " + str(t_P[0]) +
              " (" + str(t_P[0].mjd2000) + " mjd2000) ")
        print("Duration: " + str(T[0]) + "days")
        print("VINF: " + str(x[3] / 1000) + " km/sec")

        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(
            r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

        print("DSM after " + str(x[4] * T[0]) + " days")

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * DAY2SEC
        l = lambert_problem_multirev(v, lambert_problem(
            r, r_P[1], dt, self.common_mu, cw=False, max_revs=self.max_revs))
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])
        print("DSM magnitude: " + str(DV[0]) + "m/s")

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            print("\nleg no. " + str(i + 1) + ": " +
                  self._seq[i].name + " to " + self._seq[i + 1].name)
            print("Duration: " + str(T[i]) + "days")
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
            print(
                "Fly-by epoch: " + str(t_P[i]) + " (" + str(t_P[i].mjd2000) + " mjd2000) ")
            print(
                "Fly-by radius: " + str(x[7 + (i - 1) * 4]) + " planetary radii")
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            print("DSM after " + str(x[8 + (i - 1) * 4] * T[i]) + " days")
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC
            l = lambert_problem_multirev(v, lambert_problem(r, r_P[i + 1], dt,
                                self.common_mu, cw=False, max_revs=self.max_revs))
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
            print("DSM magnitude: " + str(DV[i]) + "m/s")

        # Last Delta-v
        print("\nArrival at " + self._seq[-1].name)
        DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
        print(
            "Arrival epoch: " + str(t_P[-1]) + " (" + str(t_P[-1].mjd2000) + " mjd2000) ")
        print("Arrival Vinf: " + str(DV[-1]) + "m/s")
        if self._orbit_insertion:
            # In this case we compute the insertion DV as a single pericenter
            # burn
            DVper = np.sqrt(DV[-1] * DV[-1] + 2 *
                            self._seq[-1].mu_self / self._rp_target)
            DVper2 = np.sqrt(2 * self._seq[-1].mu_self / self._rp_target -
                            self._seq[-1].mu_self / self._rp_target * (1. - self._e_target))
            DVinsertion = np.abs(DVper - DVper2)
            print("Insertion DV: " + str(DVinsertion) + "m/s")

        print("Total mission time: " + str(sum(T) / 365.25) + " years (" + str(sum(T)) + " days)")

    # Plot of the trajectory
    def plot(self, x, ax = None):
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
        t_P = list([None] * (self.n_legs + 1))
        r_P = list([None] * (self.n_legs + 1))
        v_P = list([None] * (self.n_legs + 1))
        DV = list([None] * (self.n_legs + 1))

        for i, planet in enumerate(self._seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = planet.eph(t_P[i])
            plot_planet(planet, t0=t_P[i], color=(
                0.8, 0.6, 0.8), legend=True, units=AU, axes=axis, N=150)

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = propagate_lagrangian(
            r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu)

        plot_kepler(r_P[0], v0, x[4] * T[0] * DAY2SEC, self.common_mu,
                    N=100, color='b', units=AU, axes=axis)

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * DAY2SEC
        
        l = lambert_problem_multirev(v, lambert_problem(
            r, r_P[1], dt, self.common_mu, cw=False, max_revs=self.max_revs))
        
        plot_lambert(l, sol=0, color='r', units=AU, axes=axis)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            # Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[
                            7 + (i - 1) * 4] * self._seq[i].radius, x[6 + (i - 1) * 4], self._seq[i].mu_self)
            # s/c propagation before the DSM
            r, v = propagate_lagrangian(
                r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC, self.common_mu)
            plot_kepler(r_P[i], v_out, x[8 + (i - 1) * 4] * T[i] * DAY2SEC,
                        self.common_mu, N=100, color='b', units=AU, axes=axis)
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * DAY2SEC

            l = lambert_problem_multirev(v, lambert_problem(r, r_P[i + 1], dt,
                self.common_mu, cw=False, max_revs=self.max_revs))

            plot_lambert(l, sol=0, color='r', legend=False,
                         units=AU, N=1000, axes=axis)

            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        plt.show()
        return axis

    def get_extra_info(self):
        return ("\n\t Sequence: " + [pl.name for pl in self._seq].__repr__() +
                "\n\t Add launcher vinf to the objective?: " + self._add_vinf_dep.__repr__() +
                "\n\t Add final vinf to the objective?: " + self._add_vinf_arr.__repr__())
