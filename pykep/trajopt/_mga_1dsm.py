import pykep as _pk
import numpy as _np
from math import pi, cos, sin, acos, log, sqrt
from copy import deepcopy
from typing import Any, Dict, List, Tuple
from bisect import bisect_left as _bisect_left


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
    def __init__(
        self,
        seq=[
            _pk.planet(_pk.udpla.jpl_lp("earth")),
            _pk.planet(_pk.udpla.jpl_lp("venus")),
            _pk.planet(_pk.udpla.jpl_lp("earth")),
        ],
        t0=[0, 1000],
        tof=[[30, 200], [200, 300]],
        vinf=[0.5, 2.5],
        add_vinf_dep=False,
        add_vinf_arr=True,
        tof_encoding="direct",
        multi_objective=False,
        orbit_insertion=False,
        e_target=None,
        rp_target=None,
        eta_bounds=[0.1, 0.9],
        rp_ub=30,
    ):
        """
        pykep.trajopt.mga_1dsm(seq, t0, tof, vinf_bounds = [0.5, 2.5], add_vinf_dep=False, add_vinf_arr=True, tof_encoding="direct", multi_objective=False, orbit_insertion=False, e_target=None, rp_target=None, eta_bounds=[0.01, 0.99], rp_ub=30)

            *seq* (:class:`list` [:class:`~pykep.planet`]): sequence of planetary encounters including the departure body.

            *t0* (:class:`list` [:class:`float` or :class:`~pykep.epoch`]): lower and upper bounds for the launch epoch. When floats are used MJD2000 is assumed.

            *tof* (:class:`list` or :class:`float`): defines the bounds on the time of flight. If *tof_encoding* is 'direct', this contains a list
            of 2D lists defining the upper and lower bounds on each leg. If *tof_encoding* is 'alpha',
            this contains a 2D list with the lower and upper bounds on the total time-of-flight. If *tof_encoding*
            is 'eta', tof is a float defining an upper bound for the time-of-flight.

            *vinf_bounds* (:class:`list` [:class:`float`]): the minimum and maximum allowed initial hyperbolic velocity (at launch), in km/sec.

            *add_vinf_dep* (:class:`bool`): when True the computed Dv includes the initial hyperbolic velocity (at launch).

            *add_vinf_arr* (:class:`bool`): when True the computed Dv includes the final hyperbolic velocity (at the last planet).

            *tof_encoding* (:class:`str`): one of 'direct', 'alpha' or 'eta'. Selects the encoding for the time of flights.

            *multi_objective* (:class:`bool`): when True constructs a multiobjective problem (dv, T).

            *orbit_insertion* (:class:`bool`): when True the arrival dv is computed as that required to acquire a target orbit defined by e_target and rp_target.

            *e_target* (:class:`float`): if orbit_insertion is True this defines the target orbit eccentricity around the final planet.

            *rp_target* (:class:`float`): if orbit_insertion is True this defines the target orbit pericenter around the final planet (in m).

            *eta_bounds* (:class:`list` [:class:`float`]): lower and upper bounds for the eta variable (defining the position of the DSM).

            *rp_ub* (:class:`float`): upper bound for the pericenter radius of the flybys (in planetary radii).
        """

        # Sanity checks
        # 1 - Planets need to have the same mu_central_body
        if [pla.mu_central_body for pla in seq].count(seq[0].mu_central_body) != len(
            seq
        ):
            raise ValueError(
                "All planets in the sequence need to have identical mu_central_body"
            )

        # 2 - tof encoding needs to be one of 'alpha', 'eta', 'direct'
        if tof_encoding not in ["alpha", "eta", "direct"]:
            raise TypeError("tof encoding must be one of 'alpha', 'eta', 'direct'")

        # 3 - tof is expected to have different content depending on the tof_encoding
        if tof_encoding == "direct":
            if _np.shape(_np.array(tof)) != (len(seq) - 1, 2):
                raise TypeError(
                    "tof_encoding is "
                    + tof_encoding
                    + " and tof must be a list of two dimensional lists and with length equal to the number of legs"
                )
        if tof_encoding == "alpha":
            if _np.shape(_np.array(tof)) != (2,):
                raise TypeError(
                    "tof_encoding is "
                    + tof_encoding
                    + " and tof must be a list of two floats"
                )
        if tof_encoding == "eta":
            if _np.shape(_np.array(tof)) != ():
                raise TypeError(
                    "tof_encoding is " + tof_encoding + " and tof must be a float"
                )

        # 4 - Check launch window t0. If defined in terms of floats transform into epochs
        if len(t0) != 2:
            raise TypeError(
                "t0 is " + t0 + " while should be a list of two floats or epochs"
            )
        # 4b - We try to build epochs out of the t0 list (mjd2000 by default)
        for i in range(len(t0)):
            if type(t0[i]) != type(_pk.epoch(0.0)):
                t0[i] = _pk.epoch(t0[i], _pk.epoch.julian_type.MJD2000)

        # 5 - Check that if orbit insertion is selected e_target and r_p are
        # defined
        if orbit_insertion:
            if rp_target is None:
                raise ValueError(
                    "The rp_target needs to be specified when orbit insertion is selected"
                )
            if e_target is None:
                raise ValueError(
                    "The e_target needs to be specified when orbit insertion is selected"
                )
            if add_vinf_arr is False:
                raise ValueError(
                    "When orbit insertion is selected, the add_vinf_arr must be True"
                )

        # 6 - Check that the eta bounds are a list of two floats
        if len(eta_bounds) != 2:
            raise ValueError("The eta_bounds must be a list of two floats")
        if type(eta_bounds[0]) != type(0.0) or type(eta_bounds[1]) != type(0.0):
            raise ValueError("The eta_bounds must be a list of two floats")

        # Public data members
        self.n_legs = len(seq) - 1
        self.common_mu = seq[0].mu_central_body

        # Private data members
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
        self._eta_lb = eta_bounds[0]
        self._eta_ub = eta_bounds[1]
        self._rp_ub = rp_ub

    def get_nobj(self):
        return self._multi_objective + 1

    def get_bounds(self):
        t0 = self._t0
        tof = self._tof
        vinf = self._vinf
        seq = self._seq
        # Base for all possiblities (eta encoding)
        lb = (
            [t0[0].mjd2000]
            + [0.0, 0.0, vinf[0] * 1000, self._eta_lb, 1e-3]
            + [-2 * pi, _np.nan, self._eta_lb, 1e-3] * (self.n_legs - 1)
        )
        ub = (
            [t0[1].mjd2000]
            + [1.0, 1.0, vinf[1] * 1000, self._eta_ub, 1.0 - 1e-3]
            + [2 * pi, self._rp_ub, self._eta_ub, 1.0 - 1e-3] * (self.n_legs - 1)
        )
        # Distinguishing among cases (only direct and alpha)
        if self._tof_encoding == "alpha":
            lb = lb + [tof[0]]
            ub = ub + [tof[1]]
        elif self._tof_encoding == "direct":
            for i in range(self.n_legs):
                lb[5 + 4 * i] = tof[i][0]
                ub[5 + 4 * i] = tof[i][1]

        # Setting the minimum rp/rP using the planet safe radius
        for i, pl in enumerate(seq[1:-1]):
            lb[7 + 4 * i] = pl.safe_radius / pl.radius
        return (lb, ub)

    def _decode_times_and_vinf(self, x):
        # 1 - we decode the times of flight
        if self._tof_encoding == "alpha":
            # decision vector is  [t0] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ... + [T]
            retval_T = _pk.alpha2direct(x[5::4], x[-1])
        elif self._tof_encoding == "direct":
            # decision vector is  [t0] + [u, v, Vinf, eta1, T1] + [beta, rp/rV, eta2, T2] + ...
            retval_T = x[5::4]
        elif self._tof_encoding == "eta":
            # decision vector is [t0] + [u, v, Vinf, eta1, n1] + [beta, rp/rV, eta2, n2] + ...
            retval_T = _pk.eta2direct(x[5::4], self._tof)

        # 2 - We decode the hyperbolic velocity at departure
        theta = 2 * pi * x[1]
        phi = acos(2 * x[2] - 1) - pi / 2

        Vinfx = x[3] * cos(phi) * cos(theta)
        Vinfy = x[3] * cos(phi) * sin(theta)
        Vinfz = x[3] * sin(phi)

        return (retval_T, Vinfx, Vinfy, Vinfz)

    def _decode_tofs(self, x):
        T, _, _, _ = self._decode_times_and_vinf(x)
        return T
    
    @staticmethod
    def alpha2direct(x):
        """alpha2direct(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA1DSM trajectory in the alpha encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA1DSM trajectory using the direct encoding
        """
        # decision vector is  [t0] + [u, v, Vinf, eta1, a1] + [beta, rp/rV, eta2, a2] + ... + [T]
        retval = deepcopy(x)
        retval[5::4] = _pk.alpha2direct(x[5::4], x[-1])
        retval = _np.delete(retval, -1)
        return retval

    @staticmethod
    def direct2alpha(x):
        """direct2alpha(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the alpha encoding
        """
        # decision vector is  [t0] + [u, v, Vinf, eta1, T1] + [beta, rp/rV, eta2, T2] + ...
        retval = deepcopy(x)
        retval[5::4], T = _pk.direct2alpha(x[5::4])
        retval = _np.append(retval,T)
        return retval

    @staticmethod
    def eta2direct(x, max_tof):
        """eta2direct(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the eta encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the direct encoding
        """
        # decision vector is [t0] + [u, v, Vinf, eta1, n1] + [beta, rp/rV, eta2, n2] + ...
        retval = deepcopy(x)
        retval[5::4] = _pk.eta2direct(x[5::4], max_tof)
        return retval

    @staticmethod
    def direct2eta(x, max_tof):
        """direct2eta(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the eta encoding
        """
        retval = deepcopy(x)
        retval[5::4] = _pk.direct2eta(x[5::4], max_tof)
        return retval

    def _compute_dvs(self, x: List[float]) -> Tuple[
        List[float],  # DVs
        List[Any],  # Lambert legs
        List[float],  # T
        List[Tuple[List[float], List[float]]],  # ballistic legs
        List[float],  # epochs of ballistic legs
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
            t_P[i] = _pk.epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])
        ballistic_legs: List[Tuple[List[float], List[float]]] = []
        ballistic_ep: List[float] = []
        lamberts = []

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        ballistic_legs.append((r_P[0], v0))
        ballistic_ep.append(t_P[0].mjd2000)
        r, v = _pk.propagate_lagrangian(
            [r_P[0], v0], x[4] * T[0] * _pk.DAY2SEC, self.common_mu
        )

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * _pk.DAY2SEC
        l = _pk.lambert_problem(r, r_P[1], dt, self.common_mu, cw=False, multi_revs=0)

        v_end_l = l.v1[0]
        v_beg_l = l.v0[0]
        lamberts.append(l)

        ballistic_legs.append((r, v_beg_l))
        ballistic_ep.append(t_P[0].mjd2000 + x[4] * T[0])

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            # Fly-by
            v_out = _pk.fb_vout(
                v_in=v_end_l,
                v_pla=v_P[i],
                rp=x[7 + (i - 1) * 4] * self._seq[i].radius,
                beta=x[6 + (i - 1) * 4],
                mu=self._seq[i].mu_self,
            )
            ballistic_legs.append((r_P[i], v_out))
            ballistic_ep.append(t_P[i].mjd2000)
            # s/c propagation before the DSM
            r, v = _pk.propagate_lagrangian(
                [r_P[i], v_out], x[8 + (i - 1) * 4] * T[i] * _pk.DAY2SEC, self.common_mu
            )
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * _pk.DAY2SEC
            l = _pk.lambert_problem(
                r, r_P[i + 1], dt, self.common_mu, cw=False, multi_revs=0
            )
            v_end_l = l.v1[0]
            v_beg_l = l.v0[0]
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
                DVper = _np.sqrt(
                    DV[-1] * DV[-1] + 2 * self._seq[-1].mu_self / self._rp_target
                )
                DVper2 = _np.sqrt(
                    2 * self._seq[-1].mu_self / self._rp_target
                    - self._seq[-1].mu_self / self._rp_target * (1.0 - self._e_target)
                )
                DV[-1] = _np.abs(DVper - DVper2)

        if self._add_vinf_dep:
            DV[0] += x[3]

        return (DV, lamberts, T, ballistic_legs, ballistic_ep)

    # Objective function
    def fitness(self, x):
        DV, _, T, _, _ = self._compute_dvs(x)
        if not self._multi_objective:
            return [
                sum(DV),
            ]
        else:
            return (sum(DV), sum(T))

    def pretty(self, x):
        """
        prob.pretty(x)

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
            t_P[i] = _pk.epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self._seq[i].eph(t_P[i])

        # 3 - We start with the first leg
        print("First Leg: " + self._seq[0].name + " to " + self._seq[1].name)
        print("Departure: " + str(t_P[0]) + " (" + str(t_P[0].mjd2000) + " mjd2000) ")
        print("Duration: " + str(T[0]) + "days")
        print("VINF: " + str(x[3] / 1000) + " km/sec")

        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = _pk.propagate_lagrangian(
            [r_P[0], v0], x[4] * T[0] * _pk.DAY2SEC, self.common_mu
        )

        print("DSM after " + str(x[4] * T[0]) + " days")

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * _pk.DAY2SEC
        l = _pk.lambert_problem(r, r_P[1], dt, self.common_mu, cw=False, multi_revs=0)
        v_end_l = l.v1[0]
        v_beg_l = l.v0[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])
        print("DSM magnitude: " + str(DV[0]) + "m/s")

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            print(
                "\nleg no. "
                + str(i + 1)
                + ": "
                + self._seq[i].name
                + " to "
                + self._seq[i + 1].name
            )
            print("Duration: " + str(T[i]) + "days")
            # Fly-by
            v_out = _pk.fb_vout(
                v_in=v_end_l,
                v_pla=v_P[i],
                rp=x[7 + (i - 1) * 4] * self._seq[i].radius,
                beta=x[6 + (i - 1) * 4],
                mu=self._seq[i].mu_self,
            )
            print(
                "Fly-by epoch: "
                + str(t_P[i])
                + " ("
                + str(t_P[i].mjd2000)
                + " mjd2000) "
            )
            print("Fly-by radius: " + str(x[7 + (i - 1) * 4]) + " planetary radii")
            # s/c propagation before the DSM
            r, v = _pk.propagate_lagrangian(
                [r_P[i], v_out], x[8 + (i - 1) * 4] * T[i] * _pk.DAY2SEC, self.common_mu
            )
            print("DSM after " + str(x[8 + (i - 1) * 4] * T[i]) + " days")
            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * _pk.DAY2SEC
            l = _pk.lambert_problem(
                r, r_P[i + 1], dt, self.common_mu, cw=False, multi_revs=0
            )
            v_end_l = l.v1[0]
            v_beg_l = l.v0[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
            print("DSM magnitude: " + str(DV[i]) + "m/s")

        # Last Delta-v
        print("\nArrival at " + self._seq[-1].name)
        DV[-1] = norm([a - b for a, b in zip(v_end_l, v_P[-1])])
        print(
            "Arrival epoch: "
            + str(t_P[-1])
            + " ("
            + str(t_P[-1].mjd2000)
            + " mjd2000) "
        )
        print("Arrival Vinf: " + str(DV[-1]) + "m/s")
        if self._orbit_insertion:
            # In this case we compute the insertion DV as a single pericenter
            # burn
            DVper = _np.sqrt(
                DV[-1] * DV[-1] + 2 * self._seq[-1].mu_self / self._rp_target
            )
            DVper2 = _np.sqrt(
                2 * self._seq[-1].mu_self / self._rp_target
                - self._seq[-1].mu_self / self._rp_target * (1.0 - self._e_target)
            )
            DVinsertion = _np.abs(DVper - DVper2)
            print("Insertion DV: " + str(DVinsertion) + "m/s")

        print(
            "Total mission time: "
            + str(sum(T) / 365.25)
            + " years ("
            + str(sum(T))
            + " days)"
        )

    # Plot of the trajectory
    def plot(
        self,
        x,
        ax=None,
        units=_pk.AU,
        N=60,
        c_orbit="dimgray",
        c_lambert="indianred",
        c_ballistic="royalblue",
        leg_ids = [],
        figsize=(5, 5),
        **kwargs
    ):
        """
        Plots the trajectory encoded into *x* in 3D axes.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.

            *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): the 3D axis. If None a new one is created.

            *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.

            *N* (:class:`int`, optional): number of points to be used in the plots. Defaults to 60.

            *c_orbit* (:class:`str`, optional): color of the planetary orbits. Defaults to 'dimgray'.

            *c* (:class:`str`, optional): color of the trajectory. Defaults to 'indianred'.

            *figsize* (:class:`tuple`, optional): size of the figure. Defaults to (5, 5).

            *\\*\\*kwargs*: Additional keyword arguments to pass to the trajectory plotting functions.
        """
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from pykep.plot import (
            add_ballistic_arc,
            add_lambert,
            add_planet,
            add_planet_orbit,
            make_3Daxis,
            add_sun,
        )

        # If no axis is provided we create a new one
        if ax is None:
            ax = make_3Daxis(figsize=figsize)
        else:
            ax = ax
            
        # Plot of leg unless specified
        if len(leg_ids) == 0:
            leg_ids = list(range(self.n_legs))

        # We add the sun to the plot
        add_sun(ax)

        # 1 -  we 'decode' the chromosome recording the various times of flight
        # (days) in the list T and the cartesian components of vinf
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.n_legs + 1))
        r_P = list([None] * (self.n_legs + 1))
        v_P = list([None] * (self.n_legs + 1))
        DV = list([None] * (self.n_legs + 1))

        for i, item in enumerate(self._seq):
            t_P[i] = _pk.epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = item.eph(t_P[i])
            if i in leg_ids:
                add_planet(ax, pla=item, when=t_P[i], units=units, c=c_orbit)
            add_planet_orbit(ax, pla=item, units=units, N=N, c=c_orbit)

        # 3 - We start with the first leg
        v0 = [a + b for a, b in zip(v_P[0], [Vinfx, Vinfy, Vinfz])]
        r, v = _pk.propagate_lagrangian(
            [r_P[0], v0], x[4] * T[0] * _pk.DAY2SEC, self.common_mu
        )
        if 0 in leg_ids:
            add_ballistic_arc(
                ax,
                [r_P[0], v0],
                x[4] * T[0] * _pk.DAY2SEC,
                self.common_mu,
                N=N,
                units=units,
                c=c_ballistic,
                **kwargs
            )

        # Lambert arc to reach seq[1]
        dt = (1 - x[4]) * T[0] * _pk.DAY2SEC

        l = _pk.lambert_problem(r, r_P[1], dt, self.common_mu, cw=False, multi_revs=0)
        if 0 in leg_ids:
            add_lambert(ax, lp=l, N=N, sol=0, units=units, c=c_lambert, **kwargs)
        v_end_l = l.v1[0]
        v_beg_l = l.v0[0]

        # First DSM occuring at time nu1*T1
        DV[0] = norm([a - b for a, b in zip(v_beg_l, v)])

        # 4 - And we proceed with each successive leg
        for i in range(1, self.n_legs):
            # Fly-by
            v_out = _pk.fb_vout(
                v_end_l,
                v_P[i],
                x[7 + (i - 1) * 4] * self._seq[i].radius,
                x[6 + (i - 1) * 4],
                self._seq[i].mu_self,
            )
            # s/c propagation before the DSM
            r, v = _pk.propagate_lagrangian(
                [r_P[i], v_out], x[8 + (i - 1) * 4] * T[i] * _pk.DAY2SEC, self.common_mu
            )
            if i in leg_ids:
                add_ballistic_arc(
                    ax,
                    [r_P[i], v_out],
                    x[8 + (i - 1) * 4] * T[i] * _pk.DAY2SEC,
                    self.common_mu,
                    N=N,
                    units=units,
                    c=c_ballistic,
                    **kwargs
                )

            # Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1 - x[8 + (i - 1) * 4]) * T[i] * _pk.DAY2SEC

            l = _pk.lambert_problem(
                r, r_P[i + 1], dt, self.common_mu, cw=False, multi_revs=0
            )
            if i in leg_ids:
                add_lambert(ax, lp=l, sol=0, units=units, N=N, c=c_lambert, **kwargs)

            v_end_l = l.v1[0]
            v_beg_l = l.v0[0]
            # DSM occuring at time nu2*T2
            DV[i] = norm([a - b for a, b in zip(v_beg_l, v)])
        return ax

    def get_extra_info(self):
        return (
            "\t Sequence: "
            + [pl.name for pl in self._seq].__repr__()
            + "\n\t Add launcher vinf to the objective?: "
            + self._add_vinf_dep.__repr__()
            + "\n\t Add final vinf to the objective?: "
            + self._add_vinf_arr.__repr__()
        )

    def get_name(self):
        return "MGA_1DSM Trajectory"

    def __repr__(self):
        return self.get_name()

    def get_eph_function(self, x):
        """
        For a chromosome x, returns a function object eph to compute the ephemerides of the spacecraft

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Example:

          eph = prob.get_eph_function(population.champion_x)
          pos, vel = eph(pykep.epoch(7000))

        """
        if len(x) != len(self.get_bounds()[0]):
            raise ValueError(
                "Expected chromosome of length "
                + str(len(self.get_bounds()[0]))
                + " but got length "
                + str(len(x))
            )
        _, _, _, b_legs, b_ep = self._compute_dvs(x)

        def eph(
            t: float,
        ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:

            if t < b_ep[0]:
                raise ValueError(
                    "Given epoch " + str(t) + " is before launch date " + str(b_ep[0])
                )

            if t == b_ep[0]:
                # exactly at launch
                return self._seq[0].eph(t)

            i = _bisect_left(b_ep, t)  # ballistic leg i goes from planet i to planet i+1

            assert i >= 1 and i <= len(b_ep)
            if i < len(b_ep):
                assert t <= b_ep[i]

            # get start of ballistic leg
            r_b, v_b = b_legs[i - 1]

            elapsed_seconds = (t - b_ep[i - 1]) * _pk.DAY2SEC
            assert elapsed_seconds >= 0

            # propagate the lagrangian
            r, v = _pk.propagate_lagrangian(
                [r_b, v_b], elapsed_seconds, self._seq[0].mu_central_body
            )

            return r, v

        return eph

# Cleaning the namespace
del Any, Dict, List, Tuple