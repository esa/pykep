## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)##
## This file is part of the kep3 library.##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pykep as pk
import numpy as _np

from typing import Any, Dict, List, Tuple
from bisect import bisect_left


class mga:
    r"""The Multiple Gravity Assist (MGA) encoding of an interplanetary trajectory.

    This class may be used as a User Defined Problem (UDP) for the pygmo (http://esa.github.io/pygmo/) optimisation suite.

    The decision vector (chromosome) is::

      direct encoding: z = [t0, T1, T2 ... ] in [mjd2000, days, days ... ]
      alpha encoding: z = [t0, T, a1, a2 ...] in [mjd2000, days, nd, nd ... ]
      eta encoding: z = [t0, n1, n2, n3 ...] in [mjd2000, nd, nd ...]

    .. note::

       The time of flights of a MGA trajectory (and in general) can be encoded in different ways.
       When they are directly present in the decision vector, we talk about a *direct* encoding. This is the most
       'evolvable' encoding but also the one that requires the most problem knowledge (e.g. to define the bounds on each leg)
       and is not very flexible in dealing with constraints on the total time of flight. The *alpha* and *eta* encodings,
       instead, allow to only specify bounds on the time of flight of the entire trajectory, and not on the single legs:
       a property that is attractive for multi-objective optimization, for example.

       In the *alpha* encoding each leg time-of-flight is decoded as follows,

       :math:`T_i = T \log\alpha_i / \sum_n\log\alpha_n`.

       In the *eta* encoding  each leg time-of-flight is decoded as follows,

       :math:`T_i = (T_{max} - \sum_{j=0}^{(i-1)}T_j) \eta_i`

       The chromosome dimension for the direct and eta encoding is the same, while the alpha encoding requires one more gene.

    .. note::

       The resulting problem is box-bounded (unconstrained).
    """

    def __init__(
        self,
        seq=[
            pk.planet(pk.udpla.jpl_lp("earth")),
            pk.planet(pk.udpla.jpl_lp("venus")),
            pk.planet(pk.udpla.jpl_lp("earth")),
        ],
        t0=[0, 1000],
        tof=[[30, 200], [200, 300]],
        vinf=2.5,
        multi_objective=False,
        tof_encoding="direct",
        orbit_insertion=False,
        e_target=None,
        rp_target=None,
    ):
        r"""mga(seq, t0, tof, vinf, multi_objective=False, tof_encoding="direct"", orbit_insertion=False, e_target=None, rp_target=None)

        Args:
            *seq* (:class:`list` [:class:`~pykep.planet`]): sequence of planetary encounters including the departure body.

            *t0* (:class:`list` [:class:`float` or :class:`~pykep.epoch`]): lower and upper bounds for the launch epoch. When floats are used MJD2000 is assumed.

            *tof* (:class:`list` or :class:`float`): defines the bounds on the time of flight. If *tof_encoding* is 'direct', this contains a list
            of 2D lists defining the upper and lower bounds on each leg. If *tof_encoding* is 'alpha',
            this contains a 2D list with the lower and upper bounds on the total time-of-flight. If *tof_encoding*
            is 'eta' tof is a float defining an upper bound for the time-of-flight.

            *vinf* (:class:`float`): the vinf provided at launch for free

            *multi_objective* (:class:`bool`): when True constructs a multiobjective problem (dv, T). In this case, 'alpha' or `eta` encodings are recommended

            *tof_encoding* (:class:`str`): one of 'direct', 'alpha' or 'eta'. Selects the encoding for the time of flights

            *orbit_insertion* (:class:`bool`): when True the arrival dv is computed as that required to acquire a target orbit defined by e_target and rp_target

            *e_target* (:class:`float`): if orbit_insertion is True this defines the target orbit eccentricity around the final planet

            *rp_target* (:class:`float`): if orbit_insertion is True this defines the target orbit pericenter around the final planet (in m)

            *max_revs* (:class:`int`): maximal number of revolutions for lambert transfer
        """

        # Sanity checks
        # 1 - All planets need to have the same mu_central_body
        if [r.mu_central_body for r in seq].count(seq[0].mu_central_body) != len(seq):
            raise ValueError(
                "All planets in the sequence need to have exactly the same mu_central_body"
            )

        # 2 - We try to build epochs out of the t0 list (mjd2000 by default)
        for i in range(len(t0)):
            if type(t0[i]) != type(pk.epoch(0.0)):
                t0[i] = pk.epoch(t0[i], pk.epoch.julian_type.MJD2000)

        # 3 - Check the tof bounds
        if tof_encoding == "alpha":
            if type(tof) is not type([]):
                raise TypeError(
                    r"When the tof_encoding is 'alpha', the tof must be in the form [lb, ub]."
                )
            if len(tof) != 2:
                raise TypeError(
                    r"When the tof_encoding is 'alpha', the tof must be in the form [lb, ub]."
                )
        elif tof_encoding == "direct":
            if type(tof) is not type([]):
                raise TypeError(
                    r"The selected encoding is 'direct', the tof must be of type list."
                )
            if len(tof) != (len(seq) - 1):
                raise TypeError(
                    r"The selected encoding is 'direct', the length of the tof list must be exactly "
                    + str(len(seq) - 1)
                    + r" while it seems to be "
                    + str(len(tof))
                )
            if not all(isinstance(elem, list) and len(elem) == 2 for elem in tof):
                raise TypeError(
                    r"The selected encoding is 'direct', the content of the tof list must be 2-D lists."
                )
            if not all(elem[1] >= elem[0] for elem in tof):
                raise TypeError(
                    r"The selected encoding is 'direct', the content of the tof list must contain lower bounds < upper bounds "
                )
        elif tof_encoding == "eta":
            try:
                float(tof)
            except TypeError:
                raise TypeError(
                    r"When tof_encoding is \'eta\', the tof must be a float."
                )
        if not tof_encoding in ["alpha", "eta", "direct"]:
            raise ValueError("tof_encoding must be one of 'alpha', 'eta', 'direct'")

        # 4 - Check that if orbit insertion is selected e_target and r_p are
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

        # Public data members
        self.seq = seq
        self.t0 = t0
        self.tof = tof
        self.vinf = vinf * 1000
        self.multi_objective = multi_objective
        self.tof_encoding = tof_encoding
        self.orbit_insertion = orbit_insertion
        self.e_target = e_target
        self.rp_target = rp_target

        # Private data members
        self._n_legs = len(seq) - 1
        self._common_mu = seq[0].get_mu_central_body()

    def get_nobj(self):
        return self.multi_objective + 1

    def get_bounds(self):
        t0 = self.t0
        tof = self.tof
        n_legs = self._n_legs

        if self.tof_encoding == "alpha":
            # decision vector is  [t0, T, a1, a2, ....]
            lb = [t0[0].mjd2000, tof[0]] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000, tof[1]] + [1.0 - 1e-3] * (n_legs)
        elif self.tof_encoding == "direct":
            # decision vector is  [t0, T1, T2, T3, ... ]
            lb = [t0[0].mjd2000] + [it[0] for it in self.tof]
            ub = [t0[1].mjd2000] + [it[1] for it in self.tof]
        elif self.tof_encoding == "eta":
            # decision vector is  [t0, n1, n2, ....]
            lb = [t0[0].mjd2000] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000] + [1.0 - 1e-3] * (n_legs)
        return (lb, ub)

    def _decode_tofs(self, x: List[float]) -> List[float]:
        if self.tof_encoding == "alpha":
            # decision vector is  [t0, T, a1, a2, ....]
            return pk.alpha2direct(x[2:], x[1])
        elif self.tof_encoding == "direct":
            # decision vector is  [t0, T1, T2, T3, ... ]
            return x[1:]
        elif self.tof_encoding == "eta":
            # decision vector is  [t0, n1, n2, n3, ... ]
            return pk.eta2direct(x[1:], self.tof)

    @staticmethod
    def alpha2direct(x):
        """alpha2direct(x)

        Cpnverts the full MGA chromosome into a different encoding.

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the alpha encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the direct encoding
        """
        retval = pk.alpha2direct(x[2:], x[1])
        retval = _np.insert(retval, 0, x[0])
        return retval

    @staticmethod
    def direct2alpha(x):
        """direct2alpha(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the alpha encoding
        """
        alphas, T = pk.direct2alpha(x[1:])
        retval = _np.insert(alphas, 0, [x[0], T])
        return retval

    def eta2direct(self, x):
        """eta2direct(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the eta encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the direct encoding

        Raises:
            - ValueError: when the tof_encoding is not 'eta'
        """
        if self.tof_encoding != "eta":
            raise ValueError("cannot call this method if the tof_encoding is not 'eta'")

        # decision vector is  [t0, n1, n2, n3, ... ]
        T = pk.eta2direct(x[1:], self.tof)
        T = _np.insert(T, 0, x[0])
        return T

    def direct2eta(self, x):
        """direct2eta(x)

        Args:
            *x* (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            :class:`numpy.ndarray`: a chromosome encoding the MGA trajectory using the eta encoding

        Raises:
            - ValueError: when the tof_encoding is not 'eta'
        """
        if self.tof_encoding != "eta":
            raise ValueError("cannot call this method if the tof_encoding is not 'eta'")
        from copy import deepcopy

        retval = deepcopy(x)
        retval[1:] = pk.direct2eta(x[1:], self.tof)
        return retval

    def _compute_dvs(self, x: List[float]) -> Tuple[
        float,  # DVlaunch
        List[float],  # DVs
        float,  # DVarrival,
        List[Any],  # Lambert legs
        float,  # DVlaunch_tot
        List[float],  # T
        List[Tuple[List[float], List[float]]],  # ballistic legs
        List[float],  # epochs of ballistic legs
    ]:
        # 1 -  we 'decode' the times of flights and compute all epochs (mjd2000)
        # of the various planetary encounters
        T: List[float] = self._decode_tofs(x)  # T is [T1, T2 ...]
        ep = _np.insert(T, 0, x[0])  # T is [t0, T1, T2 ...]
        ep = _np.cumsum(ep)  # [t0, t1, t2, ...]

        # 2 - we compute the ephemerides
        r = [0] * len(self.seq)
        v = [0] * len(self.seq)
        for i in range(len(self.seq)):
            r[i], v[i] = self.seq[i].eph(ep[i])

        # 3 - we solve the lambert problems (and store trajectory r,v)
        lps = list()
        for i in range(self._n_legs):
            lp = pk.lambert_problem(
                r0=r[i],
                r1=r[i + 1],
                tof=T[i] * pk.DAY2SEC,
                mu=self._common_mu,
                cw=False,
                multi_revs=0,
            )
            lps.append(lp)

        # 4 - we compute the various dVs needed at fly-bys to match incoming
        # and outcoming
        DVfb = list()
        for i in range(len(lps) - 1):
            v_rel_in = [a - b for a, b in zip(lps[i].v1[0], v[i + 1])]
            v_rel_out = [a - b for a, b in zip(lps[i + 1].v0[0], v[i + 1])]
            DVfb.append(pk.fb_dv(v_rel_in, v_rel_out, self.seq[i + 1]))

        # 5 - we add the departure and arrival dVs
        DVlaunch_tot = _np.linalg.norm([a - b for a, b in zip(v[0], lps[0].v0[0])])
        DVlaunch = max(0, DVlaunch_tot - self.vinf)
        DVarrival = _np.linalg.norm([a - b for a, b in zip(v[-1], lps[-1].v1[0])])
        if self.orbit_insertion:
            # In this case we compute the insertion DV as a single pericenter
            # burn
            MU_SELF = self.seq[-1].get_mu_self()
            DVper = _np.sqrt(DVarrival * DVarrival + 2 * MU_SELF / self.rp_target)
            DVper2 = _np.sqrt(
                2 * MU_SELF / self.rp_target
                - MU_SELF / self.rp_target * (1.0 - self.e_target)
            )
            DVarrival = _np.abs(DVper - DVper2)
        return (DVlaunch, DVfb, DVarrival, lps, DVlaunch_tot, ep, T)

    # Objective function
    def fitness(self, x):
        DVlaunch, DVfb, DVarrival, _, _, _, _ = self._compute_dvs(x)
        if self.tof_encoding == "direct":
            T = sum(x[1:])
        elif self.tof_encoding == "alpha":
            T = x[1]
        elif self.tof_encoding == "eta":
            T = sum(self.eta2direct(x)[1:])
        if self.multi_objective:
            return [DVlaunch + _np.sum(DVfb) + DVarrival, T]
        else:
            return [DVlaunch + _np.sum(DVfb) + DVarrival]

    def to_planet(self, x: List[float]):
        """
        Returns a :class:`~pykep.planet` representing the spacecraft trajectory encoded into *x*.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.

        Example:
          sc = mga_udp.to_planet(population.champion_x)
          r, v = sc.eph(7000.)
        """

        class mga_udpla:
            def __init__(self, seq, lambert_legs, mjd2000s):
                self.seq = seq  # p0, p1, .., pn
                self.lambert_legs = lambert_legs  # l0, l1, .., l(n-1)
                self.mjd2000s = mjd2000s  # ep1, ep2, .., epn

            def eph(self, mjd2000: float):
                if mjd2000 < self.mjd2000s[0]:
                    raise ValueError(
                        "Ephemerides out-of-bounds, requested at mjd2000 "
                        + str(mjd2000)
                        + " which is before the launch date "
                        + str(self.mjd2000s[0])
                    )
                if mjd2000 > self.mjd2000s[-1]:
                    raise ValueError(
                        "Ephemerides out-of-bounds, requested at mjd2000 "
                        + str(mjd2000)
                        + " which is after the arrival date "
                        + str(self.mjd2000s[0])
                    )

                if mjd2000 == self.mjd2000s[0]:
                    # exactly at launch
                    return self.lambert_legs[0].r0, self.lambert_legs[0].v0[0]

                i = bisect_left(
                    self.mjd2000s, mjd2000
                )  # finds the corresponding leg from seq[i-1] to seq[i]

                # the starting conditions of the leg
                r0, v0 = self.lambert_legs[i - 1].r0, self.lambert_legs[i - 1].v0[0]
                elapsed_seconds = (mjd2000 - mjd2000s[i - 1]) * pk.DAY2SEC

                # propagate ballistically the starting conditions
                r1, v1 = pk.propagate_lagrangian(
                    rv=[r0, v0],
                    tof=elapsed_seconds,
                    mu=self.seq[0].get_mu_central_body(),
                    stm=False,
                )
                return r1, v1

            def get_name(self):
                return "MGA spacecraft"

        if len(x) != len(self.get_bounds()[0]):
            raise ValueError(
                "Expected chromosome of length "
                + str(len(self.get_bounds()[0]))
                + " but got length "
                + str(len(x))
            )

        _, _, _, lambert_legs, _, mjd2000s, _ = self._compute_dvs(x)
        return pk.planet(mga_udpla(self.seq, lambert_legs, mjd2000s))

    def pretty(self, x):
        """
        Prints a human readable representation of the transfer.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.
        """
        DVlaunch, DVfb, DVarrival, lambert_legs, DVlaunch_tot, mjd2000s, T = (
            self._compute_dvs(x)
        )
        print("Multiple Gravity Assist (MGA) problem: ")
        print("\tPlanet sequence: ", [pl.get_name() for pl in self.seq])
        print("\tEncoding for tofs: ", self.tof_encoding)
        print("\tOrbit Insertion: ", self.orbit_insertion)

        print("Departure: ", self.seq[0].get_name())
        print("\tEpoch: ", mjd2000s[0], " [mjd2000]")
        print("\tSpacecraft velocity: ", lambert_legs[0].v0[0], "[m/s]")
        print("\tHyperbolic velocity: ", DVlaunch_tot, "[m/s]")
        print("\tInitial DV: ", DVlaunch, "[m/s]")

        for pl, e, dv in zip(self.seq[1:-1], mjd2000s[1:-1], DVfb):
            print("Fly-by: ", pl.get_name())
            print("\tEpoch: ", e, " [mjd2000]")
            print("\tDV: ", dv, "[m/s]")

        print("Arrival: ", self.seq[-1].get_name())
        print("\tEpoch: ", mjd2000s[-1], " [mjd2000]")
        print("\tSpacecraft velocity: ", lambert_legs[-1].v1[0], "[m/s]")
        print("\tArrival DV: ", DVarrival, "[m/s]")

        print("Time of flights: ", T, "[days]")

        print("\nTotal DV: ", DVlaunch + _np.sum(DVfb) + DVarrival)

    def plot(
        self,
        x,
        ax=None,
        units=pk.AU,
        N=60,
        c_orbit="dimgray",
        c_lambert="indianred",
        leg_ids=[],
        figsize=(5, 5),
        **kwargs
    ):
        """
        Plots the trajectory encoded into *x* in 3D axes.

        Args:
            *x* (:class:`list`): The decision vector in the correct tof encoding.

            *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`, optional): The 3D axis to plot on. Defaults to None.

            *units* (:class:`float`, optional): The unit scale for the plot. Defaults to pk.AU.

            *N* (:class:`int`, optional): The number of points to use when plotting the trajectory. Defaults to 60.
            
            *c_orbit* (:class:`str`, optional): The color of the planet orbits. Defaults to 'dimgray'.
            
            *c* (:class:`str`, optional): The color of the trajectory. Defaults to 'indianred'.

            *figsize* (:class:`tuple`): The figure size (only used if a*ax* is None and axis have to be created.), Defaults to (5, 5).

            *leg_ids* (:class:`list`): selects the legs to plot. Optional, defaults to all legs.

            *\\*\\*kwargs*: Additional keyword arguments to pass to the trajectory plot (all Lambert arcs)

        Returns:
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: The 3D axis where the trajectory was plotted.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = pk.plot.make_3Daxis(figsize=figsize)
            
        # Plot of leg unless specified
        if len(leg_ids) == 0:
            leg_ids = list(range(self._n_legs))
            
        _, _, _, lps, _, mjd2000s, _ = self._compute_dvs(x)
        for i, item in enumerate(self.seq):
            if i in leg_ids:
                pk.plot.add_planet(pla=item, ax=ax, when=mjd2000s[i], c=c_orbit, units=units)
                pk.plot.add_lambert(
                    ax, lps[i], N=60, sol=0, units=units, c=c_lambert, **kwargs
                )
            pk.plot.add_planet_orbit(pla=item, ax=ax, units=units, N=N, c=c_orbit)
        pk.plot.add_sun(ax=ax)

        return ax


del Any, Dict, List, Tuple
