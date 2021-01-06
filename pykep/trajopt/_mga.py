from pykep.core import epoch, lambert_problem, propagate_lagrangian, DAY2SEC, fb_vel, AU
from pykep.planet import jpl_lp
from pykep.trajopt._lambert import lambert_problem_multirev

import numpy as np
from typing import Any, Dict, List, Tuple
from bisect import bisect_left


class mga:
    r"""
    This class transcribes a Multiple Gravity Assist (MGA) trajectory with no deep space maneuvers into an optimisation problem.
    It may be used as a User Defined Problem (UDP) for the pygmo (http://esa.github.io/pygmo/) optimisation suite.

    - Izzo, Dario. "Global optimization and space pruning for spacecraft trajectory design." Spacecraft Trajectory Optimization 1 (2010): 178-200.

    The decision vector (chromosome) is::

      direct encoding: [t0, T1, T2 ... ] in [mjd2000, days, days ... ]
      alpha encoding:  [t0, T, a1, a2 ...] in [mjd2000, days, nd, nd ... ]
      eta encoding:    [t0, n1, n2, n3 ...] in [mjd2000, nd, nd ...]

    .. note::

       The time of flights of a MGA trajectory (and in general) can be encoded in different ways.
       When they are directly present in the decision vector, we have the *direct* encoding. This is the most 'evolvable' encoding
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
                 seq=[jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')],
                 t0=[0, 1000],
                 tof=[[30, 200], [200, 300]],
                 vinf=2.5,
                 multi_objective=False,
                 tof_encoding='direct',
                 orbit_insertion=False,
                 e_target=None,
                 rp_target=None,
                 max_revs = 0
                 ):
        """mga(seq=[pk.planet.jpl_lp('earth'), pk.planet.jpl_lp('venus'), pk.planet.jpl_lp('earth')], t0=[0, 1000], tof=[100, 500], vinf=2.5, multi_objective=False, alpha_encoding=False, orbit_insertion=False, e_target=None, rp_target=None)

        Args:
            - seq (``list of pk.planet``): sequence of body encounters including the starting object
            - t0 (``list of pk.epoch``): the launch window
            - tof (``list`` or ``float``): bounds on the time of flight. If *tof_encoding* is 'direct', this contains a list
              of 2D lists defining the upper and lower bounds on each leg. If *tof_encoding* is 'alpha',
              this contains a 2D list with the lower and upper bounds on the time-of-flight. If *tof_encoding*
              is 'eta' tof is a float defining the upper bound on the time-of-flight
            - vinf (``float``): the vinf provided at launch for free
            - multi_objective (``bool``): when True constructs a multiobjective problem (dv, T). In this case, 'alpha' or `eta` encodings are recommended
            - tof_encoding (``str``): one of 'direct', 'alpha' or 'eta'. Selects the encoding for the time of flights
            - orbit_insertion (``bool``): when True the arrival dv is computed as that required to acquire a target orbit defined by e_target and rp_target
            - e_target (``float``): if orbit_insertion is True this defines the target orbit eccentricity around the final planet
            - rp_target (``float``): if orbit_insertion is True this defines the target orbit pericenter around the final planet (in m)
            - max_revs (``int``): maximal number of revolutions for lambert transfer

        Raises:
            - ValueError: if *planets* do not share the same central body (checked on the mu_central_body attribute)
            - ValueError: if *t0* does not contain objects able to construct a epoch (e.g. pk. epoch or floats)
            - ValueError: if *tof* is badly defined
            - ValueError: it the target orbit is not defined and *orbit_insertion* is True
        """

        # Sanity checks
        # 1 - All planets need to have the same mu_central_body
        if ([r.mu_central_body for r in seq].count(seq[0].mu_central_body) != len(seq)):
            raise ValueError(
                'All planets in the sequence need to have exactly the same mu_central_body')
        # 2 - We try to build epochs out of the t0 list (mjd2000 by default)
        for i in range(len(t0)):
            if (type(t0[i]) != type(epoch(0))):
                t0[i] = epoch(t0[i])
        # 3 - Check the tof bounds
        if tof_encoding == 'alpha':
            if len(tof) != 2:
                raise ValueError(
                    r'When the tof_encoding is \'alpha\', tof is expected to be something like [lb, ub]')
        elif tof_encoding == 'direct':
            if len(tof) != (len(seq) - 1):
                raise ValueError(
                    'When tof_encoding is direct, the tof must be a float (upper bound on the time of flight)' + str(len(seq) - 1))
        elif tof_encoding == 'eta':
            try:
                float(tof)
            except TypeError:
                raise ValueError(
                    'The tof needs to be have len equal to  ' + str(len(seq) - 1))
        if not tof_encoding in ['alpha', 'eta', 'direct']:
            raise ValueError(
                "tof_encoding must be one of 'alpha', 'eta', 'direct'")

        # 4 - Check that if orbit insertion is selected e_target and r_p are
        # defined
        if orbit_insertion:
            if rp_target is None:
                raise ValueError(
                    'The rp_target needs to be specified when orbit insertion is selected')
            if e_target is None:
                raise ValueError(
                    'The e_target needs to be specified when orbit insertion is selected')

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
        self._common_mu = seq[0].mu_central_body
        self.max_revs = max_revs

    def get_nobj(self):
        return self.multi_objective + 1

    def get_bounds(self):
        t0 = self.t0
        tof = self.tof
        n_legs = self._n_legs

        if self.tof_encoding == 'alpha':
            # decision vector is  [t0, T, a1, a2, ....]
            lb = [t0[0].mjd2000, tof[0]] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000, tof[1]] + [1.0 - 1e-3] * (n_legs)
        elif self.tof_encoding == 'direct':
            # decision vector is  [t0, T1, T2, T3, ... ]
            lb = [t0[0].mjd2000] + [it[0] for it in self.tof]
            ub = [t0[1].mjd2000] + [it[1] for it in self.tof]
        elif self.tof_encoding == 'eta':
            # decision vector is  [t0, n1, n2, ....]
            lb = [t0[0].mjd2000] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000] + [1.0 - 1e-3] * (n_legs)
        return (lb, ub)

    def _decode_tofs(self, x: List[float]) -> List[float]:
        if self.tof_encoding == 'alpha':
            # decision vector is  [t0, T, a1, a2, ....]
            T = np.log(x[2:])
            return T / sum(T) * x[1]
        elif self.tof_encoding == 'direct':
            # decision vector is  [t0, T1, T2, T3, ... ]
            return x[1:]
        elif self.tof_encoding == 'eta':
            # decision vector is  [t0, n1, n2, n3, ... ]
            dt = self.tof
            T = [0] * self._n_legs
            T[0] = dt * x[1]
            for i in range(1, len(T)):
                T[i] = (dt - sum(T[:i])) * x[i + 1]
            return T

    def alpha2direct(self, x):
        """alpha2direct(x)

        Args:
            - x (``array-like``): a chromosome encoding an MGA trajectory in the alpha encoding

        Returns:
            ``numpy.array``: a chromosome encoding the MGA trajectory using the direct encoding
        """
        T = np.log(x[2:])
        retval = T / sum(T) * x[1]
        retval = np.insert(retval, 0, x[0])
        return retval

    def direct2alpha(self, x):
        """direct2alpha(x)

        Args:
            - x (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            ``numpy.array``: a chromosome encoding the MGA trajectory using the alpha encoding
        """
        T = np.sum(x[1:])
        alphas = np.exp(x[1:] / (-T))
        retval = np.insert(alphas, 0, [x[0], T])
        return retval

    def eta2direct(self, x):
        """eta2direct(x)

        Args:
            - x (``array-like``): a chromosome encoding an MGA trajectory in the eta encoding

        Returns:
            ``numpy.array``: a chromosome encoding the MGA trajectory using the direct encoding

        Raises:
            - ValueError: when the tof_encoding is not 'eta'
        """
        if self.tof_encoding != 'eta':
            raise ValueError(
                "cannot call this method if the tof_encoding is not 'eta'")

        # decision vector is  [t0, n1, n2, n3, ... ]
        n = len(x) - 1
        dt = self.tof
        T = [0] * n
        T[0] = dt * x[1]
        for i in range(1, len(T)):
            T[i] = (dt - sum(T[:i])) * x[i + 1]
        np.insert(T, 0, [0])
        return T

    def direct2eta(self, x):
        """direct2eta(x)

        Args:
            - x (``array-like``): a chromosome encoding an MGA trajectory in the direct encoding

        Returns:
            ``numpy.array``: a chromosome encoding the MGA trajectory using the eta encoding

        Raises:
            - ValueError: when the tof_encoding is not 'eta'
        """
        if self.tof_encoding != 'eta':
            raise ValueError(
                "cannot call this method if the tof_encoding is not 'eta'")
        from copy import deepcopy
        retval = deepcopy(x)
        retval[1] = x[1] / self.tof
        for i in range(2, len(x)):
            retval[i] = x[i] / (self.tof - sum(x[1:i]))
        return retval

    def _compute_dvs(self, x: List[float]) -> Tuple[
        float, # DVlaunch
        List[float], # DVs
        float, # DVarrival,
        List[Any], # Lambert legs
        float, #DVlaunch_tot
        List[float], # T
        List[Tuple[List[float], List[float]]], # ballistic legs
        List[float], # epochs of ballistic legs
    ]:
        # 1 -  we 'decode' the times of flights and compute epochs (mjd2000)
        T: List[float] = self._decode_tofs(x)  # [T1, T2 ...]
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        # 2 - we compute the ephemerides
        r = [0] * len(self.seq)
        v = [0] * len(self.seq)
        for i in range(len(self.seq)):
            r[i], v[i] = self.seq[i].eph(float(ep[i]))

        l = list()
        ballistic_legs: List[Tuple[List[float],List[float]]] = []
        ballistic_ep: List[float] = []

        # 3 - we solve the lambert problems
        vi = v[0]
        for i in range(self._n_legs):
            lp = lambert_problem_multirev(
                vi, lambert_problem(
                    r[i], r[i + 1], T[i] * DAY2SEC, self._common_mu, False, self.max_revs))
            l.append(lp)
            vi = lp.get_v2()[0]
            ballistic_legs.append((r[i], lp.get_v1()[0]))
            ballistic_ep.append(ep[i])
        # 4 - we compute the various dVs needed at fly-bys to match incoming
        # and outcoming
        DVfb = list()
        for i in range(len(l) - 1):
            vin = [a - b for a, b in zip(l[i].get_v2()[0], v[i + 1])]
            vout = [a - b for a, b in zip(l[i + 1].get_v1()[0], v[i + 1])]
            DVfb.append(fb_vel(vin, vout, self.seq[i + 1]))
        # 5 - we add the departure and arrival dVs
        DVlaunch_tot = np.linalg.norm(
            [a - b for a, b in zip(v[0], l[0].get_v1()[0])])
        DVlaunch = max(0, DVlaunch_tot - self.vinf)
        DVarrival = np.linalg.norm(
            [a - b for a, b in zip(v[-1], l[-1].get_v2()[0])])
        if self.orbit_insertion:
            # In this case we compute the insertion DV as a single pericenter
            # burn
            DVper = np.sqrt(DVarrival * DVarrival + 2 *
                            self.seq[-1].mu_self / self.rp_target)
            DVper2 = np.sqrt(2 * self.seq[-1].mu_self / self.rp_target -
                             self.seq[-1].mu_self / self.rp_target * (1. - self.e_target))
            DVarrival = np.abs(DVper - DVper2)
        return (DVlaunch, DVfb, DVarrival, l, DVlaunch_tot, T, ballistic_legs, ballistic_ep)

    # Objective function
    def fitness(self, x):
        DVlaunch, DVfb, DVarrival, _, _, _, _, _ = self._compute_dvs(x)
        if self.tof_encoding == 'direct':
            T = sum(x[1:])
        elif self.tof_encoding == 'alpha':
            T = x[1]
        elif self.tof_encoding == 'eta':
            T = sum(self.eta2direct(x)[1:])
        if self.multi_objective:
            return [DVlaunch + np.sum(DVfb) + DVarrival, T]
        else:
            return [DVlaunch + np.sum(DVfb) + DVarrival]

    def pretty(self, x):
        """pretty(x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Prints human readable information on the trajectory represented by the decision vector x
        """
        T = self._decode_tofs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        DVlaunch, DVfb, DVarrival, l, DVlaunch_tot, _, _, _ = self._compute_dvs(x)
        print("Multiple Gravity Assist (MGA) problem: ")
        print("Planet sequence: ", [pl.name for pl in self.seq])

        print("Departure: ", self.seq[0].name)
        print("\tEpoch: ", ep[0], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[0].get_v1()[0], "[m/s]")
        print("\tHyperbolic velocity: ", DVlaunch_tot, "[m/s]")
        print("\tInitial DV: ", DVlaunch, "[m/s]")

        for pl, e, dv in zip(self.seq[1:-1], ep[1:-1], DVfb):
            print("Fly-by: ", pl.name)
            print("\tEpoch: ", e, " [mjd2000]")
            print("\tDV: ", dv, "[m/s]")

        print("Arrival: ", self.seq[-1].name)
        print("\tEpoch: ", ep[-1], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[-1].get_v2()[0], "[m/s]")
        print("\tArrival DV: ", DVarrival, "[m/s]")

        print("Time of flights: ", T, "[days]")

    def plot(self, x, axes=None, units=AU, N=60):
        """plot(self, x, axes=None, units=pk.AU, N=60)

        Plots the spacecraft trajectory.

        Args:
            - x (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.
            - N (``float``): Number of points to plot per leg
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        from pykep.orbit_plots import plot_planet, plot_lambert

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        _, _, _, l, _, T, _, _ = self._compute_dvs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        
        for pl, e in zip(self.seq, ep):
            plot_planet(pl, epoch(e), units=units, legend=True,
                        color=(0.7, 0.7, 1), axes=axes)
        for lamb in l:
            plot_lambert(lamb, N=N, sol=0, units=units, color='k',
                         legend=False, axes=axes, alpha=0.8)
        return axes

    def get_eph_function(self, x: List[float]):
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

        _, _, _, _, _, _, b_legs, b_ep = self._compute_dvs(x)
        
        def eph(
            t: float
        ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:

            if t < b_ep[0]:
                raise ValueError(
                    "Given epoch " + str(t) + " is before launch date " + str(b_ep[0])
                )

            if t == b_ep[0]:
                # exactly at launch
                return self.seq[0].eph(t)

            i = bisect_left(b_ep, t)  # ballistic leg i goes from planet i to planet i+1

            assert i >= 1 and i <= len(b_ep)
            if i < len(b_ep):
                assert t <= b_ep[i]

            # get start of ballistic leg
            r_b, v_b = b_legs[i - 1]

            elapsed_seconds = (t - b_ep[i - 1]) * DAY2SEC
            assert elapsed_seconds >= 0

            # propagate the lagrangian
            r, v = propagate_lagrangian(r_b, v_b, elapsed_seconds, self.seq[0].mu_central_body)

            return r, v
        
        return eph



if __name__ == "__main__":
    import pygmo as pg
    seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp(
        'venus'), jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('saturn')]
    udp = mga(seq=seq,
              t0=[-1000., 0.],
              tof=[4000., 7000.],
              vinf=3.,
              tof_encoding='alpha',
              orbit_insertion=True,
              e_target=0.98,
              rp_target=108950000)

    udp = mga(seq=seq,
              t0=[-1000., 0.],
              tof=[[30, 400], [100, 470], [30, 400], [400, 2000], [1000, 6000]],
              vinf=3.,
              tof_encoding='direct',
              orbit_insertion=True,
              e_target=0.98,
              rp_target=108950000)

    udp = mga(seq=seq,
              t0=[-1000., 0.],
              tof=7000.,
              vinf=3.,
              tof_encoding='eta',
              orbit_insertion=True,
              e_target=0.98,
              rp_target=108950000)

    #udp = mga(seq=seq, t0=[-1000., 0.], tof=[[130,200], [430,470], [30, 70], [900, 1200], [4000, 5000]], vinf=0., alpha_encoding=False)
    prob = pg.problem(udp)
    uda = pg.cmaes(1500, force_bounds=True, sigma0=0.5, ftol=1e-4)
    #uda = pg.sade(4500)
    algo = pg.algorithm(uda)
    algo.set_verbosity(10)
    res = list()
    # for i in range(100):
    pop = pg.population(prob, 100)
    pop = algo.evolve(pop)
    res.append(pop.champion_f)
