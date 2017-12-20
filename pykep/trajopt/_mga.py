import pykep as pk
import numpy as np


class mga:
    """
    This class is a pygmo (http://esa.github.io/pygmo/) problem representing a Multiple Gravity Assist
    trajectory with no deep space manouvres.

    - Izzo, Dario. "Global optimization and space pruning for spacecraft trajectory design." Spacecraft Trajectory Optimization 1 (2010): 178-200.

    The decision vector (chromosome) is::

      direct encoding: [t0, T1, T2 ... ] in [mjd2000, days, days ... ]
      alpha encoding: [t0, T, a1, a2 ...] in [mjd2000, days, nd, nd ... ]

    .. note::

       The time of flights of the MGA trajectory can be encoded in two different ways: directly or using the alpha encoding.
       The alpha encoding (introduced by our previous research), allows to only specify bounds on the entire trajectory, and not
       on the single legs and was introduced for multi-objective cases where the total time-of-flight is the second objective
       and benefits from being box bounded.
       In the alpha encoding each leg time-of-flight is decoded as follows, T_n = T log(alpha_n) / \sum_i(log(alpha_i)).

    .. note::

       The resulting problem is box-bounded (unconstrained).
    """

    def __init__(self,
                 seq=[pk.planet.jpl_lp('earth'), pk.planet.jpl_lp(
                     'venus'), pk.planet.jpl_lp('earth')],
                 t0=[pk.epoch(0), pk.epoch(1000)],
                 tof=[100, 500],
                 vinf=2.5,
                 multi_objective=False,
                 alpha_encoding=False,
                 orbit_insertion=False,
                 e_target=None,
                 rp_target=None
                 ):
        """
        pk.trajopt.mga(seq=[pk.planet.jpl_lp('earth'), pk.planet.jpl_lp('venus'), pk.planet.jpl_lp('earth')], t0=[pk.epoch(0), pk.epoch(1000)], tof=[100, 500], vinf=2.5, multi_objective=False, alpha_encoding=False, orbit_insertion=False, e_target=None, rp_target=None)

        Args:
            - seq (``list of pk.planet``): sequence of body encounters including the starting object
            - t0 (``list of pk.epoch``): the launch window
            - tof (``list``): list of pairs defining lower and upper bounds ind days on the various legs. If alpha_encoding
                 is True, then the list contains only two floats defining the upper and lower bounds on the total tof 
            - vinf (``float``): the vinf provided at launch for free
            - multi_objective (``bool``): when True constructs a multiobjective problem (dv, T). In this case alpha encoding is reccomended
            - orbit_insertion (``bool``): when True the arrival dv is computed as that required to acquire a target orbit defined by e_target and rp_target
            - e_target (``float``): if orbit_insertion is True this defines the target orbit eccentricity around the final planet
            - rp_target (``float``): if orbit_insertion is True this defines the target orbit pericenter around the final planet

        Raises:
            - ValueError: if *planets* do not share the same central body (checked on the mu_central_body attribute)
            - ValueError: if *t0* does not contain objects able to construct a pk.epoch (e.g. pk. epoch or floats)
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
            if (type(t0[i]) != type(pk.epoch(0))):
                t0[i] = pk.epoch(t0[i])
        # 3 - Check the tof bounds
        if alpha_encoding:
            if len(tof) != 2:
                raise ValueError(
                    'The tof needs to be have len equal to 2')
        else:
            if len(tof) != (len(seq) - 1):
                raise ValueError(
                    'The tof needs to be have len equal to  ' + str(len(seq) - 1))
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
        self.alpha_encoding = alpha_encoding
        self.orbit_insertion = orbit_insertion
        self.e_target = e_target
        self.rp_target = rp_target

        # Private data members
        self._n_legs = len(seq) - 1
        self._common_mu = seq[0].mu_central_body

    def get_nobj(self):
        return self.multi_objective + 1

    def get_bounds(self):
        t0 = self.t0
        tof = self.tof
        seq = self.seq
        n_legs = self._n_legs

        if self.alpha_encoding:
            # decision vector is  [t0, T, a1, a2, ....]
            lb = [t0[0].mjd2000, tof[0]] + [1e-3] * (n_legs)
            ub = [t0[1].mjd2000, tof[1]] + [1.0 - 1e-3] * (n_legs)
        else:
            # decision vector is  [t0, T1, T2, T3, ... ]
            lb = [t0[0].mjd2000] + [it[0] for it in self.tof]
            ub = [t0[1].mjd2000] + [it[1] for it in self.tof]
        return (lb, ub)

    def _decode_tofs(self, x):
        if self.alpha_encoding:
            # decision vector is  [t0, T, a1, a2, ....]
            T = np.log(x[2:])
            return T / sum(T) * x[1]
        else:
            # decision vector is  [t0, T1, T2, T3, ... ]
            return x[1:]

    def alpha2direct(self, x):
        T = np.log(x[2:])
        retval = T / sum(T) * x[1]
        retval = np.insert(retval, 0, x[0])
        return retval

    def direct2alpha(self, x):
        T = np.sum(x[1:])
        alphas = np.exp(x[1:] / (-T))
        retval = np.insert(alphas,0, [x[0], T])
        return retval

    def _compute_dvs(self, x):
        # 1 -  we 'decode' the times of flights and compute epochs (mjd2000)
        T = self._decode_tofs(x)  # [T1, T2 ...]
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        # 2 - we compute the ephemerides
        r = [0] * len(self.seq)
        v = [0] * len(self.seq)
        for i in range(len(self.seq)):
            r[i], v[i] = self.seq[i].eph(ep[i])
        # 3 - we solve the lambert problems
        l = list()
        for i in range(self._n_legs):
            l.append(pk.lambert_problem(
                r[i], r[i + 1], T[i] * pk.DAY2SEC, self._common_mu, False, 0))
        # 4 - we compute the various dVs needed at fly-bys to match incoming
        # and outcoming
        DVfb = list()
        for i in range(len(l) - 1):
            vin = [a - b for a, b in zip(l[i].get_v2()[0], v[i + 1])]
            vout = [a - b for a, b in zip(l[i + 1].get_v1()[0], v[i + 1])]
            DVfb.append(pk.fb_vel(vin, vout, self.seq[i + 1]))
        # 5 - we add the departure and arrival dVs
        DVlaunch_tot = np.linalg.norm(
            [a - b for a, b in zip(v[0], l[0].get_v1()[0])])
        DVlaunch = max(0,DVlaunch_tot - self.vinf)
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
        return (DVlaunch, DVfb, DVarrival, l, DVlaunch_tot)

    # Objective function
    def fitness(self, x):
        DVlaunch, DVfb, DVarrival, _, _= self._compute_dvs(x)
        if self.multi_objective:
            return [DVlaunch + np.sum(DVfb) + DVarrival, x[1]]
        else:
            return [DVlaunch + np.sum(DVfb) + DVarrival]

    def pretty(self, x):
        """prob.pretty(x)

        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).

        Prints human readable information on the trajectory represented by the decision vector x
        """
        T = self._decode_tofs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        DVlaunch, DVfb, DVarrival, l, DVlaunch_tot = self._compute_dvs(x)
        print("Multiple Gravity Assist (MGA) problem: ")
        print("Planet sequence: ", [pl.name for pl in self.seq])

        print("Departure: ", seq[0].name)
        print("\tEpoch: ", ep[0], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[0].get_v1()[0], "[m/s]")
        print("\tHyperbolic velocity: ", DVlaunch_tot, "[m/s]")
        print("\tInitial DV: ", DVlaunch, "[m/s]")

        for pl,e, dv in zip(seq[1:-1], ep[1:-1], DVfb):
            print("Fly-by: ", pl.name)
            print("\tEpoch: ", e, " [mjd2000]")
            print("\tDV: ", dv, "[m/s]")

        print("Arrival: ", seq[-1].name)
        print("\tEpoch: ", ep[-1], " [mjd2000]")
        print("\tSpacecraft velocity: ", l[-1].get_v2()[0], "[m/s]")
        print("\tArrival DV: ", DVarrival, "[m/s]")

        print("Time of flights: ", T, "[days]")

    def plot(self, x, axes=None, units=pk.AU, N=60):
        """Plots the spacecraft trajectory.

        Args:
            - x (``tuple``, ``list``, ``numpy.ndarray``): Decision chromosome.
            - axes (``matplotlib.axes._subplots.Axes3DSubplot``): 3D axes to use for the plot
            - units (``float``, ``int``): Length unit by which to normalise data.
            - N (``float``): Number of points to plot per leg

        Examples:
            >>> prob.extract(pykep.trajopt.indirect_pt2or).plot_traj(pop.champion_x)
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        T = self._decode_tofs(x)
        ep = np.insert(T, 0, x[0])  # [t0, T1, T2 ...]
        ep = np.cumsum(ep)  # [t0, t1, t2, ...]
        DVlaunch, DVfb, DVarrival, l, _ = self._compute_dvs(x)
        for pl, e in zip(self.seq, ep):
            pk.orbit_plots.plot_planet(pl, pk.epoch(e), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)
        for lamb in l:
            pk.orbit_plots.plot_lambert(
                lamb, N=N, sol=0, units=units, color='k', legend=False, ax=axes, alpha=0.8)
        return axes

if __name__ == "__main__":
    import pygmo as pg
    seq = [pk.planet.jpl_lp('earth'), pk.planet.jpl_lp('venus'), pk.planet.jpl_lp(
        'venus'), pk.planet.jpl_lp('earth'), pk.planet.jpl_lp('jupiter'), pk.planet.jpl_lp('saturn')]
    udp = mga(seq=seq,
              t0=[-1000., 0.],
              tof=[4000., 7000.],
              vinf=0.,
              alpha_encoding=True,
              orbit_insertion=True,
              e_target=0.98,
              rp_target=108950000)
    
    udp = mga(seq=seq, 
            t0=[-1000., 0.], 
            tof=[[30,400], [100,470], [30, 400], [400, 2000], [1000, 6000]], 
            vinf=3., 
            alpha_encoding=False,
            orbit_insertion=True,
            e_target = 0.98,
            rp_target = 108950000)   
    
    #udp = mga(seq=seq, t0=[-1000., 0.], tof=[[130,200], [430,470], [30, 70], [900, 1200], [4000, 5000]], vinf=0., alpha_encoding=False)
    prob = pg.problem(udp)
    #uda = pg.cmaes(1500, force_bounds = True, sigma0 = 0.5, ftol = 1e-4)
    uda = pg.sade(2500)
    algo = pg.algorithm(uda)
    algo.set_verbosity(10)
    res = list()
    for i in range(100):
        pop = pg.population(prob, 20)
        pop = algo.evolve(pop)
        res.append(pop.champion_f)
