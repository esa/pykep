from pykep.core import epoch, DAY2SEC, lambert_problem, propagate_lagrangian, SEC2DAY, AU, ic2par
from pykep.planet import jpl_lp
from math import pi, cos, sin, log, acos
from scipy.linalg import norm


class pl2pl_N_impulses(object):
    """
    This class is a pygmo (http://esa.github.io/pygmo/) problem representing a single leg transfer
    between two planets allowing up to a maximum number of impulsive Deep Space Manouvres.

    The decision vector is::

      [t0,T] + [alpha,u,v,V_inf] * (N-2) + [alpha] + ([tf])

    ... in the units: [mjd2000, days] + [nd, nd, m/sec, nd] + [nd] + [mjd2000]

    Each leg time-of-flight can be decoded as follows, T_n = T log(alpha_n) / \sum_i(log(alpha_i))

    .. note::

       The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
    """

    def __init__(self,
                 start=jpl_lp('earth'),
                 target=jpl_lp('venus'),
                 N_max=3,
                 tof=[20., 400.],
                 vinf=[0., 4.],
                 phase_free=True,
                 multi_objective=False,
                 t0=None
                 ):
        """
        prob = pykep.trajopt.pl2pl_N_impulses(start=jpl_lp('earth'), target=jpl_lp('venus'), N_max=3, tof=[20., 400.], vinf=[0., 4.], phase_free=True, multi_objective=False, t0=None)

        Args: 
            - start (``pykep.planet``): the starting planet
            - target (``pykep.planet``): the target planet
            - N_max (``int``): maximum number of impulses
            - tof (``list``): the box bounds [lower,upper] for the time of flight (days)
            - vinf (``list``): the box bounds [lower,upper] for each DV magnitude (km/sec)
            - phase_free (``bool``): when True, no randezvous condition are enforced and start and arrival anomalies will be free
            - multi_objective (``bool``):  when True, a multi-objective problem is constructed with DV and time of flight as objectives
            - t0 (``list``):  the box bounds on the launch window containing two pykep.epoch. This is not needed if phase_free is True.
        """

        # Sanity checks
        # 1) all planets need to have the same mu_central_body
        if (start.mu_central_body != target.mu_central_body):
            raise ValueError('Starting and ending pykep.planet must have the same mu_central_body')
        # 2) Number of impulses must be at least 2
        if N_max < 2:
            raise ValueError('Number of impulses N is less than 2')
        # 3) If phase_free is True, t0 does not make sense
        if (t0 is None and not phase_free):
            t0 = [epoch(0), epoch(1000)]
        if (t0 is not None and phase_free):
            raise ValueError('When phase_free is True no t0 can be specified')

        self.obj_dim = multi_objective + 1
        # We then define all class data members
        self.start = start
        self.target = target
        self.N_max = N_max
        self.phase_free = phase_free
        self.multi_objective = multi_objective
        self.vinf = [s*1000 for s in vinf]

        self.__common_mu = start.mu_central_body
        
        # And we compute the bounds
        if phase_free:
            self._lb = [0, tof[0]] + [0.0, 0.0, 0.0, vinf[0] * 1000] * (N_max - 2) + [0.0] + [0]
            self._ub = [2 * start.compute_period(epoch(0)) * SEC2DAY, tof[1]] + [1.0, 1.0, 1.0, vinf[1] * 1000] * (N_max - 2) + [1.0] + [2 * target.compute_period(epoch(0)) * SEC2DAY]
        else:
            self._lb = [t0[0].mjd2000, tof[0]] + [0.0, 0.0, 0.0, vinf[0] * 1000] * (N_max - 2) + [0.0]
            self._ub = [t0[1].mjd2000, tof[1]] + [1.0, 1.0, 1.0, vinf[1] * 1000] * (N_max - 2) + [1.0]

    def get_nobj(self):
        return self.obj_dim

    def get_bounds(self):
        return (self._lb, self._ub)

    def fitness(self, x):
        # 1 -  we 'decode' the chromosome recording the various deep space
        # manouvres timing (days) in the list T
        T = list([0] * (self.N_max - 1))

        for i in range(len(T)):
            T[i] = log(x[2 + 4 * i])
        total = sum(T)
        T = [x[1] * time / total for time in T]

        # 2 - We compute the starting and ending position
        r_start, v_start = self.start.eph(epoch(x[0]))
        if self.phase_free:
            r_target, v_target = self.target.eph(epoch(x[-1]))
        else:
            r_target, v_target = self.target.eph(epoch(x[0] + x[1]))

        # 3 - We loop across inner impulses
        rsc = r_start
        vsc = v_start
        for i, time in enumerate(T[:-1]):
            theta = 2 * pi * x[3 + 4 * i]
            phi = acos(2 * x[4 + 4 * i] - 1) - pi / 2

            Vinfx = x[5 + 4 * i] * cos(phi) * cos(theta)
            Vinfy = x[5 + 4 * i] * cos(phi) * sin(theta)
            Vinfz = x[5 + 4 * i] * sin(phi)

            # We apply the (i+1)-th impulse
            vsc = [a + b for a, b in zip(vsc, [Vinfx, Vinfy, Vinfz])]
            rsc, vsc = propagate_lagrangian(rsc, vsc, T[i] * DAY2SEC, self.__common_mu)
        cw = (ic2par(rsc, vsc, self.start.mu_central_body)[2] > pi / 2)

        # We now compute the remaining two final impulses
        # Lambert arc to reach seq[1]
        dt = T[-1] * DAY2SEC
        l = lambert_problem(rsc, r_target, dt, self.__common_mu, cw, False)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        DV1 = norm([a - b for a, b in zip(v_beg_l, vsc)])
        DV2 = norm([a - b for a, b in zip(v_end_l, v_target)])

        DV_others = sum(x[5::4])
        if self.obj_dim == 1:
            return (DV1 + DV2 + DV_others,)
        else:
            return (DV1 + DV2 + DV_others, x[1])

    def plot_trajectory(self, x, axes=None):
        """
        ax = prob.plot_trajectory(x, axes=None)

        - x: encoded trajectory
        - axes: matplotlib axis where to plot. If None figure and axis will be created
        - [out] ax: matplotlib axis where to plot

        Plots the trajectory represented by a decision vector x on the 3d axis ax

        Example::

          ax = prob.plot(x)
        """
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from pykep.orbit_plots import plot_planet, plot_lambert, plot_kepler

        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')


        axes.scatter(0, 0, 0, color='y')

        # 1 -  we 'decode' the chromosome recording the various deep space
        # manouvres timing (days) in the list T
        T = list([0] * (self.N_max - 1))

        for i in range(len(T)):
            T[i] = log(x[2 + 4 * i])
        total = sum(T)
        T = [x[1] * time / total for time in T]

        # 2 - We compute the starting and ending position
        r_start, v_start = self.start.eph(epoch(x[0]))
        if self.phase_free:
            r_target, v_target = self.target.eph(epoch(x[-1]))
        else:
            r_target, v_target = self.target.eph(epoch(x[0] + x[1]))
        plot_planet(self.start, t0=epoch(x[0]), color=(0.8, 0.6, 0.8), legend=True, units = AU, ax=axes, s=0)
        plot_planet(self.target, t0=epoch(x[0] + x[1]), color=(0.8, 0.6, 0.8), legend=True, units = AU, ax=axes, s=0)

        DV_list = x[5::4]
        maxDV = max(DV_list)
        DV_list = [s / maxDV * 30 for s in DV_list]
        colors = ['b','r'] * (len(DV_list) + 1)

        # 3 - We loop across inner impulses
        rsc = r_start
        vsc = v_start
        for i, time in enumerate(T[:-1]):
            theta = 2 * pi * x[3 + 4 * i]
            phi = acos(2 * x[4 + 4 * i] - 1) - pi / 2

            Vinfx = x[5 + 4 * i] * cos(phi) * cos(theta)
            Vinfy = x[5 + 4 * i] * cos(phi) * sin(theta)
            Vinfz = x[5 + 4 * i] * sin(phi)

            # We apply the (i+1)-th impulse
            vsc = [a + b for a, b in zip(vsc, [Vinfx, Vinfy, Vinfz])]
            axes.scatter(rsc[0]/AU, rsc[1]/AU, rsc[2]/AU, color='k', s = DV_list[i])
            plot_kepler(rsc, vsc, T[i] * DAY2SEC, self.__common_mu, N=200, color=colors[i], legend=False, units=AU, ax=axes)
            rsc, vsc = propagate_lagrangian(rsc, vsc, T[i] * DAY2SEC, self.__common_mu)

        cw = (ic2par(rsc, vsc, self.start.mu_central_body)[2] > pi / 2)
        # We now compute the remaining two final impulses
        # Lambert arc to reach seq[1]
        dt = T[-1] * DAY2SEC
        l = lambert_problem(rsc, r_target, dt, self.__common_mu, cw, False)
        plot_lambert(l, sol=0, color=colors[i+1], legend=False, units=AU, ax=axes, N=200)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]
        DV1 = norm([a - b for a, b in zip(v_beg_l, vsc)])
        DV2 = norm([a - b for a, b in zip(v_end_l, v_target)])

        axes.scatter(rsc[0]/AU, rsc[1]/AU, rsc[2]/AU, color='k', s = min(DV1/ maxDV *30, 40))
        axes.scatter(r_target[0]/AU, r_target[1]/AU, r_target[2]/AU, color='k', s = min(DV2/ maxDV *30, 40))



        return axes

    def pretty(self, x):
        # 1 -  we 'decode' the chromosome recording the various deep space
        # manouvres timing (days) in the list T
        T = list([0] * (self.N_max - 1))

        for i in range(len(T)):
            T[i] = log(x[2 + 4 * i])
        total = sum(T)
        T = [x[1] * time / total for time in T]

        # 2 - We compute the starting and ending position
        r_start, v_start = self.start.eph(epoch(x[0]))
        if self.phase_free:
            r_target, v_target = self.target.eph(epoch(x[-1]))
        else:
            r_target, v_target = self.target.eph(epoch(x[0] + x[1]))

        # 3 - We loop across inner impulses
        rsc = r_start
        vsc = v_start
        for i, time in enumerate(T[:-1]):
            theta = 2 * pi * x[3 + 4 * i]
            phi = acos(2 * x[4 + 4 * i] - 1) - pi / 2

            Vinfx = x[5 + 4 * i] * cos(phi) * cos(theta)
            Vinfy = x[5 + 4 * i] * cos(phi) * sin(theta)
            Vinfz = x[5 + 4 * i] * sin(phi)

            # We apply the (i+1)-th impulse
            vsc = [a + b for a, b in zip(vsc, [Vinfx, Vinfy, Vinfz])]
            rsc, vsc = propagate_lagrangian(rsc, vsc, T[i] * DAY2SEC, self.__common_mu)
        cw = (ic2par(rsc, vsc, self.start.mu_central_body)[2] > pi / 2)

        # We now compute the remaining two final impulses
        # Lambert arc to reach seq[1]
        dt = T[-1] * DAY2SEC
        l = lambert_problem(rsc, r_target, dt, self.__common_mu, cw, False)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        DV1 = norm([a - b for a, b in zip(v_beg_l, vsc)])
        DV2 = norm([a - b for a, b in zip(v_end_l, v_target)])

        DV_others = list(x[5::4])
        DV_others.extend([DV1, DV2])

        print("Dvs (m/s): ", DV_others)
        print("Tofs (days): ", T)