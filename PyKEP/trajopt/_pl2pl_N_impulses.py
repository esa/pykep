from PyGMO.problem import base as base_problem
from PyKEP.core import epoch, DAY2SEC, planet_ss, lambert_problem, propagate_lagrangian, SEC2DAY, AU, ic2par
from math import pi, cos, sin, log, acos
from scipy.linalg import norm


class pl2pl_N_impulses(base_problem):
    """
    This class is a PyGMO (http://esa.github.io/pygmo/) problem representing a single leg transfer
    between two planets allowing up to a maximum number of impulsive Deep Space Manouvres.

    The decision vector is::

      [t0,T] + [alpha,u,v,V_inf]*(N-2) +[alpha] + ([tf])

    ... in the units: [mjd2000, days] + [nd, nd, m/sec, nd] + [nd] + [mjd2000]

    Each leg time-of-flight can be decoded as follows, T_n = T log(alpha_n) / \sum_i(log(alpha_i))

    .. note::

       The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.
    """

    def __init__(self,
                 start=planet_ss('earth'),
                 target=planet_ss('venus'),
                 N_max=3,
                 tof=[20., 400.],
                 vinf=[0., 4.],
                 phase_free=True,
                 multi_objective=False,
                 t0=None
                 ):
        """
        prob = PyKEP.trajopt.pl2pl_N_impulses(start=planet_ss('earth'), target=planet_ss('venus'), N_max=3, tof=[20., 400.], vinf=[0., 4.], phase_free=True, multi_objective=False, t0=None)

        - start:            a PyKEP planet defining the starting orbit
        - target:           a PyKEP planet defining the target orbit
        - N_max:            maximum number of impulses
        - tof:              a list containing the box bounds [lower,upper] for the time of flight (days)
        - vinf:             a list containing the box bounds [lower,upper] for each DV magnitude (km/sec)
        - phase_free:       when True, no randezvous condition is enforced and the final orbit will be reached at an optimal true anomaly
        - multi_objective:  when True, a multi-objective problem is constructed with DV and time of flight as objectives
        - t0:               launch window defined as a list of two epochs [epoch,epoch]
        """

        # Sanity checks
        # 1) all planets need to have the same mu_central_body
        if (start.mu_central_body != target.mu_central_body):
            raise ValueError('Starting and ending PyKEP.planet must have the same mu_central_body')
        # 2) Number of impulses must be at least 2
        if N_max < 2:
            raise ValueError('Number of impulses N is less than 2')
        # 3) If phase_free is True, t0 does not make sense
        if (t0 is None and not phase_free):
            t0 = [epoch(0), epoch(1000)]
        if (t0 is not None and phase_free):
            raise ValueError('When phase_free is True no t0 can be specified')

        # We compute the PyGMO problem dimensions
        dim = 2 + 4 * (N_max - 2) + 1 + phase_free
        obj_dim = multi_objective + 1
        # First we call the constructor for the base PyGMO problem
        # As our problem is n dimensional, box-bounded (may be multi-objective), we write
        # (dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
        super(pl2pl_N_impulses, self).__init__(dim, 0, obj_dim, 0, 0, 0)

        # We then define all class data members
        self.start = start
        self.target = target
        self.N_max = N_max
        self.phase_free = phase_free
        self.multi_objective = multi_objective

        self.__common_mu = start.mu_central_body

        # And we compute the bounds
        if phase_free:
            lb = [start.ref_epoch.mjd2000, tof[0]] + [0.0, 0.0, 0.0, vinf[0] * 1000] * (N_max - 2) + [0.0] + [target.ref_epoch.mjd2000]
            ub = [start.ref_epoch.mjd2000 + 2 * start.period * SEC2DAY, tof[1]] + [1.0, 1.0, 1.0, vinf[1] * 1000] * (N_max - 2) + [1.0] + [target.ref_epoch.mjd2000 + 2 * target.period * SEC2DAY]
        else:
            lb = [t0[0].mjd2000, tof[0]] + [0.0, 0.0, 0.0, vinf[0] * 1000] * (N_max - 2) + [0.0]
            ub = [t0[1].mjd2000, tof[1]] + [1.0, 1.0, 1.0, vinf[1] * 1000] * (N_max - 2) + [1.0]

        # And we set them
        self.set_bounds(lb, ub)

    # Objective function
    def _objfun_impl(self, x):
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
            rsc, vsc = propagate_lagrangian(
                rsc, vsc, T[i] * DAY2SEC, self.__common_mu)
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
        if self.f_dimension == 1:
            return (DV1 + DV2 + DV_others,)
        else:
            return (DV1 + DV2 + DV_others, x[1])

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
        from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler

        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.gca(projection='3d')
        else:
            axis = ax

        axis.scatter(0, 0, 0, color='y')

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
        plot_planet(self.start, t0=epoch(x[0]), color=(0.8, 0.6, 0.8), legend=True, units = AU, ax=axis)
        plot_planet(self.target, t0=epoch(x[0] + x[1]), color=(0.8, 0.6, 0.8), legend=True, units = AU, ax=axis)

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
            plot_kepler(rsc, vsc, T[
                        i] * DAY2SEC, self.__common_mu, N=200, color='b', legend=False, units=AU, ax=axis)
            rsc, vsc = propagate_lagrangian(
                rsc, vsc, T[i] * DAY2SEC, self.__common_mu)

        cw = (ic2par(rsc, vsc, self.start.mu_central_body)[2] > pi / 2)
        # We now compute the remaining two final impulses
        # Lambert arc to reach seq[1]
        dt = T[-1] * DAY2SEC
        l = lambert_problem(rsc, r_target, dt, self.__common_mu, cw, False)
        plot_lambert(
            l, sol=0, color='r', legend=False, units=AU, ax=axis, N=200)
        plt.show()
        return axis
