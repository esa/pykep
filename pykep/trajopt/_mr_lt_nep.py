from PyGMO.problem._base import base as base_problem
from pykep.planet import gtoc7


class mr_lt_nep(base_problem):

    """
    This class represents, as a global optimization problem (linearly constrained,
    high diemensional), a Multiple Randezvous trajectory of a low-thrust spacecraft equipped
    with a nuclear electric propulsion engine.

    - Izzo, D. et al., GTOC7 - Solution Returned by the ACT/ISAS team

    The decision vector (chromosome) is:
    [t1, tof1, rest1, m_f1] + [throttles1] +
    [t2, tof2, rest2, m_f2] + [throttles2] + ....
    .... + [total_tof]

    where the units are [mjd2000, days, days, kg] + [n/a] + .... + [days]

    .. note::

       The resulting problem is non linearly constrained. The resulting trajectory is not time-bounded.
    """

    def __init__(
            self,
            seq=[gtoc7(3413), gtoc7(234), gtoc7(11432)],
            n_seg=5,
            t0=[13000, 13200],
            leg_tof=[1, 365.25 * 3],
            rest=[30., 365.25],
            mass=[800, 2000],
            Tmax=0.3,
            Isp=3000.,
            traj_tof=365.25 * 6,
            objective='mass',
            c_tol=1e-05
    ):
        """
        prob = mr_lt_nep(seq=[pykep.gtoc7(3413),pykep.gtoc7(234), pykep.gtoc7(11432)], n_seg=5, t0=[13000, 13200],
                    leg_tof=[1, 365.25 * 3], rest=[30., 365.25], mass=[800, 2000], Tmax=0.3,
                    Isp=3000., traj_tof=365.25 * 6, objective='mass', c_tol=1e-05)

        * seq: list of pykep.planet defining the encounter sequence for the trajectoty (including the initial planet)
        * n_seg: list of integers containing the number of segments to be used for each leg (len(n_seg) = len(seq)-1)
        * t0: list of two pykep epochs defining the launch window
        * leg_tof: list of two floats defining the minimum and maximum time of each leg (days)
        * rest: list of two floats defining the minimum and maximum time the spacecraft can rest at one planet (days)
        * mass: list of two floats defining the minimum final spacecraft mass and the starting spacecraft mass (kg)
        * Tmax: maximum thrust (N)
        * Isp: engine specific impulse (sec)
        * traj_tof maximum total mission duration (days)
        * c_tol: tolerance on the constraints
        """
        import pykep
        # Number of legs
        n = len(seq) - 1
        # Problem dimension
        dim = (4 + n_seg * 3) * n + 1
        # Number of equality constraints
        dim_eq = 7 * n
        # Number of Inequality constraints
        dim_ineq = n * n_seg + n

        # We call the base problem constaructor
        super(mr_lt_nep, self).__init__(dim, 0, 1,
                                        dim_eq + dim_ineq,  # constraint dimension
                                        dim_ineq,  # inequality constraints
                                        1e-5)  # constrain tolerance

        # We define data members
        self.__seq = seq
        self.__num_legs = n
        self.__nseg = n_seg
        self.__dim_leg = 4 + n_seg * 3
        self.__start_mass = mass[1]
        self.__max_total_time = traj_tof

        # We create n distinct legs objects
        self.__legs = []
        for i in range(n):
            self.__legs.append(pykep.sims_flanagan.leg())
        for leg in self.__legs:
            leg.high_fidelity = True
            leg.set_mu(pykep.MU_SUN)

        if objective not in ['mass', 'time']:
            raise ValueError("Error in defining the objective. Was it one of mass or time?")
        self.__objective = objective

        # We set the ptoblem box-bounds
        # set leg bounds
        lb_leg = [t0[0], leg_tof[0], rest[0], mass[0]] + [-1] * n_seg * 3
        ub_leg = [t0[1] + traj_tof * n, leg_tof[1], rest[1], mass[1]] + [1] * n_seg * 3

        # set n leg bounds
        lb = lb_leg * n
        ub = [t0[1], leg_tof[1], rest[1], mass[1]] + [1] * n_seg * 3 + ub_leg * (n - 1)

        # set total time bounds
        lb += [1.]
        ub += [self.__max_total_time]

        self.set_bounds(lb, ub)

    def _objfun_impl(self, x):
        if self.__objective == 'mass':
            final_mass = x[-1 - self.__dim_leg + 3]
            return (-final_mass, )
        elif self.__objective == 'time':
            tof = x[-1 - self.__dim_leg] - x[0]
        return (tof, )

    def _compute_constraints_impl(self, x_full):
        import pykep
        sc_mass = self.__start_mass
        eqs = []
        ineqs = []

        for i in range(self.__num_legs):
            x = x_full[i * self.__dim_leg:(i + 1) * self.__dim_leg]

            start = pykep.epoch(x[0])
            end = pykep.epoch(x[0] + x[1])

            # Computing starting spaceraft state
            r, v = self.__seq[i].eph(start)
            x0 = pykep.sims_flanagan.sc_state(r, v, sc_mass)

            # Computing ending spaceraft state
            r, v = self.__seq[i + 1].eph(end)
            xe = pykep.sims_flanagan.sc_state(r, v, x[3])

            # Building the SF leg
            self.__legs[i].set_spacecraft(pykep.sims_flanagan.spacecraft(sc_mass, .3, 3000.))
            self.__legs[i].set(start, x0, x[-3 * self.__nseg:], end, xe)

            # Setting all constraints
            eqs.extend(self.__legs[i].mismatch_constraints())
            ineqs.extend(self.__legs[i].throttles_constraints())

            eqs[-7] /= pykep.AU
            eqs[-6] /= pykep.AU
            eqs[-5] /= pykep.AU
            eqs[-4] /= pykep.EARTH_VELOCITY
            eqs[-3] /= pykep.EARTH_VELOCITY
            eqs[-2] /= pykep.EARTH_VELOCITY
            eqs[-1] /= self.__start_mass

            sc_mass = x[3]  # update mass to final mass of leg

            if i < self.__num_legs - 1:
                x_next = x_full[(i + 1) * self.__dim_leg:(i + 2) * self.__dim_leg]
                time_ineq = x[0] + x[1] + x[2] - x_next[0]
                ineqs.append(time_ineq / 365.25)
            else:
                final_time_ineq = x[0] + x[1] + x[2] - x_full[0] - x_full[-1]  # <- total time
                ineqs.append(final_time_ineq / 365.25)

        retval = eqs + ineqs
        return retval

    def resting_times(self, x):
        return list(x[2::self.__dim_leg])

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

        # Computing the legs
        self._compute_constraints_impl(x)

        # Plotting the legs
        for leg in self.__legs:
            plot_sf_leg(leg, units=AU, N=10, ax=axis)

        # Plotting the pykep.planet both at departure and arrival dates
        for i in range(self.__num_legs):
            idx = i * self.__dim_leg
            plot_planet(self.__seq[i], epoch(x[idx]), units=AU, legend=True, color=(0.7, 0.7, 1), s=30, ax=axis)
            plot_planet(self.__seq[i + 1], epoch(x[idx] + x[idx + 1]), units=AU, legend=False, color=(0.7, 0.7, 1), s=30, ax=axis)
        plt.show()
        return axis
