import pykep as pk


class mr_lt_nep:
    """
    This class represents, as a global optimization problem (linearly constrained,
    high dimensional), a Multiple Randezvous trajectory of a low-thrust spacecraft equipped
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
            seq=[pk.planet.gtoc7(3413), pk.planet.gtoc7(
                234), pk.planet.gtoc7(11432)],
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

        # Number of legs
        n = len(seq) - 1
        # Problem dimension
        dim = (4 + n_seg * 3) * n + 1
        # Number of equality constraints
        self.__dim_eq = 7 * n
        # Number of Inequality constraints
        self.__dim_ineq = n * n_seg + n

        # We define data members
        self.__seq = seq
        self.__num_legs = n
        self.__nseg = n_seg
        self.__dim_leg = 4 + n_seg * 3
        self.__start_mass = mass[1]
        self.__max_total_time = traj_tof
        self.__t0 = t0
        self.__leg_tof = leg_tof
        self.__rest = rest
        self.__mass = mass

        # We create n distinct legs objects
        self.__legs = []
        for i in range(n):
            self.__legs.append(pk.sims_flanagan.leg())
        for leg in self.__legs:
            leg.high_fidelity = True
            leg.set_mu(pk.MU_SUN)

        if objective not in ['mass', 'time']:
            raise ValueError(
                "Error in defining the objective. Was it one of mass or time?")
        self.__objective = objective

    def fitness(self, x_full):
        retval = []
        # 1 - obj fun
        if self.__objective == 'mass':
            retval.append(x_full[-1 - self.__dim_leg + 3])
        elif self.__objective == 'time':
            retval.append(x_full[-1 - self.__dim_leg] - x_full[0])

        sc_mass = self.__start_mass
        eqs = []
        ineqs = []

        # 2 - constraints
        for i in range(self.__num_legs):
            x = x_full[i * self.__dim_leg:(i + 1) * self.__dim_leg]

            start = pk.epoch(x[0])
            end = pk.epoch(x[0] + x[1])

            # Computing starting spaceraft state
            r, v = self.__seq[i].eph(start)
            x0 = pk.sims_flanagan.sc_state(r, v, sc_mass)

            # Computing ending spaceraft state
            r, v = self.__seq[i + 1].eph(end)
            xe = pk.sims_flanagan.sc_state(r, v, x[3])

            # Building the SF leg
            self.__legs[i].set_spacecraft(
                pk.sims_flanagan.spacecraft(sc_mass, .3, 3000.))
            self.__legs[i].set(start, x0, x[-3 * self.__nseg:], end, xe)

            # Setting all constraints
            eqs.extend(self.__legs[i].mismatch_constraints())
            ineqs.extend(self.__legs[i].throttles_constraints())

            eqs[-7] /= pk.AU
            eqs[-6] /= pk.AU
            eqs[-5] /= pk.AU
            eqs[-4] /= pk.EARTH_VELOCITY
            eqs[-3] /= pk.EARTH_VELOCITY
            eqs[-2] /= pk.EARTH_VELOCITY
            eqs[-1] /= self.__start_mass

            sc_mass = x[3]  # update mass to final mass of leg

            if i < self.__num_legs - 1:
                x_next = x_full[
                    (i + 1) * self.__dim_leg:(i + 2) * self.__dim_leg]
                time_ineq = x[0] + x[1] + x[2] - x_next[0]
                ineqs.append(time_ineq / 365.25)
            else:
                final_time_ineq = x[0] + x[1] + x[2] - \
                    x_full[0] - x_full[-1]  # <- total time
                ineqs.append(final_time_ineq / 365.25)

        retval = retval + eqs + ineqs
        return retval

    def get_bounds(self):
        t0 = self.__t0
        leg_tof = self.__leg_tof
        rest = self.__rest
        mass = self.__mass
        nseg = self.__nseg
        traj_tof = self.__max_total_time
        n = self.__num_legs
        # We set the ptoblem box-bounds
        # set leg bounds
        lb_leg = [t0[0], leg_tof[0], rest[0], mass[0]] + [-1] * nseg * 3
        ub_leg = [t0[1] + traj_tof * n, leg_tof[1],
                  rest[1], mass[1]] + [1] * nseg * 3

        # set n leg bounds
        lb = lb_leg * n
        ub = [t0[1], leg_tof[1], rest[1], mass[1]] + \
            [1] * nseg * 3 + ub_leg * (n - 1)

        # set total time bounds
        lb += [1.]
        ub += [self.__max_total_time]

        return (lb, ub)

    def get_nic(self):
        return self.__dim_ineq

    def get_nec(self):
        return self.__dim_eq

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

        # Plotting the pykep.planet both at departure and arrival dates
        for i in range(self.__num_legs):
            idx = i * self.__dim_leg
            plot_planet(self.__seq[i], epoch(x[idx]), units=AU, legend=True, color=(
                0.7, 0.7, 0.7), s=30, axes=axis)
            plot_planet(self.__seq[i + 1], epoch(x[idx] + x[idx + 1]),
                        units=AU, legend=False, color=(0.7, 0.7, 0.7), s=30, axes=axis)

        # Computing the legs
        self.fitness(x)

        # Plotting the legs
        for leg in self.__legs:
            plot_sf_leg(leg, units=AU, N=10, axes=axis, legend=False)

        return axis
