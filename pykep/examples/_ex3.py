import pykep as pk

class mga_lt_EVMe(object):

    """
    This constructs a pagmo udp that represents a low-thrust transfer between Earth and Mercury with a Venus
    fly-by. The decision vector contains
    [
    t0,
    T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1,
    T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2,
    [throttles1], [throttles2]
    ]
    in the following units: [mjd2000, days, kg, m/s,m/s,m/s, [non-dimensional]]
    """

    def __init__(self, mass=2000, Tmax=0.5, Isp=3500, Vinf_dep=3, Vinf_arr=2, nseg1=5, nseg2=20):
        # We define some data members (we use the double underscore to
        # indicate they are private)
        self.__earth = pk.planet.jpl_lp('earth')
        self.__venus = pk.planet.jpl_lp('venus')
        self.__mercury = pk.planet.jpl_lp('mercury')
        self.__sc = pk.sims_flanagan.spacecraft(mass, Tmax, Isp)
        self.__Vinf_dep = Vinf_dep * 1000
        self.__Vinf_arr = Vinf_arr * 1000
        self.__leg1 = pk.sims_flanagan.leg()
        self.__leg2 = pk.sims_flanagan.leg()
        self.__leg1.set_mu(pk.MU_SUN)
        self.__leg1.set_spacecraft(self.__sc)
        self.__leg2.set_mu(pk.MU_SUN)
        self.__leg2.set_spacecraft(self.__sc)
        self.__nseg1 = nseg1
        self.__nseg2 = nseg2
        self.__mass = mass

    def get_nec(self):
        return 15

    def get_nic(self):
        return self.__nseg1 + self.__nseg2 + 3

    def get_bounds(self):
        lb = [3000, 100, self.__mass / 2] + [-self.__Vinf_dep] * 3 + [-6000] * 3 + [200, self.__mass / 9] + [-6000] * 3 + [-self.__Vinf_arr] * 3 + [-1, -1, -1] * (self.__nseg1 + self.__nseg2)
        ub = [4000, 1000, self.__mass] + [self.__Vinf_dep] * 3 + [6000] * 3 + [2000, self.__mass] + [6000] * 3 + [self.__Vinf_arr] * 3 + [1, 1, 1] * (self.__nseg1 + self.__nseg2)
        return (lb, ub)

    # This is the objective function
    def fitness(self, x):
        from pykep import epoch, AU, EARTH_VELOCITY, fb_con
        from pykep.sims_flanagan import leg, sc_state
        from numpy.linalg import norm
        from math import sqrt, asin, acos

        retval = [-x[10]]

        # Ephemerides
        t_E = epoch(x[0])
        t_V = epoch(x[0] + x[1])
        t_M = epoch(x[0] + x[1] + x[9])
        rE, vE = self.__earth.eph(t_E)
        rV, vV = self.__venus.eph(t_V)
        rM, vM = self.__mercury.eph(t_M)

        # First Leg
        v = [a + b for a, b in zip(vE, x[3:6])]
        x0 = sc_state(rE, v, self.__sc.mass)
        v = [a + b for a, b in zip(vV, x[6:9])]
        xe = sc_state(rV, v, x[2])
        self.__leg1.set(
            t_E, x0, x[-3 * (self.__nseg1 + self.__nseg2):-self.__nseg2 * 3], t_V, xe)

        # Second leg
        v = [a + b for a, b in zip(vV, x[11:14])]
        x0 = sc_state(rV, v, x[2])
        v = [a + b for a, b in zip(vM, x[14:17])]
        xe = sc_state(rM, v, x[10])
        self.__leg2.set(t_E, x0, x[(-3 * self.__nseg2):], t_V, xe)

        # Defining the constraints
        # departure
        v_dep_con = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5] -
                        self.__Vinf_dep * self.__Vinf_dep) / (EARTH_VELOCITY * EARTH_VELOCITY)
        # arrival
        v_arr_con = (x[14] * x[14] + x[15] * x[15] + x[16] * x[16] -
                        self.__Vinf_arr * self.__Vinf_arr) / (EARTH_VELOCITY * EARTH_VELOCITY)
        # fly-by at Venus
        DV_eq, alpha_ineq = fb_con(x[6:9], x[11:14], self.__venus)

        # Assembling the constraints
        constraints = list(self.__leg1.mismatch_constraints() + self.__leg2.mismatch_constraints()) + [DV_eq] + list(
            self.__leg1.throttles_constraints() + self.__leg2.throttles_constraints()) + [v_dep_con] + [v_arr_con] + [alpha_ineq]

        # We then scale all constraints to non-dimensional values
        # leg 1
        constraints[0] /= AU
        constraints[1] /= AU
        constraints[2] /= AU
        constraints[3] /= EARTH_VELOCITY
        constraints[4] /= EARTH_VELOCITY
        constraints[5] /= EARTH_VELOCITY
        constraints[6] /= self.__sc.mass
        # leg 2
        constraints[7] /= AU
        constraints[8] /= AU
        constraints[9] /= AU
        constraints[10] /= EARTH_VELOCITY
        constraints[11] /= EARTH_VELOCITY
        constraints[12] /= EARTH_VELOCITY
        constraints[13] /= self.__sc.mass
        # fly-by at Venus
        constraints[14] /= (EARTH_VELOCITY * EARTH_VELOCITY)

        return retval + constraints

    # And this helps to visualize the trajectory
    def plot(self, x):
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from pykep import epoch, AU
        from pykep.sims_flanagan import sc_state
        from pykep.orbit_plots import plot_planet, plot_sf_leg

        t_E = epoch(x[0])
        t_V = epoch(x[0] + x[1])
        t_M = epoch(x[0] + x[1] + x[9])
        rE, vE = self.__earth.eph(t_E)
        rV, vV = self.__venus.eph(t_V)
        rM, vM = self.__mercury.eph(t_M)

        # First Leg
        v = [a + b for a, b in zip(vE, x[3:6])]
        x0 = sc_state(rE, v, self.__sc.mass)
        v = [a + b for a, b in zip(vV, x[6:9])]
        xe = sc_state(rV, v, x[2])
        self.__leg1.set(
            t_E, x0, x[-3 * (self.__nseg1 + self.__nseg2):-self.__nseg2 * 3], t_V, xe)

        # Second leg
        v = [a + b for a, b in zip(vV, x[11:14])]
        x0 = sc_state(rV, v, x[2])
        v = [a + b for a, b in zip(vM, x[14:17])]
        xe = sc_state(rM, v, x[10])
        self.__leg2.set(t_E, x0, x[(-3 * self.__nseg2):], t_V, xe)

        fig = plt.figure()
        axis = fig.gca(projection='3d')

        # The Sun
        axis.scatter([0], [0], [0], color='y')
        # The legs
        plot_sf_leg(self.__leg1, units=AU, N=10, ax=axis)
        plot_sf_leg(self.__leg2, units=AU, N=10, ax=axis)
        # The planets
        plot_planet(
            self.__earth, t_E, units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)
        plot_planet(
            self.__venus, t_V, units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)
        plot_planet(
            self.__mercury, t_M, units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)
        plt.show()

"""
This example constructs, using pygmo for optimization, an interplanetary low-thrust optimization
problem that can then be solved using one of the available pygmo solvers. The problem is a non-linear constrained
problem that uses the Sims-Flanagan transcription to model the low-thrust trajectory. pykep plotting capabilities
are also demonstrated via the plot method. The interplanetary mission modelled is an LT-MGA Earth-Venus-Mercury mission.

"""
def run_example3():
    import pygmo as pg
    from pykep.examples import add_gradient, algo_factory

    # problem
    udp = add_gradient(mga_lt_EVMe(), with_grad=False)
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # algorithm
    uda = algo_factory("snopt7", False)
    uda2 = pg.mbh(uda, 5, 0.05)
    algo = pg.algorithm(uda2)
    algo.set_verbosity(1)

    # 3 - Population
    pop = pg.population(prob, 1)

    # 4 - Solve the problem (evolve)
    print("Running Monotonic Basin Hopping ....")
    pop = algo.evolve(pop)

    print("Is the solution found a feasible trajectory? " +
            str(prob.feasibility_x(pop.champion_x)))
    udp.udp_inner.plot(pop.champion_x)

if __name__ == "__main__":
    run_example3()
