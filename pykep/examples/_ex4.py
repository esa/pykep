import pykep as pk


class mga_lt_earth_mars_sundmann(object):

    """
    This constructs a pagmo.problem object that represents a low-thrust transfer between Earth and Jupiter.
    The interplantary leg uses Sundman's variable. The decision vector
    contains [t0,T,sf,mf,Vx,Vy,Vz,[throttles]] in the following units: [mjd2000, days, days, kg, m/s,m/s,m/s, [non-dimensional]]
    """

    def __init__(self, mass=1000, Tmax=0.1, Isp=2500, Vinf=3.0, nseg=20):
        # We define some data members (we use the double underscore to
        # indicate they are private)
        self.__earth = pk.planet.jpl_lp('earth')
        self.__mars = pk.planet.jpl_lp('jupiter')
        self.__sc = pk.sims_flanagan.spacecraft(mass, Tmax, Isp)
        self.__Vinf = Vinf * 1000
        # here we construct the trajectory leg in Sundman's variable t =
        # (r/10AU)^1.5 s
        self.__leg = pk.sims_flanagan.leg_s(
            nseg, 1.0 / (100 * pk.AU) ** 1.0, 1.0)
        self.__leg.set_mu(pk.MU_SUN)
        self.__leg.set_spacecraft(self.__sc)
        # This is needed to use the plotting function plot_sf_leg
        self.__leg.high_fidelity = False
        self.__nseg = nseg

    def get_nic(self):
        return self.__nseg + 1

    def get_nec(self):
        return 8

    def get_bounds(self):
        lb = [5000, 2400, 10000, self.__sc.mass / 10, -self.__Vinf, -
              self.__Vinf, -self.__Vinf] + [-1] * 3 * self.__nseg
        ub = [8000, 2500, 150000, self.__sc.mass, self.__Vinf,
              self.__Vinf, self.__Vinf] + [1] * 3 * self.__nseg
        return (lb, ub)

    def fitness(self, x):
        from pykep import epoch, AU, EARTH_VELOCITY, DAY2SEC
        from pykep.sims_flanagan import sc_state
        # This is the objective function
        objfun = [-x[3]]

        # And these are the constraints
        start = epoch(x[0])
        end = epoch(x[0] + x[1])

        r, v = self.__earth.eph(start)
        v = [a + b for a, b in zip(v, x[4:7])]
        x0 = sc_state(r, v, self.__sc.mass)

        r, v = self.__mars.eph(end)
        xe = sc_state(r, v, x[3])
        self.__leg.set(start, x0, x[-3 * self.__nseg:],
                       end, xe, x[2] * DAY2SEC)
        v_inf_con = (x[4] * x[4] + x[5] * x[5] + x[6] * x[6] -
                     self.__Vinf * self.__Vinf) / (EARTH_VELOCITY * EARTH_VELOCITY)
        try:
            constraints = list(self.__leg.mismatch_constraints(
            ) + self.__leg.throttles_constraints()) + [v_inf_con]
        except:
            print(
                "warning: CANNOT EVALUATE constraints .... possible problem in the Taylor integration in the Sundmann variable")
            constraints = (1e14,) * (8 + 1 + self.__nseg + 2)
        # We then scale all constraints to non-dimensional values
        constraints[0] /= AU
        constraints[1] /= AU
        constraints[2] /= AU
        constraints[3] /= EARTH_VELOCITY
        constraints[4] /= EARTH_VELOCITY
        constraints[5] /= EARTH_VELOCITY
        constraints[6] /= self.__sc.mass
        constraints[7] /= 365.25 * DAY2SEC
        return objfun + constraints

    # And this helps to visualize the trajectory
    def plot(self, x):
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from pykep import epoch, AU, DAY2SEC
        from pykep.sims_flanagan import sc_state
        from pykep.orbit_plots import plot_planet, plot_sf_leg

        start = epoch(x[0])
        end = epoch(x[0] + x[1])

        r, v = self.__earth.eph(start)
        v = [a + b for a, b in zip(v, x[4:7])]
        x0 = sc_state(r, v, self.__sc.mass)

        r, v = self.__mars.eph(end)
        xe = sc_state(r, v, x[3])
        self.__leg.set(
            start, x0, x[-3 * self.__nseg:], end, xe, x[2] * DAY2SEC)

        fig = plt.figure()
        axis = fig.gca(projection='3d')
        # The Sun
        axis.scatter([0], [0], [0], color='y')
        # The leg
        plot_sf_leg(self.__leg, units=AU, N=10, ax=axis)
        # The planets
        plot_planet(
            self.__earth, start, units=AU, legend=True, color=(0.8, 0.8, 1), ax=axis)
        plot_planet(
            self.__mars, end, units=AU, legend=True, color=(0.8, 0.8, 1), ax=axis)
        plt.show()


"""
This example demonstrates the use of the interplanetary leg in Sundman's variable to obtain automated mesh optimization
"""


def run_example4():
    import pygmo as pg
    from pykep.examples import add_gradient, algo_factory

    N = 20

    # problem
    udp = add_gradient(mga_lt_earth_mars_sundmann(nseg=N), with_grad=False)
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
    run_example4()
