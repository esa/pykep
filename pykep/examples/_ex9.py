def run_example9():
    """
    This example demonstrates the indirect method (cartesian) on a point to planet variable time scenario.
    The starting conditions are taken from a run of the indirect method.
    """
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory

    # 1 - Algorithm
    algo = algo_factory("snopt7")

    # 2 - Problem
    udp = add_gradient(pk.trajopt.indirect_pt2pl(
        x0 = [-41243271661.730042, -146069872650.01047, 1169367.2384967851, 31030.968708889388, -8841.9478971260651, -309.95724794416344, 1000],
        t0 = 1251.3387358157067,
        pf = "mars",
        thrust = 0.3,
        isp = 3000,
        mu = pk.MU_SUN,
        tof=[240, 280],
        alpha=0,    # quadratic control
        bound = True),
        with_grad=False
    )
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # 3 - Population
    pop = pg.population(prob)
    z = np.hstack(([np.random.uniform(udp.udp_inner.tof[0], udp.udp_inner.tof[1])], 10 * np.random.randn(7)))
    pop.push_back(z)

    # 4 - Solve the problem (evolve)
    pop = algo.evolve(pop)

    # 5 - Extract the solution (quadratic control)
    z =  pop.champion_x

    # 6 - Define a mass optimal problem
    udp = add_gradient(pk.trajopt.indirect_pt2pl(
        x0 = [-41243271661.730042, -146069872650.01047, 1169367.2384967851, 31030.968708889388, -8841.9478971260651, -309.95724794416344, 1000],
        t0 = 1251.3387358157067,
        pf = "mars",
        thrust = 0.3,
        isp = 3000,
        mu = pk.MU_SUN,
        tof=[240, 280],
        alpha=1,    # mass optimal
        bound = True),
        with_grad=False
    )
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # 7 - Solve it
    pop = pg.population(prob)
    pop.push_back(z)
    pop = algo.evolve(pop)

    # 8 - Inspect the solution
    print("Feasible?:", prob.feasibility_x(pop.champion_x) )

    # plot trajectory
    axis = udp.udp_inner.plot_traj(pop.champion_x, quiver=True, mark='k')
    plt.title("The trajectory in the heliocentric frame")

    # plot control
    udp.udp_inner.plot_control(pop.champion_x)
    plt.title("The control profile (throttle)")   

    plt.ion()
    plt.show()

    udp.udp_inner.pretty(pop.champion_x)

    print("\nDecision vector: ", list(pop.champion_x))



if __name__ == "__main__":
    run_example9()
