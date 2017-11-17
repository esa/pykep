def run_example1(impulses=4):
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory

    # problem
    udp = add_gradient(pk.trajopt.pl2pl_N_impulses(
        start=pk.planet.jpl_lp('earth'),
        target=pk.planet.jpl_lp('venus'),
        N_max=impulses,
        tof=[100., 1000.],
        vinf=[0., 4],
        phase_free=False,
        multi_objective=False,
        t0=[pk.epoch(0), pk.epoch(1000)]), with_grad=False)
    prob = pg.problem(udp)

    # algorithm
    uda = pg.cmaes(gen=1000)
    algo = pg.algorithm(uda)
    algo.set_verbosity(10)

    # population
    pop = pg.population(prob, 20)

    # solve the problem
    pop = algo.evolve(pop)

    # inspect the solution
    udp.udp_inner.plot_trajectory(pop.champion_x)
    plt.ion()
    plt.show()

    udp.udp_inner.pretty(pop.champion_x)

if __name__ == "__main__":
    run_example1()
