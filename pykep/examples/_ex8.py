def run_example8():
    """
    This example demonstrates the direct method (sims-flanagan) on a planet to planet scenario.
    """
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory


    # 1 - Algorithm
    algo = algo_factory("snopt7")

    # 2 - Problem
    udp = add_gradient(pk.trajopt.direct_pl2pl(
        p0="earth",
        pf="mars",
        mass=1000,
        thrust=0.3,
        isp=3000,
        vinf_arr=1e-6,
        vinf_dep=3.5,
        hf=True,
        nseg=20,
        t0=[1100, 1400],
        tof=[200, 290]),
        with_grad=False
    )

    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # 3 - Population
    pop = pg.population(prob, 1)


    # 4 - Solve the problem (evolve)
    pop = algo.evolve(pop)

    # 5 - Inspect the solution
    if prob.feasibility_x(pop.champion_x):
        print("Optimal Found!!")
    else:
        print("No solution found, try again :)")

    udp.udp_inner.pretty(pop.champion_x)

    axis = udp.udp_inner.plot_traj(pop.champion_x)
    plt.title("The trajectory in the heliocentric frame")
    axis = udp.udp_inner.plot_control(pop.champion_x)
    plt.title("The control profile (throttle)")   

    plt.ion()
    plt.show()


if __name__ == "__main__":
    run_example8()
