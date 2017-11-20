def run_example10():
    """
    This example demonstrates the indirect method (cartesian) on a point to point fixed time scenario.
    The boundary conditions are taken from a run of the indirect method.
    """
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory

    # 1 - Algorithm
    algo = algo_factory("snopt7")

    # 2 - Problem
    udp = add_gradient(pk.trajopt.indirect_pt2pt(
        x0=[44914296854.488266, -145307873786.94177, 1194292.6437741749,
            31252.149474878544, 9873.214642584162, -317.08718075574404, 1000],
        xf=[-30143999066.728119, -218155987244.44385, -3829753551.2279921,
            24917.707565772216, -1235.74045124602, -638.05209482866155, 905.47894037275546],
        thrust=0.1,
        isp=3000,
        mu=pk.MU_SUN,
        tof=[616.77087591237546, 616.77087591237546],
        freetime=False,
        alpha=0,    # quadratic control
        bound=True),
        with_grad=False
    )
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # 3 - Population
    pop = pg.population(prob)
    z = np.hstack(([np.random.uniform(udp.udp_inner.tof[0],
                                      udp.udp_inner.tof[1])], 10 * np.random.randn(7)))
    pop.push_back(z)

    # 4 - Solve the problem (evolve)
    pop = algo.evolve(pop)

    # 5 - Continue the solution to mass optimal
    homotopy_path = [0.5, 0.75, 0.9, 1]
    for alpha in homotopy_path:
        z = pop.champion_x
        print("alpha: ", alpha)
        udp = add_gradient(pk.trajopt.indirect_pt2pt(
            x0=[44914296854.488266, -145307873786.94177, 1194292.6437741749,
                31252.149474878544, 9873.214642584162, -317.08718075574404, 1000],
            xf=[-30143999066.728119, -218155987244.44385, -3829753551.2279921,
                24917.707565772216, -1235.74045124602, -638.05209482866155, 905.47894037275546],
            thrust=0.1,
            isp=3000,
            mu=pk.MU_SUN,
            tof=[616.77087591237546, 616.77087591237546],
            freetime=False,
            alpha=alpha,    # quadratic control
            bound=True),
            with_grad=False
        )
        prob = pg.problem(udp)
        prob.c_tol = [1e-5] * prob.get_nc()

        # 7 - Solve it
        pop = pg.population(prob)
        pop.push_back(z)
        pop = algo.evolve(pop)

    # 8 - Inspect the solution
    print("Feasible?:", prob.feasibility_x(pop.champion_x))

    # plot trajectory
    axis = udp.udp_inner.plot_traj(pop.champion_x, quiver=True, mark="k")
    plt.title("The trajectory in the heliocentric frame")

    # plot control
    udp.udp_inner.plot_control(pop.champion_x)
    plt.title("The control profile (throttle)")

    plt.ion()
    plt.show()

    udp.udp_inner.pretty(pop.champion_x)

    print("\nDecision vector: ", list(pop.champion_x))


if __name__ == "__main__":
    run_example10()
