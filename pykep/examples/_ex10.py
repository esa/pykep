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
        x0 = [-51051524893.335152, -142842795180.97464, 1139935.2553601924, 30488.847061907356, -10612.482697050367, -204.23284335657095, 1000],
        xf = [24753885674.871033, 231247560000.17883, 4236305010.4256544, -23171.900670190855, 4635.6817290400222, 666.44019588506023, 910.48383959441833],
        thrust = 0.3,
        isp = 3000,
        mu = pk.MU_SUN,
        tof=[276.15166075931495, 276.15166075931495],
        freetime=False, 
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
    
    # 5 - Continue the solution to mass optimal
    homotopy_path = [0.5,0.75, 0.9,1]
    for alpha in homotopy_path:
        z =  pop.champion_x
        udp = add_gradient(pk.trajopt.indirect_pt2pt(
            x0 = [-51051524893.335152, -142842795180.97464, 1139935.2553601924, 30488.847061907356, -10612.482697050367, -204.23284335657095, 1000],
            xf = [24753885674.871033, 231247560000.17883, 4236305010.4256544, -23171.900670190855, 4635.6817290400222, 666.44019588506023, 910.48383959441833],
            thrust = 0.3,
            isp = 3000,
            mu = pk.MU_SUN,
            tof=[276.15166075931495, 276.15166075931495],
            freetime=False, 
            alpha=alpha,    # continuation paramter
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
