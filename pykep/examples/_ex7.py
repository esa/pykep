def run_example7(solver="snopt7"):
    """
    This example demonstrates the indirect method (cartesian) on a orbit to orbit scenario.
    The orbits are those of Earth and Mars.
    """
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory

    # Some pre-computed solutions (obtaining spending some iterations and
    # restarts from random initial guesses)
    z_mass_optimal = [397.88267909228767, 2.0719343674552215, 2.7313941033119407, 11.882803539732214, 10.567639551625298,
                      0.50803671389927796, -11.056641527923768, 12.176151455434058, 4.1269457596809245, 3.4434953247725324]
    z_quadratic_control = [381.32031472240106, 1.0102363292172423, 1.8134352964367946, 19.522442569527868, -
                           6.7894353762521105, -3.3749783899165928, 7.0438655057343054, 19.923912672512174, 0.93896446800741751, 3.5483645070393743]
    z_quadratic_control2 = [459.51623108767666, -1.616057488705803, 0.33049652475302532, -17.735981532357027, -
                            3.2374905349904912, 2.2249621531880934, 2.9550456430212937, -20.226761256676323, -2.9684113654904061, 3.1471248891703905]
    z_quadratic_control3 = [519.45371103815569, 0.39617485433378341, 2.7008977766929818, 7.9136210333255468, -
                            11.03747486077437, -3.0776988186969136, 9.1796869310249747, 6.5013311040515687, -0.2054349910826633, 3.0084671211666865]
    # A random initial solution
    z_random = np.hstack(
        [[np.random.uniform(100, 700)], 2 * np.random.randn(9)])

    # We use an initial guess within 10% of a known optima, you can experiment what happens
    # with a different choice
    z = z_quadratic_control + z_quadratic_control * np.random.randn(10) * 0.1
    #z = z_random

    # 1 - Algorithm
    algo = algo_factory(solver)

    # 2 - Problem. We define a minimum quadratic control problem (alpha=0) with free time
    # (hamiltonian will be forced to be 0). We provide the option for estimating the gradient numerically for
    # algorithms that require it.
    udp = add_gradient(pk.trajopt.indirect_or2or(
        elem0=[149598261129.93335, 0.016711230601231957,
               2.640492490927786e-07, 3.141592653589793, 4.938194050401601, 0],
        elemf=[227943822376.03537, 0.09339409892101332,
               0.032283207367640024, 0.8649771996521327, 5.000312830124232, 0],
        mass=1000,
        thrust=0.3,
        isp=2500,
        atol=1e-12,
        rtol=1e-12,
        tof=[100, 700],
        freetime=True,
        alpha=0,
        bound=True,
        mu=pk.MU_SUN),
        with_grad=True)

    prob = pg.problem(udp)
    prob.c_tol = [1e-7] * 10

    # 3 - Population (i.e. initial guess)
    pop = pg.population(prob)
    pop.push_back(z)

    # 4 - Solve the problem (evolve)
    pop = algo.evolve(pop)

    # 5 - Inspect the solution
    if prob.feasibility_x(pop.champion_x):
        print("Optimal Found!!")
        # We call the fitness to set the leg
        udp.udp_inner.fitness(pop.champion_x)
        arr = udp.udp_inner.leg.get_states(1e-12, 1e-12)
        print("Final mass is: ", arr[-1, 7])
    else:
        print("No solution found, try again :)")
    # plot trajectory
    axis = udp.udp_inner.plot_traj(
        pop.champion_x, quiver=True, mark="k", length=1)
    plt.title("The trajectory in the heliocentric frame")

    # plot control
    udp.udp_inner.plot_control(pop.champion_x)
    plt.title("The control profile (throttle)")
    plt.ion()
    plt.show()

    # Show the trajectory data
    udp.udp_inner.pretty(pop.champion_x)

if __name__ == "__main__":
    run_example7()
