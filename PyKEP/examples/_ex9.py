def run_example9():
    """
    This example demonstrates the optimization of a transfer from Earth to Mars orbit (with phasing). It makes
    use of a direct method based on the class trajopt.direct_pl2pl which defines an optimization problem,
    in the pagmo formalism a user defined problem (UDP), and solves it calling some pagmo algorithm (UDA).

    An algorithm factory is provided to experiment with, for example, ipopt or slsqp (scipy) or, if available
    snopt7
    """
    import PyKEP as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from ._ex_utilities import add_gradient, algo_factory
 
    # We start defining a minimum quadratic control problem (alpha=0). The initial conditions (x0 and t0) are
    # taken from the solution coming from an indirect method with the final time of flight left free to vary slightly
    udp = add_gradient(pk.trajopt.indirect_pt2pl(
        x0 = [-41243271661.730042, -146069872650.01047, 1169367.2384967851, 31030.968708889388, -8841.9478971260651, -309.95724794416344, 1000],
        t0 = 1251.3387358157067,
        pf = "mars",
        thrust = 0.3,
        isp = 3000,
        mu = pk.MU_SUN,
        tof=[240, 280],
        alpha=1,
        bound = True),
        with_grad=True
    )
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # population
    pop = pg.population(prob)
    z = np.hstack(([np.random.uniform(udp.udp_inner.tof[0], udp.udp_inner.tof[1])], 10 * np.random.randn(7)))
    #z = [270.0, 1.2312118579592024, -6.1390780691032347, 6.0238456490963834, 11.812464867695153, -6.4053652792034494, 0.72804010269853625, 1.854339042824044]
    #z = np.array(z) + np.random.randn(8) * 0.4
    pop.push_back(z)

    # algorithm
    algo = algo_factory("ipopt")
    #algo.set_verbosity(1)

    # evolve (solve the problem)
    pop = algo.evolve(pop)

    # print some information
    print("Feasible?:", prob.feasibility_x(pop.champion_x) )
    udp.udp_inner.pretty(pop.champion_x)
    print(udp.udp_inner.pf.eph(pop.champion_x[0]+udp.udp_inner.t0.mjd2000))
    # plot result
    axis = udp.udp_inner.plot_traj(pop.champion_x)
    plt.title("The trajectory in the heliocentric frame")
    axis = udp.udp_inner.plot_control(pop.champion_x)

    plt.ion()
    plt.show()

    print("F:", list(pop.champion_f))
    print("x:", list(pop.champion_x))



if __name__ == "__main__":
    run_example9()
