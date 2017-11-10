def run_example8():
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

    # Problem. We start defining a minimum quadratic control problem (alpha=0)
    # with free time (hamiltonian will be foced to be 0)
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

    # population
    pop = pg.population(prob, 1)

    # algorithm
    algo = algo_factory("snopt7")

    # evolve
    pop = algo.evolve(pop)

    if prob.feasibility_x(pop.champion_x):
        print("Optimal Found!!")
        udp.udp_inner.pretty(pop.champion_x)
    else:
        print("No solution found, try again :)")

    # plot trajectory
    axis = udp.udp_inner.plot_traj(pop.champion_x)
    plt.title("The trajectory in the heliocentric frame")

    axis = udp.udp_inner.plot_control(pop.champion_x)

    plt.ion()
    plt.show()


if __name__ == "__main__":
    run_example8()
