def run_example6(n_seg=5):
    """
    This example demonstrates the optimization of a multiple rendezvous mission (low-thrust).
    Such a mission (including more asteroids) is also called asteroid hopping

    The spacecraft performances, as well as the three asteroids visited, are taken from the GTOC7 problem description.
    """
    import pygmo as pg
    from pykep.trajopt import mr_lt_nep
    from pykep.planet import gtoc7
    from pykep.examples import add_gradient, algo_factory
    from matplotlib import pyplot as plt

    # If you have no access to snopt7, try slsqp (multiple starts may be
    # necessary)
    algo = algo_factory("snopt7")

    udp = add_gradient(
        mr_lt_nep(
            t0=[9600., 9700.],
            seq=[gtoc7(5318), gtoc7(14254), gtoc7(7422), gtoc7(5028)],
            n_seg=n_seg,
            mass=[800., 2000.],
            leg_tof=[100., 365.25],
            rest=[30., 365.25],
            Tmax=0.3,
            Isp=3000.,
            traj_tof=365.25 * 3.,
            objective='mass',
            c_tol=1e-05
        ),
        with_grad=False)

    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    pop = pg.population(prob, 1)
    pop = algo.evolve(pop)

    solution = pop.champion_x
    if prob.feasibility_x(solution):
        print("FEASIBILE!!!")
        ax = udp.udp_inner.plot(solution)
    else:
        print("INFEASIBLE :(")
        ax = None

    plt.ion()
    plt.show()

if __name__ == "__main__":
    run_example6()
