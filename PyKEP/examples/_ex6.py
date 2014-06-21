def run_example6(n_seg = 5 ):
    """
    This example demonstrates the optimization of a multiple randezvous mission (low-thrust). 
    Such a mission (including more asteroids) is also called asteroid hopping

    The spacecraft performances, as well as the three asteroids visited, are taken from the GTOC7 problem description.
    """
    from PyGMO import algorithm, population
    from PyKEP.trajopt import mr_lt_nep
    from PyKEP import planet_gtoc7

    algo = algorithm.scipy_slsqp(max_iter = 500, acc=1e-5,screen_output=True)

    prob = mr_lt_nep(
                        t0=[9600.,9700.],
                        seq=[planet_gtoc7(5318), planet_gtoc7(14254), planet_gtoc7(7422), planet_gtoc7(5028)], 
                        n_seg=n_seg, 
                        mass=[800.,2000.], 
                        leg_tof=[100., 365.25], 
                        rest=[30., 365.25], 
                        Tmax=0.3, 
                        Isp=3000., 
                        traj_tof=365.25 * 3., 
                        objective='mass',
                        c_tol = 1e-05
                    )

    pop = population(prob, 1)
    pop = algo.evolve(pop)
    
    solution = pop.champion.x
    if prob.feasibility_x(solution):
        print "FEASIBILE!!!"
        ax = prob.plot(solution)
    else:
        print "INFEASIBLE :("
        ax = None

    return prob,solution, ax



