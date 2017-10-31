def run_example7():
    import PyKEP as pk
    import pygmo as pg
    import numpy as np

    # We define a meta-problem augmenting a genric UDP with a numercial  gradient. 
    class add_gradient:
        def __init__(self, udp, with_grad = False):
            self.udp_inner = udp
            self.prob = pg.problem(udp)
            self.with_grad = with_grad
        def fitness(self, x):
            return self.prob.fitness(x)
        def get_bounds(self):
            return self.prob.get_bounds()
        def get_nec(self):
            return self.prob.get_nec()
        def get_nic(self):
            return self.prob.get_nic()
        def get_nobj(self):
            return self.prob.get_nobj()
        def gradient(self, x):
            return pg.estimate_gradient(lambda x: self.fitness(x), x, 1e-8) # we here use the low precision gradient
        def has_gradient(self):
            return self.with_grad

    # We define a utility class that sets the uda.
    def algo_factory(name):
        if name is "slsqp":
            uda = pg.nlopt('slsqp')
            uda.xtol_rel = 1e-6
            uda.ftol_rel = 1e-10
            algo = pg.algorithm(uda)
            algo.set_verbosity(1)
            return algo
        elif name is "ipopt":
            uda = pg.ipopt()
            uda.set_integer_option("print_level", 5)
            uda.set_integer_option("acceptable_iter", 4)

            uda.set_numeric_option("tol", 1e-8)
            uda.set_numeric_option("dual_inf_tol", 1e-8)
            uda.set_numeric_option("constr_viol_tol", 1e-8)
            uda.set_numeric_option("compl_inf_tol", 1e-8)

            uda.set_numeric_option("acceptable_tol", 1e-3)
            uda.set_numeric_option("acceptable_dual_inf_tol", 1e-2)
            uda.set_numeric_option("acceptable_constr_viol_tol", 1e-6)
            uda.set_numeric_option("acceptable_compl_inf_tol", 1e-6)

            algo = pg.algorithm(uda)
            return algo
        elif name is "snopt7":
            import pygmo_plugins_nonfree as pg7
            uda = pg7.snopt7(True, "/usr/local/lib/libsnopt7_c.so")
            uda.set_integer_option("Major iterations limit", 2000)
            uda.set_integer_option("Iterations limit", 2000)
            uda.set_numeric_option("Major optimality tolerance", 1e-4)
            uda.set_numeric_option("Major feasibility tolerance", 1e-6)

            algo = pg.algorithm(uda)
            return algo

    # ------------------------- EXPERIMENT SET-UP -----------------------------
    # time of flight bounds [days]
    Tlb = 350
    Tub = Tlb

    # eccentric anomaly bounds
    pi = 3.14159265359
    E0lb = -4 * pi
    E0ub = 4 * pi
    Eflb = -4 * pi
    Efub = 4 * pi

    # Keplerian elements of the strating and target orbit (EARTH-MARS in this case)
    elem0 = [
        149598261129.93335,
        0.016711230601231957,
        2.640492490927786e-07,
        3.141592653589793,
        4.938194050401601,
        0
    ]
    elemf = [
        227943822376.03537,
        0.09339409892101332,
        0.032283207367640024,
        0.8649771996521327,
        5.000312830124232,
        0
    ]

    # spacecraft parameters
    mass = 1000
    thrust = 0.3
    isp = 2500

    # gravitational parameter
    mu = pk.MU_SUN

    # guess
    z = [350.0, 1.2051321359957445, 1.6409933927327405, 23.812518969857098, -2.4681125404961812, -3.6779605780216844, 4.5442844912881384, 25.436534571988112, 1.4661964586914011, 3.9279048602759508]
    z = z + np.random.randn(10) * 0.01
    z[0] = Tlb

    # problem
    udp = add_gradient(pk.trajopt.indirect_or2or(elem0, elemf, mass, thrust, isp, 1e-12, 1e-12, Tlb, Tub,
                                    E0lb, E0ub, Eflb, Efub, freetime=False, alpha=0, bound=False, mu=mu), with_grad = True)
    prob = pg.problem(udp)
    prob.c_tol = [1e-7] * 7

    # population
    pop = pg.population(prob)
    pop.push_back(z)

    # algorithm
    algo = algo_factory("ipopt")

    # evolve
    pop = algo.evolve(pop)

    print(list(pop.champion_x))

    # plot trajectory
    udp.udp_inner.plot_traj(pop.champion_x)

    # plot control
    udp.udp_inner.plot_control(pop.champion_x)

    # 
    udp.fitness(pop.champion_x)
    arr = udp.udp_inner.leg.get_states(1e-12,1e-12)
    print("Final mass is: ", arr[-1,7])

if __name__ == "__main__":
    run_example7()
