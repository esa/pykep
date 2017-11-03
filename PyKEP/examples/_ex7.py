def run_example7():
    """
    This example demonstrates the optimization of a transfer from Earth to Mars orbit (no phasing). It makes
    use of an indirect method based on the class trajopt.indirect_or2or which defines an optimization problem
    in the pagmo formalism a user defined problem (UDP) and solves it calling some pagmo algorithm (UDA).

    An algorithm factory is provided to experiment with, for example, ipopt or slsqp (scipy) or, if available
    snopt7
    """
    import PyKEP as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt

    # We define a meta-problem augmenting a generic UDP with a numercial gradient so that algorithms
    # needing the gradient can be used. Note that algorithm such as snopt7 have an internal system to compute gradients
    # in case the problem is fed to it without. Hence the kwarg with_grad is provided to deactivate the gradient.
    class add_gradient:
        def __init__(self, udp, with_grad=False):
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
            # we here use the low precision gradient
            return pg.estimate_gradient(lambda x: self.fitness(x), x, 1e-8)

        def has_gradient(self):
            return self.with_grad

    # We define a utility class that sets the algorithm to be used (in pagmo's terminology the UDA or user defined algorithm).
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

    # ------------------------------------ EXPERIMENT SET-UP -----------------------------------
    # time of flight bounds [days]
    Tlb = 100
    Tub = 700

    # eccentric anomaly bounds
    pi = 3.14159265359
    E0lb = -4 * pi
    E0ub = 4 * pi
    Eflb = -4 * pi
    Efub = 4 * pi

    # Keplerian elements of the starting and target orbit (EARTH-MARS in this case)
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

    # Some pre-computed solutions (obtaining spending some iterations and restarts from radnom initial guesses)
    z_mass_optimal = [397.88267909228767, 2.0719343674552215, 2.7313941033119407, 11.882803539732214, 10.567639551625298, 0.50803671389927796, -11.056641527923768, 12.176151455434058, 4.1269457596809245, 3.4434953247725324]
    z_quadratic_control = [381.32031472240106, 1.0102363292172423, 1.8134352964367946, 19.522442569527868, -6.7894353762521105, -3.3749783899165928, 7.0438655057343054, 19.923912672512174, 0.93896446800741751, 3.5483645070393743]
    z_quadratic_control2 = [459.51623108767666, -1.616057488705803, 0.33049652475302532, -17.735981532357027, -3.2374905349904912, 2.2249621531880934, 2.9550456430212937, -20.226761256676323, -2.9684113654904061, 3.1471248891703905]
    z_quadratic_control3 = [519.45371103815569, 0.39617485433378341, 2.7008977766929818, 7.9136210333255468, -11.03747486077437, -3.0776988186969136, 9.1796869310249747, 6.5013311040515687, -0.2054349910826633, 3.0084671211666865]
    # A random initial guess
    z_random = np.hstack([[np.random.uniform(Tub, Tlb)], 2 * np.random.randn(9)])

    # We use an initial guess within 10% of a known optima, you can experiment what happens with a differet choice
    z = z_quadratic_control + z_quadratic_control * np.random.randn(10) * 0.1
    #z = z_random

    # Problem. We start defining a minimum quadratic control problem (alpha=0) with free time (hamiltonian will be foced to be 0)
    udp = add_gradient(pk.trajopt.indirect_or2or(elem0, elemf, mass, thrust, isp, 1e-12, 1e-12, Tlb, Tub,
                                                 E0lb, E0ub, Eflb, Efub, freetime=True, alpha=0, bound=True, mu=mu), with_grad=True)
    prob = pg.problem(udp)
    prob.c_tol = [1e-7] * 10

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
    arr = udp.udp_inner.leg.get_states(1e-12, 1e-12)
    print("Final mass is: ", arr[-1, 7])

    plt.ion()
    plt.show()

if __name__ == "__main__":
    run_example7()
