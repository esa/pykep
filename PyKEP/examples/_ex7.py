def run_example7():
    import PyKEP as pk
    import pygmo as pg
    import pygmo_plugins_nonfree as pg7
    import numpy as np

    # planets [names]
    p0_name, pf_name = "earth", "mars"

    # departure time bounds [mjd2000]
    t0lb = 0
    t0ub = 0

    # time of flight bounds [days]
    Tlb = 200
    Tub = 500

    # eccentric anomoly bounds
    pi = 3.14159265359
    E0lb = -4*pi
    E0ub = 4*pi
    Eflb = -4*pi
    Efub = 4*pi

    # planets
    p0 = pk.planet.jpl_lp(p0_name)
    pf = pk.planet.jpl_lp(pf_name)

    # Keplerian elements
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

    # spacecraft parametres
    mass = 1000
    thrust = 0.3
    isp = 2500

    # gravitational parameter
    mu = pk.MU_SUN

    # high_fidelity transcription
    hf = False

    # guess
    z = np.array(
        [  3.53926540e+02  , 3.06963808e-01  , 1.00956769e+00  , 8.24089827e+00,
  -1.36968187e+01 , -3.32237224e+00 ,  1.36158790e+01,   7.94030531e+00,
  -2.51714050e+00 ,  3.42226543e+00]
)

    #z = z + z * np.random.randn(len(z)) * 0.1

    # indirect orbit to orbit problem
    print("Solving indirect transcription orbit to orbit problem with quadratic control.")
    udp = pk.trajopt.indirect_or2or(elem0, elemf, mass, thrust, isp, 1e-12, 1e-12, Tlb, Tub, E0lb, E0ub, Eflb, Efub,  freemass=True, freetime=True, alpha=1, bound=True, mu = mu)
    prob = pg.problem(udp)
    prob.c_tol = [1e-6] * 8

    # population
    pop = pg.population(prob)
    pop.push_back(z)

    # algorithm
    #uda = pg.nlopt("auglag")
    #algo = pg.algorithm(uda)
    #algo.extract(pg.nlopt).local_optimizer = pg.nlopt('var2')
    uda = pg7.snopt7(True, "/usr/local/lib/libsnopt7_c.so")
    uda.set_integer_option("Major iterations limit", 4000)
    uda.set_integer_option("Iterations limit", 40000)
    uda.set_numeric_option("Major optimality tolerance", 1e-1)
    uda.set_numeric_option("Major feasibility tolerance", 1e-7)
    algo = pg.algorithm(uda)
    # algo.set_verbosity(1)

    # evolve
    #uda.set_numeric_option("Major feasibility tolerance", 1e-8)

    pop = algo.evolve(pop)

    print(pop.champion_x)

    # plot trajectory
    udp.plot_traj(pop.champion_x)

    # plot control
    udp.plot_control(pop.champion_x)


if __name__ == "__main__":
    run_example7()
