import PyKEP as pk
import pygmo as pg
import pygmo_plugins_nonfree as pg7
import numpy as np

# algorithm
uda = pg7.snopt7(True, "/usr/lib/libsnopt7_c.so")
uda.set_integer_option("Major iterations limit", 4000)
uda.set_integer_option("Iterations limit", 40000)
uda.set_numeric_option("Major optimality tolerance", 1e-2)
uda.set_numeric_option("Major feasibility tolerance", 1e-8)
algo = pg.algorithm(uda)


def run_example7():

    # planets [names]
    p0_name, pf_name = "earth", "mars"

    # departure time bounds [mjd2000]
    t0lb = 0
    t0ub = 0

    # time of flight bounds [days]
    Tlb = 200
    Tub = 1000

    # mean anomoly bounds
    pi = 3.14159265359
    M0lb = 0
    M0ub = 4*pi
    Mflb = 0
    Mfub = 4*pi

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

    # number of segements
    nseg = 50

    # gravitational parametre
    mu = pk.MU_SUN

    # high_fidelity transcription
    hf = True

    # direct orbit to orbit problem
    print("First solving direct transcription orbit to orbit problem.")
    udp = pk.trajopt.direct_or2or(
        elem0, elemf, mass, thrust, isp, nseg, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, mu, hf)
    prob = pg.problem(udp)

    # popoulation
    pop = pg.population(prob, 1)

    # optimise
    algo = pg.algorithm(uda)
    pop = algo.evolve(pop)
    print(pop.champion_x)

    # plot trajectory
    udp.plot_traj(pop.champion_x)

    # plot control
    udp.plot_control(pop.champion_x)

    # get useful information
    T = pop.champion_x[0]
    mf = pop.champion_x[1]
    M0 = pop.champion_x[2]
    Mf = pop.champion_x[3]

    # indirect orbit to orbit problem
    print("Secondly solving indirect transcription orbit to orbit problem with quadratic control.")
    udp = pk.trajopt.indirect_or2or(elem0, elemf, mass, thrust, isp, 1e-10, 1e-10, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, True, True, 1, True, mu)
    prob = pg.problem(udp)

    # guess
    z = np.hstack(([T, M0, Mf, np.random.randn(7)]))

    # population
    pop = pg.population(prob)
    pop.push_back(z)

    # evolve
    algo = pg.algorithm(uda)
    pop = algo.evolve(pop)

    # plot trajectory
    udp.plot_traj(pop.champion_x)

    # plot control
    udp.plot_control(pop.champion_x)


if __name__ == "__main__":
    run_example7()
