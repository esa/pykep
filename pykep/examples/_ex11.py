def run_example11(n_seg=30):
    """
    This example demonstrates the use of the class lt_margo developed for the internal ESA CDF study on the 
    interplanetary mission named MARGO. The class was used to produce the preliminary traget selection
    for the mission resulting in 88 selected possible targets 
    (http://www.esa.int/spaceinimages/Images/2017/11/Deep-space_CubeSat)
    (http://www.esa.int/gsp/ACT/mad/projects/lt_to_NEA.html)
    """
    import pykep as pk
    import pygmo as pg
    import numpy as np
    from matplotlib import pyplot as plt
    from pykep.examples import add_gradient, algo_factory

    # 0 - Target asteroid from MPCORB
    mpcorbline = "K14Y00D 24.3   0.15 K1794 105.44160   34.12337  117.64264    1.73560  0.0865962  0.88781021   1.0721510  2 MPO369254   104   1  194 days 0.24 M-v 3Eh MPCALB     2803          2014 YD            20150618"

    # 1 - Algorithm
    algo = algo_factory("snopt7")

    # 2 - Problem
    udp = add_gradient(pk.trajopt.lt_margo(
        target=pk.planet.mpcorb(mpcorbline),
        n_seg=n_seg,
        grid_type="uniform",
        t0=[pk.epoch(8000), pk.epoch(9000)],
        tof=[200, 365.25 * 3],
        m0=20.0,
        Tmax=0.0017,
        Isp=3000.0,
        earth_gravity=False,
        sep=True,
        start="earth"),
        with_grad=False
    )
    prob = pg.problem(udp)
    prob.c_tol = [1e-5] * prob.get_nc()

    # 3 - Population
    pop = pg.population(prob, 1)

    # 4 - Solve the problem (evolve)
    pop = algo.evolve(pop)

    # 5 - Inspect the solution
    print("Feasible?:", prob.feasibility_x(pop.champion_x))

    # 6 - plot trajectory
    axis = udp.udp_inner.plot_traj(pop.champion_x, plot_segments=True)
    plt.title("The trajectory in the heliocentric frame")

    # 7 - plot control
    udp.udp_inner.plot_dists_thrust(pop.champion_x)

    # 8 - pretty
    udp.udp_inner.pretty(pop.champion_x)

    plt.ion()
    plt.show()


if __name__ == "__main__":
    run_example11()
