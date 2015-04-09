def run_example5():
    from PyGMO import archipelago, problem
    from PyGMO.algorithm import jde
    from PyGMO.topology import ring
    from PyKEP import epoch
    from PyKEP.planet import jpl_lp
    from PyKEP.trajopt import mga_1dsm

    # We define an Earth-Venus-Earth problem (single-objective)
    seq = [jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')]
    prob = mga_1dsm(seq=seq)

    prob.set_tof(0.7, 3)
    prob.set_vinf(2.5)
    prob.set_launch_window(epoch(5844), epoch(6209))
    prob.set_tof(0.7, 3)

    # We solve it!!
    algo = jde(100)
    topo = ring()
    archi = archipelago(algo, prob, 8, 20, topology=topo)
    print(
        "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands")
    archi.evolve(10)
    archi.join()
    isl = min(archi, key=lambda x: x.population.champion.f[0])
    print("Done!! Best solution found is: " +
          str(isl.population.champion.f[0] / 1000) + " km / sec")
    prob.pretty(isl.population.champion.x)
    prob.plot(isl.population.champion.x)
