def run_example2():
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D

    import matplotlib.pyplot as plt
    from PyKEP import epoch, DAY2SEC, AU, MU_SUN, lambert_problem
    from PyKEP.planet import jpl_lp
    from PyKEP.orbit_plots import plot_planet, plot_lambert

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    axis = fig.gca(projection='3d')

    t1 = epoch(0)
    t2 = epoch(640)
    dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

    axis.scatter([0], [0], [0], color='y')

    pl = jpl_lp('earth')
    plot_planet(
        pl, t0=t1, color=(0.8, 0.8, 1), legend=True, units = AU, ax = axis)
    rE, vE = pl.eph(t1)

    pl = jpl_lp('mars')
    plot_planet(
        pl, t0=t2, color=(0.8, 0.8, 1), legend=True, units = AU, ax = axis)
    rM, vM = pl.eph(t2)

    l = lambert_problem(rE, rM, dt, MU_SUN)
    plot_lambert(l, color='b', legend=True, units=AU, ax=axis)
    plot_lambert(l, sol=1, color='g', legend=True, units=AU, ax=axis)
    plot_lambert(l, sol=2, color='g', legend=True, units=AU, ax=axis)

    plt.show()
