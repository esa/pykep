def plot_planet(plnt, t0='pykep.epoch(0)', N=60, units=1.0, color='k', alpha=1.0, s=40, legend=False, axes=None):
    """
    ax = plot_planet(plnt, t0='pykep.epoch(0)', N=60, units=1.0, color='k', s=40, legend=False, axes=None)

    - axes:      3D axis object created using fig.gca(projection='3d')
    - plnt:      pykep.planet object we want to plot
    - t0:        pykep.epoch object indicating when we want to plot the planet position
    - units:     the length unit to be used in the plot
    - color:     matplotlib color to use to plot the line
    - s:         planet size (pixel^2)
    - legend     when True plots the legend with the planet name and the epoch

    Plots the planet position and its orbit

    Example::

	import pykep as pk
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.gca(projection = '3d')
	pl = pk.planet.jpl_lp('earth')
	t_plot = pk.epoch(219)
	ax = pk.orbit_plots.plot_planet(pl, ax = ax, color='b')
    """
    from pykep import MU_SUN, SEC2DAY, epoch, AU, RAD2DEG
    from pykep.planet import keplerian
    from math import pi, sqrt
    import numpy as np
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    if t0 == 'pykep.epoch(0)':
        t0 = epoch(0)
    if type(t0) is not epoch:
        t0 = epoch(t0)

    # orbit period at epoch
    T = plnt.compute_period(t0) * SEC2DAY

    # make an osculating copy of the planet
    el = list(plnt.osculating_elements(t0))
    plnt_k = keplerian(t0, el, plnt.mu_central_body, plnt.mu_self, plnt.radius, plnt.safe_radius, plnt.name)

    # points where the orbit will be plotted
    when = np.linspace(0, T, N)

    # Ephemerides Calculation for the given planet
    x = np.array([0.0] * N)
    y = np.array([0.0] * N)
    z = np.array([0.0] * N)

    for i, day in enumerate(when):
        r, v = plnt_k.eph(epoch(t0.mjd2000 + day))
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units

    # Actual plot commands
    if legend:
        label = plnt.name + " " + t0.__repr__()[0:11]
    else:
        label = None
    ax.plot(x, y, z, label=label, c=color, alpha=alpha)
    ax.scatter([x[0]], [y[0]], [z[0]], s=s, marker='o', alpha=0.8, c=[color])

    if legend:
        ax.legend()
    return ax


def plot_lambert(l, N=60, sol=0, units=1.0, color='b', legend=False, axes=None, alpha=1.):
    """
    ax = plot_lambert(l, N=60, sol=0, units='pykep.AU', legend='False', axes=None, alpha=1.)

    - axes:       3D axis object created using fig.gca(projection='3d')
    - l:        pykep.lambert_problem object
    - N:		number of points to be plotted along one arc
    - sol:      solution to the Lambert's problem we want to plot (must be in 0..Nmax*2)
                where Nmax is the maximum number of revolutions for which there exist a solution.
    - units:    the length unit to be used in the plot
    - color:    matplotlib color to use to plot the line
    - legend:   when True it plots also the legend with info on the Lambert's solution chosen

    Plots a particular solution to a Lambert's problem

    Example::

      import pykep as pk
      import matplotlib.pyplot as plt

      fig = plt.figure()
      ax = fig.gca(projection='3d')

      t1 = pk.epoch(0)
      t2 = pk.epoch(640)
      dt = (t2.mjd2000 - t1.mjd2000) * pk.DAY2SEC

      pl = pk.planet.jpl_lp('earth')
      pk.orbit_plots.plot_planet(pl, t0=t1, axes=ax, color='k')
      rE,vE = pl.eph(t1)

      pl = pk.planet.jpl_lp('mars')
      pk.orbit_plots.plot_planet(pl, t0=t2, axes=ax, color='r')
      rM, vM = pl.eph(t2)

      l = lambert_problem(rE,rM,dt,pk.MU_SUN)
      pk.orbit_plots.plot_lambert(l, ax=ax, color='b')
      pk.orbit_plots.plot_lambert(l, sol=1, axes=ax, color='g')
      pk.orbit_plots.plot_lambert(l, sol=2, axes=ax, color='g', legend = True)

      plt.show()
    """
    from pykep import propagate_lagrangian, AU
    import numpy as np
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    if sol > l.get_Nmax() * 2:
        raise ValueError("sol must be in 0 .. NMax*2 \n * Nmax is the maximum number of revolutions for which there exist a solution to the Lambert's problem \n * You can compute Nmax calling the get_Nmax() method of the lambert_problem object")

    # We extract the relevant information from the Lambert's problem
    r = l.get_r1()
    v = l.get_v1()[sol]
    T = l.get_tof()
    mu = l.get_mu()

    # We define the integration time ...
    dt = T / (N - 1)

    # ... and alocate the cartesian components for r
    x = np.array([0.0] * N)
    y = np.array([0.0] * N)
    z = np.array([0.0] * N)

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v = propagate_lagrangian(r, v, dt, mu)

    # And we plot
    if legend:
        label = 'Lambert solution (' + str((sol + 1) // 2) + ' revs.)'
    else:
        label = None
    ax.plot(x, y, z, c=color, label=label, alpha=alpha)

    if legend:
        ax.legend()

    if axes is None:  # show only if axis is not set
        plt.show()
    return ax


def plot_kepler(r, v, t, mu, N=60, units=1, color='b', legend=False, axes=None):
    """
    ax = plot_kepler(r, v, t, mu, N=60, units=1, color='b', legend=False, axes=None):

    - axes:     3D axis object created using fig.gca(projection='3d')
    - r:        initial position (cartesian coordinates)
    - v:		initial velocity (cartesian coordinates)
    - t:		propagation time
    - mu:		gravitational parameter
    - N:		number of points to be plotted along one arc
    - units:	the length unit to be used in the plot
    - color:	matplotlib color to use to plot the line
    - legend	when True it plots also the legend

    Plots the result of a keplerian propagation
    """

    from pykep import propagate_lagrangian
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    # We define the integration time ...
    dt = t / (N - 1)

    # ... and calculate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v = propagate_lagrangian(r, v, dt, mu)

    # And we plot
    if legend:
        label = 'ballistic arc'
    else:
        label = None
    ax.plot(x, y, z, c=color, label=label)

    if legend:
        ax.legend()

    if axes is None:  # show only if axis is not set
        plt.show()
    return ax


def plot_taylor(r, v, m, u, t, mu, veff, N=60, units=1, color='b', legend=False, axes=None):
    """
    ax = plot_taylor(r, v, m, u, t, mu, veff, N=60, units=1, color='b', legend=False, axes=None):

    - axes:		3D axis object created using fig.gca(projection='3d')
    - r:		initial position (cartesian coordinates)
    - v:		initial velocity (cartesian coordinates)
    - m: 		initial mass
    - u:		cartesian components for the constant thrust
    - t:		propagation time
    - mu:		gravitational parameter
    - veff:	the product Isp * g0
    - N:		number of points to be plotted along one arc
    - units:	the length unit to be used in the plot
    - color:	matplotlib color to use to plot the line
    - legend:	when True it plots also the legend

    Plots the result of a taylor propagation of constant thrust
    """

    from pykep import propagate_taylor
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    # We define the integration time ...
    dt = t / (N - 1)

    # ... and calcuate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v, m = propagate_taylor(r, v, m, u, dt, mu, veff, -10, -10)

    # And we plot
    if legend:
        label = 'constant thrust arc'
    else:
        label = None
    ax.plot(x, y, z, c=color, label=label)

    if legend:
        ax.legend()

    if axes is None:  # show only if axis is not set
        plt.show()
    return ax


def plot_taylor_disturbance(r, v, m, thrust, disturbance, t, mu, veff, N=60, units=1, color='b', legend=False, axes=None):
    """
    ax = plot_taylor_disturbance(r, v, m, thrust, disturbance, t, mu, veff, N=60, units=1, color='b', legend=False, axes=None):

    - axes:			3D axis object created using fig.gca(projection='3d')
    - r:			initial position (cartesian coordinates)
    - v:			initial velocity (cartesian coordinates)
    - m: 			initial mass
    - thrust:		cartesian components for the constant thrust
    - disturbance:	cartesian components for the constant disturbance
    - t:			propagation time
    - mu:			gravitational parameter
    - veff:			the product Isp * g0
    - N:			number of points to be plotted along one arc
    - units:		the length unit to be used in the plot
    - color:		matplotlib color to use to plot the line
    - legend:		when True it plots also the legend

    Plots the result of a taylor propagation of constant thrust
    """

    from pykep import propagate_taylor_disturbance
    import matplotlib.pyplot as plt

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    # We define the integration time ...
    dt = t / (N - 1)

    # ... and calcuate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v, m = propagate_taylor_disturbance(
            r, v, m, thrust, disturbance, dt, mu, veff, -10, -10)

    # And we plot
    if legend:
        label = 'constant thrust arc'
    else:
        label = None
    ax.plot(x, y, z, c=color, label=label)

    if legend:
        ax.legend()

    if axes is None:  # show only if axis is not set
        plt.show()
    return ax


def plot_sf_leg(leg, N=5, units=1, color='b', legend=False, plot_line=True, plot_segments=False, axes=None):
    """
    ax = plot_sf_leg(leg, N=5, units=1, color='b', legend=False, no_trajectory=False, axes=None):

    - axes:		    3D axis object created using fig.gca(projection='3d')
    - leg:	        a pykep.sims_flanagan.leg
    - N:		    number of points to be plotted along one arc
    - units:	    the length unit to be used in the plot
    - color:	    matplotlib color to use to plot the trajectory and the grid points
    - legend	    when True it plots also the legend
    - plot_line: 	when True plots also the trajectory (between mid-points and grid points)

    Plots a Sims-Flanagan leg

    Example::

        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        t1 = epoch(0)
        pl = planet_ss('earth')
        rE,vE = pl.eph(t1)
        plot_planet(pl,t0=t1, units=AU, axes=ax)

        t2 = epoch(440)
        pl = planet_ss('mars')
        rM, vM = pl.eph(t2)
        plot_planet(pl,t0=t2, units=AU, axes=ax)

        sc = sims_flanagan.spacecraft(4500,0.5,2500)
        x0 = sims_flanagan.sc_state(rE,vE,sc.mass)
        xe = sims_flanagan.sc_state(rM,vM,sc.mass)
        l = sims_flanagan.leg(t1,x0,[1,0,0]*5,t2,xe,sc,MU_SUN)

        plot_sf_leg(l, units=AU, axes=ax)
    """
    from pykep import propagate_lagrangian, AU, DAY2SEC, G0, propagate_taylor
    import numpy as np
    from scipy.linalg import norm
    from math import exp
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D

    if axes is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        ax = axes

    # We compute the number of segments for forward and backward propagation
    n_seg = len(leg.get_throttles())
    fwd_seg = (n_seg + 1) // 2
    back_seg = n_seg // 2

    # We extract information on the spacecraft
    sc = leg.get_spacecraft()
    isp = sc.isp
    max_thrust = sc.thrust

    # And on the leg
    throttles = leg.get_throttles()
    mu = leg.get_mu()

    # Forward propagation

    # x,y,z contain the cartesian components of all points (grid+midpints)
    x = [0.0] * (fwd_seg * 2 + 1)
    y = [0.0] * (fwd_seg * 2 + 1)
    z = [0.0] * (fwd_seg * 2 + 1)

    state = leg.get_xi()

    # Initial conditions
    r = state.r
    v = state.v
    m = state.m
    x[0] = r[0] / units
    y[0] = r[1] / units
    z[0] = r[2] / units

    # We compute all points by propagation
    for i, t in enumerate(throttles[:fwd_seg]):
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        # Keplerian propagation and dV application
        if leg.high_fidelity is False:
            dV = [max_thrust / m * dt * dumb for dumb in t.value]
            if plot_line:
                plot_kepler(r, v, dt / 2, mu, N=N, units=units,
                            color=(alpha, 0, 1 - alpha), axes=ax)
            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            x[2 * i + 1] = r[0] / units
            y[2 * i + 1] = r[1] / units
            z[2 * i + 1] = r[2] / units
            # v= v+dV
            v = [a + b for a, b in zip(v, dV)]
            if plot_line:
                plot_kepler(r, v, dt / 2, mu, N=N, units=units,
                            color=(alpha, 0, 1 - alpha), axes=ax)
            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            x[2 * i + 2] = r[0] / units
            y[2 * i + 2] = r[1] / units
            z[2 * i + 2] = r[2] / units
            m *= exp(-norm(dV) / isp / G0)
        # Taylor propagation of constant thrust u
        else:
            u = [max_thrust * dumb for dumb in t.value]
            if plot_line:
                plot_taylor(r, v, m, u, dt / 2, mu, isp * G0, N=N,
                            units=units, color=(alpha, 0, 1 - alpha), axes=ax)
            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp * G0, -12, -12)
            x[2 * i + 1] = r[0] / units
            y[2 * i + 1] = r[1] / units
            z[2 * i + 1] = r[2] / units
            if plot_line:
                plot_taylor(r, v, m, u, dt / 2, mu, isp * G0, N=N,
                            units=units, color=(alpha, 0, 1 - alpha), axes=ax)
            r, v, m = propagate_taylor(
                r, v, m, u, dt / 2, mu, isp * G0, -12, -12)
            x[2 * i + 2] = r[0] / units
            y[2 * i + 2] = r[1] / units
            z[2 * i + 2] = r[2] / units

    x_grid = x[::2]
    y_grid = y[::2]
    z_grid = z[::2]
    x_midpoint = x[1::2]
    y_midpoint = y[1::2]
    z_midpoint = z[1::2]
    if plot_segments:
        ax.scatter(x_grid[:-1], y_grid[:-1], z_grid[:-1],
                     label='nodes', marker='o')
        ax.scatter(x_midpoint, y_midpoint, z_midpoint,
                     label='mid-points', marker='x')
        ax.scatter(x_grid[-1], y_grid[-1], z_grid[-1],
                     marker='^', c='y', label='mismatch point')

    # Backward propagation

    # x,y,z will contain the cartesian components of
    x = [0.0] * (back_seg * 2 + 1)
    y = [0.0] * (back_seg * 2 + 1)
    z = [0.0] * (back_seg * 2 + 1)

    state = leg.get_xf()

    # Final conditions
    r = state.r
    v = state.v
    m = state.m
    x[-1] = r[0] / units
    y[-1] = r[1] / units
    z[-1] = r[2] / units

    for i, t in enumerate(throttles[-1:-back_seg - 1:-1]):
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        if leg.high_fidelity is False:
            dV = [max_thrust / m * dt * dumb for dumb in t.value]
            if plot_line:
                plot_kepler(r, v, -dt / 2, mu, N=N, units=units,
                            color=(alpha, 0, 1 - alpha), axes=ax)
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            x[-2 * i - 2] = r[0] / units
            y[-2 * i - 2] = r[1] / units
            z[-2 * i - 2] = r[2] / units
            # v= v+dV
            v = [a - b for a, b in zip(v, dV)]
            if plot_line:
                plot_kepler(r, v, -dt / 2, mu, N=N, units=units,
                            color=(alpha, 0, 1 - alpha), axes=ax)
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            x[-2 * i - 3] = r[0] / units
            y[-2 * i - 3] = r[1] / units
            z[-2 * i - 3] = r[2] / units
            m *= exp(norm(dV) / isp / G0)
        else:
            u = [max_thrust * dumb for dumb in t.value]
            if plot_line:
                plot_taylor(r, v, m, u, -dt / 2, mu, isp * G0, N=N,
                            units=units, color=(alpha, 0, 1 - alpha), axes=ax)
            r, v, m = propagate_taylor(
                r, v, m, u, -dt / 2, mu, isp * G0, -12, -12)
            x[-2 * i - 2] = r[0] / units
            y[-2 * i - 2] = r[1] / units
            z[-2 * i - 2] = r[2] / units
            if plot_line:
                plot_taylor(r, v, m, u, -dt / 2, mu, isp * G0, N=N,
                            units=units, color=(alpha, 0, 1 - alpha), axes=ax)
            r, v, m = propagate_taylor(
                r, v, m, u, -dt / 2, mu, isp * G0, -12, -12)
            x[-2 * i - 3] = r[0] / units
            y[-2 * i - 3] = r[1] / units
            z[-2 * i - 3] = r[2] / units

    x_grid = x[::2]
    y_grid = y[::2]
    z_grid = z[::2]
    x_midpoint = x[1::2]
    y_midpoint = y[1::2]
    z_midpoint = z[1::2]

    if plot_segments:
        ax.scatter(x_grid[1:], y_grid[1:], z_grid[
                     1:], marker='o', label='nodes')
        ax.scatter(x_midpoint, y_midpoint, z_midpoint,
                     marker='x', label='mid-points')
        ax.scatter(x_grid[0], y_grid[0], z_grid[0],
                     marker='^', c='y', label='mismatch point')

    if legend:
        ax.legend()

    if axes is None:  # show only if axis is not set
        plt.show()
    return ax
