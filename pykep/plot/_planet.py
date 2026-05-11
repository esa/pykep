import pykep as _pk
import numpy as _np


def add_planet(ax, pla: _pk.planet, when, units=_pk.AU, **kwargs):
    """Adds a planet to *ax*.  All kwargs are forwarded to the scatter method of *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *pla* (:class:`~pykep.planet`): the planet.

        *when* (:class:`~pykep.epoch` or :class:`float`): the epoch (in mjd2000 if float).

        *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.
        
        *\\*\\*kwargs*: Additional keyword arguments to pass to the Axes3D.plot function.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.

    Examples:
        >>> import pykep as pk
        >>> ax = pk.plot.make_3Daxis()
        >>> pk.plot.add_sun(ax, label="Sun")
        >>> udpla = pk.udpla.jpl_lp(body="EARTH")
        >>> earth = pk.planet(udpla)
        >>> pk.plot.add_planet_orbit(ax, earth, plot_range = [0, 365.25], c = "royalblue", label = "")
        >>> pk.plot.add_planet(ax, earth, when = pk.epoch(0), c = "royalblue")
    """
    # We compute the ephemerides
    r, _ = pla.eph(when)

    # And plot them as a scatter3D
    ax.scatter(r[0] / units, r[1] / units, r[2] / units, **kwargs)

    # Returning the axes.
    return ax

def add_planets(ax, plas: _pk.planet, when, units=_pk.AU, **kwargs):
    """Adds multiple planets to *ax*.  All kwargs are forwarded to the scatter method of *ax*.
    This method is significantly faster than using :class:`~pykep.planet.add_planet` in a loop.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *plas* (:class:list:[:class:`~pykep.planet`]): a list of planets.

        *when* (:class:`~pykep.epoch` or :class:`float`): the epoch (in mjd2000 if float).

        *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.
        
        *\\*\\*kwargs*: Additional keyword arguments to pass to the Axes3D.plot function.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.

    Examples:
        >>> import pykep as pk
        >>> ax = pk.plot.make_3Daxis()
        >>> pk.plot.add_sun(ax, label="Sun")
        >>> udpla = pk.udpla.jpl_lp(body="EARTH")
        >>> earth = pk.planet(udpla)
        >>> planets = [earth, ....]
        >>> pk.plot.add_planets(ax, planets, when = pk.epoch(0), c = "royalblue")
    """
    # We compute the ephemerides 
    pos = []
    for pla in plas:
        r, _ = pla.eph(when)
        pos.append(r)
    pos = _np.array(pos)

    # And plot them as a scatter3D
    ax.scatter(pos[:, 0] / units, pos[:, 1] / units, pos[:, 2] / units, **kwargs)

    # Returning the axes.
    return ax


def add_planet_orbit(
    ax,
    pla: _pk.planet,
    plot_range=None,
    units=_pk.AU,
    N=60,
    **kwargs
):
    """Adds a planet orbit to *ax*. All kwargs are forwarded to the plot method of *ax*.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *pla* (:class:`~pykep.planet`): the planet.

        *plot_range* (:class:`list`, optional): the starting and end mjd2000 to be plotted. Defaults to [0., one_orbital_period] if a period can be computed.

        *units* (:class:`float`, optional): length units to be used. Defaults to pk.AU.

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.

    Examples:
        >>> import pykep as pk
        >>> ax = pk.plot.make_3Daxis()
        >>> pk.plot.add_sun(ax, label="Sun")
        >>> udpla = pk.udpla.jpl_lp(body="EARTH")
        >>> earth = pk.planet(udpla)
        >>> pk.plot.add_planet_orbit(ax, earth, plot_range = [0, 365.25], c = "royalblue", label = "")
        >>> pk.plot.add_planet(ax, earth, when = pk.epoch(0), c = "royalblue")
    """

    # If the plot_range is not defined by the user, then a defult is attempted [0,T]
    if plot_range is None:
        try:
            T = pla.period() * _pk.SEC2DAY
        except NotImplementedError as e:
            print(
                f"PyKEP ERROR: Cannot compute the orbital period when plotting the planet {pla.get_name()}, please define an explicit plot_range."
            )
            raise e
        epochs = _np.linspace(0.0, T, N)
    else:
        epochs = _np.linspace(plot_range[0], plot_range[1], N)

    # Compute the ephemerides
    rvs = pla.eph_v(epochs)[:, :3] / units

    # Plots the planet in the range.
    ax.plot(rvs[:, 0], rvs[:, 1], rvs[:, 2], **kwargs)

    # Returning the axes.
    return ax


def add_solar_system(
    ax, bodies=[1, 2, 3, 4, 5, 6], when=_pk.epoch(0), s=[15, 2, 3, 3, 2, 8, 8, 5, 5]
):
    """Adds to *ax* the selected solar system planets.

    Args:
        *ax* (:class:`mpl_toolkits.mplot3d.axes3d.Axes3D`): the 3D axis.

        *bodies* (:class:`list`, optional): the solar system planet id (i.e. 3 for Earth). Defaults to [1, 2, 3, 4, 5, 6].

        *when* (:class:`~pykep.epoch`): the epoch.
        
        *s* (:class:`list`, optional): the size of the Sun and all the 8 solar system planets. Defaults to [15, 2, 3, 3, 2, 8, 8, 5, 5].

    Returns:
        :class:`mpl_toolkits.mplot3d.axes3d.Axes3D`: the 3D axis.
    """
    # Sun
    _pk.plot.add_sun(ax, s=s[0])

    # Mercury
    if 1 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Mercury")
        mercury = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, mercury, c="tomato", label="")
        _pk.plot.add_planet(ax, mercury, when=when, c="tomato", s=s[1], label = "Mercury")

    # Venus
    if 2 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Venus")
        venus = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, venus, c="forestgreen", label="")
        _pk.plot.add_planet(ax, venus, when=when, c="forestgreen", s=s[2], label = "Venus")

    # Earth
    if 3 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Earth")
        earth = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, earth, c="royalblue", label="")
        _pk.plot.add_planet(ax, earth, when=when, c="royalblue", s=s[3], label = "Earth")

    # Mars
    if 4 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Mars")
        mars = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, mars, c="indianred", label="")
        _pk.plot.add_planet(ax, mars, when=when, c="indianred", s=s[4], label = "Mars")

    # Jupiter
    if 5 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Jupiter")
        jupiter = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, jupiter, c="tan", label="")
        _pk.plot.add_planet(ax, jupiter, when=when, c="tan", s=s[5], label = "Jupiter")

    # Saturn
    if 6 in bodies:
        udpla = _pk.udpla.jpl_lp(body="Saturn")
        saturn = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(ax, saturn, c="darkcyan", label="")
        _pk.plot.add_planet(ax, saturn, when=when, c="darkcyan", s=s[6], label = "Saturn")

    # For these outer we use the de440s ephs, else no data would be available.
    # Keplerian planets would likely make more sense everywere in this function
    # Uranus
    if 7 in bodies:
        udpla = _pk.udpla.jpl_lp(body="uranus")
        uranus = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(
            ax,
            uranus,
            c="slateblue",
            label="",
            plot_range=[-40000, -40000 + 84.8 * 365.25],
        )
        _pk.plot.add_planet(ax, uranus, when=when, c="slateblue", s=s[7], label = "Uranus")

    # Neptune
    if 8 in bodies:
        udpla = _pk.udpla.jpl_lp(body="neptune")
        neptune = _pk.planet_to_keplerian(_pk.planet(udpla), when=when)
        _pk.plot.add_planet_orbit(
            ax,
            neptune,
            c="orchid",
            label="",
            plot_range=[-40000, -40000 + 164.8 * 365.25],
        )
        _pk.plot.add_planet(ax, neptune, when=ep, c="orchid", s=s[7], label = "Neptune")

    return ax
