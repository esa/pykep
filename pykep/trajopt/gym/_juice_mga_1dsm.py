import pykep as _pk
from pykep.trajopt import mga_1dsm as _mga_1dsm
from pykep.trajopt import launchers as _launchers


import numpy as np
from numpy.linalg import norm
from math import log, cos, sin, asin, exp


class _juice_udp(_mga_1dsm):
    """
    This class represents a rendezvous mission to Jupiter modelled as an MGA-1DSM transfer. The selected fly-by sequence,
    E-EVEME-J, and other parameters are inspired by the ESA Juice mission. A launcher model is included, namely an Ariane5
    launch from Kourou.
    JUICE - JUpiter ICy moons Explorer - is the first large-class mission in ESA's Cosmic Vision 2015-2025 programme.
    Launched on the 14th of April 2023, ESA’s Jupiter Icy Moons Explorer, Juice, will make detailed observations
    of the giant gas planet and its three large ocean-bearing moons – Ganymede, Callisto and Europa – with a suite of
    remote sensing, geophysical and in situ instruments. The mission will characterise these moons as both
    planetary objects and possible habitats, explore Jupiter’s complex environment in depth, and study the wider
    Jupiter system as an archetype for gas giants across the Universe.
    """

    def __init__(self, multi_objective, tof_encoding, tof):
        """
        Args:
            - multi_objective (``bool``): when True the problem fitness will return also the time of flight as an added objective
            - tof_encoding (``str``): one of 'direct', 'eta' or 'alpha'. Selects the encoding for the time of flights
            - tof (``list`` or ``list`` of ``list``): time of flight bounds. As documented in ``pykep.mga_1dsm``

        """
        # Redefining the planets as to change their safe radius
        udpla = _pk.udpla.jpl_lp("earth")
        udpla.safe_radius = 1.05 * udpla.radius
        earth = _pk.planet(udpla)
        # We need the Earth eph in the fitness
        udpla = _pk.udpla.jpl_lp("venus")
        udpla.safe_radius = 1.05 * udpla.radius
        venus = _pk.planet(udpla)
        udpla = _pk.udpla.jpl_lp("mars")
        udpla.safe_radius = 1.05 * udpla.radius
        mars = _pk.planet(udpla)
        udpla = _pk.udpla.jpl_lp("jupiter")
        jupiter = _pk.planet(udpla)

        super().__init__(
            seq=[earth, earth, venus, earth, mars, earth, jupiter],
            t0=[8000, 8400],
            tof=tof,
            vinf=[1.0, 4.0],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding=tof_encoding,
            multi_objective=multi_objective,
            orbit_insertion=True,
            e_target=0.98531407996358,
            rp_target=1070400000,
            eta_bounds=[0.01, 0.99],
            rp_ub=10,
        )

    def fitness(self, x):
        T, _, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = -Vinfy * sin(earth_axis_inclination) + Vinfz * cos(
            earth_axis_inclination
        )
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.0
        # We now have the initial mass of the spacecraft
        m_initial = _launchers.ariane5(x[3] / 1000.0, declination)[0, 0]
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.0
        g0 = 9.80665

        if self._multi_objective:
            DV, T = super().fitness(x)
        else:
            (DV,) = super().fitness(x)

        DV = DV + 275.0  # losses for 5 swingbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        # Numerical guard for the exponential
        if m_final == 0:
            m_final = 1e-320

        encoded_m_final = -log(m_final)

        if self._multi_objective:
            return (encoded_m_final, T)

        return (encoded_m_final,)

    def get_name(self):
        return "Juice (Trajectory Optimization Gym P13-14)"

    def __repr__(self):
        return self.get_name()

    def get_extra_info(self):
        retval = "\t Sequence: " + [pl.name for pl in self._seq].__repr__()
        return retval

    def pretty(self, x):
        """
        prob.plot(x)

        - x: encoded trajectory

        Prints human readable information on the trajectory represented by the decision vector x

        Example::

          print(prob.pretty(x))
        """
        super().pretty(x)
        T, _, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = -Vinfy * sin(earth_axis_inclination) + Vinfz * cos(
            earth_axis_inclination
        )
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.0
        m_initial = _launchers.ariane5(x[3] / 1000.0, declination)[0, 0]
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.0
        g0 = 9.80665
        DV = super().fitness(x)[0]
        DV = DV + 275.0  # losses for 5 swgbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        print("\nInitial mass:", m_initial)
        print("Final mass:", m_final)
        print("Declination:", declination)


# Problem P13: JUICE mission MGA1DSM, single objective, direct encoding
juice = _juice_udp(
    multi_objective=False,
    tof_encoding="direct",
    tof=[[200, 500], [30, 300], [200, 500], [30, 300], [500, 800], [900, 1200]],
)

# Problem P14: JUICE mission MGA1DSM, multiple objective, alpha encoding
juice_mo = _juice_udp(multi_objective=True, tof_encoding="alpha", tof=[2000, 3000])
