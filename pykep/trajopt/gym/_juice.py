from pykep.trajopt import mga_1dsm, launchers
from pykep.planet import jpl_lp
from pykep import epoch_from_string

import numpy as np
from numpy.linalg import norm
from math import log, acos, cos, sin, asin, exp
from copy import deepcopy

class _juice_udp(mga_1dsm):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        # Redefining the planets as to change their safe radius
        earth = jpl_lp('earth')
        earth.safe_radius = 1.05
        # We need the Earth eph in the fitnes
        venus = jpl_lp('venus')
        venus.safe_radius = 1.05
        mars = jpl_lp('mars')
        mars.safe_radius = 1.05
        jupiter = jpl_lp('jupiter')

        super(_juice_udp, self).__init__(
            seq=[earth, earth, venus, earth, mars, earth, jupiter],
            t0=[8000, 8400],
            tof=[[300, 400], [100, 200], [300, 400], [100, 200], [600, 700], [1000, 1100]],
            vinf=[1., 4.],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding='direct',
            multi_objective=False,
            orbit_insertion=True,
            e_target=0.98531407996358,
            rp_target=1070400000,
            eta_lb=0.01,
            eta_ub=0.99,
            rp_ub=10
        )

    def fitness(self, x):
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x 
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = - Vinfy * sin(earth_axis_inclination) + Vinfz * cos(earth_axis_inclination)
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.
        # We now have the initial mass of the spacecraft
        m_initial = launchers.ariane5(x[3] / 1000., declination)
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.
        g0 = 9.80665
        DV = super(_juice_udp, self).fitness(x)[0]
        DV = DV + 275.  # losses for 5 swingbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        # Numerical guard for the exponential
        if m_final == 0:
            m_final = 1e-320
        return [-log(m_final), ]

    def get_name(self):
        return "JUICE"

    def get_extra_info(self):
        retval = "\t Sequence: " + \
            [pl.name for pl in self._seq].__repr__()
        return retval

    def pretty(self, x):
        """
        prob.plot(x)

        - x: encoded trajectory

        Prints human readable information on the trajectory represented by the decision vector x

        Example::

          print(prob.pretty(x))
        """
        super(_juice_udp, self).pretty(x)
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x 
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = - Vinfy * sin(earth_axis_inclination) + Vinfz * cos(earth_axis_inclination)
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.
        m_initial = launchers.ariane5(x[3] / 1000., declination)
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.
        g0 = 9.80665
        DV = super(_juice_udp, self).fitness(x)[0]
        DV = DV + 275.  # losses for 5 swgbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        print("\nInitial mass:", m_initial)
        print("Final mass:", m_final)
        print("Declination:", declination)


# Problem P13: JUICE mission MGA1DSM, single objective, direct encoding
juice = _juice_udp()