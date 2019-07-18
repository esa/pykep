from pykep.trajopt import mga_1dsm, launchers
from pykep.planet import jpl_lp
from pykep import epoch_from_string

import numpy as np
from numpy.linalg import norm
from math import log, acos, cos, sin, asin, exp
from copy import deepcopy

"""
IMPORTANT: This is not the same as the old TandEM (GTOP database) problem and thus values should not be compared
"""
class _tandem_udp(mga_1dsm):
    def __init__(self, prob_id = 1, constrained = True):
        # Redefining the planets as to change their safe radius
        earth = jpl_lp('earth')
        earth.safe_radius = 1.05
        # We need the Earth eph in the fitnes
        venus = jpl_lp('venus')
        venus.safe_radius = 1.05
        mars = jpl_lp('mars')
        mars.safe_radius = 1.05
        jupiter = jpl_lp('jupiter')
        jupiter.safe_radius = 1.7
        saturn = jpl_lp('saturn')

        # Defining the different sequences
        seq_tandem = []
        seq_tandem.append([earth, venus, venus, venus, saturn])
        seq_tandem.append([earth, venus, venus, earth, saturn])
        seq_tandem.append([earth, venus, venus, mars, saturn])
        seq_tandem.append([earth, venus, venus, jupiter, saturn])

        seq_tandem.append([earth, venus, earth, venus, saturn])
        seq_tandem.append([earth, venus, earth, earth, saturn])
        seq_tandem.append([earth, venus, earth, mars, saturn])
        seq_tandem.append([earth, venus, earth, jupiter, saturn])

        seq_tandem.append([earth, venus, mars, venus, saturn])
        seq_tandem.append([earth, venus, mars, earth, saturn])
        seq_tandem.append([earth, venus, mars, mars, saturn])
        seq_tandem.append([earth, venus, mars, jupiter, saturn])

        seq_tandem.append([earth, earth, venus, venus, saturn])
        seq_tandem.append([earth, earth, venus, earth, saturn])
        seq_tandem.append([earth, earth, venus, mars, saturn])
        seq_tandem.append([earth, earth, venus, jupiter, saturn])

        seq_tandem.append([earth, earth, earth, venus, saturn])
        seq_tandem.append([earth, earth, earth, earth, saturn])
        seq_tandem.append([earth, earth, earth, mars, saturn])
        seq_tandem.append([earth, earth, earth, jupiter, saturn])

        seq_tandem.append([earth, earth, mars, venus, saturn])
        seq_tandem.append([earth, earth, mars, earth, saturn])
        seq_tandem.append([earth, earth, mars, mars, saturn])
        seq_tandem.append([earth, earth, mars, jupiter, saturn])

        if prob_id > 24 or type(prob_id) != int or prob_id < 1:
            raise ValueError("TandEM problem id must be an integer in [1, 24]")

        super(_tandem_udp, self).__init__(
            seq = seq_tandem[prob_id - 1],
            t0 = [5475, 9132],
            tof = [[20, 2500], [20, 2500], [20, 2500], [20, 2500]],
            vinf = [2.5, 4.9],
            add_vinf_dep = False,
            add_vinf_arr = True,
            tof_encoding = 'direct',
            multi_objective = False,
            orbit_insertion = True,
            e_target = 0.98531407996358,
            rp_target = 80330000,
            eta_lb = 0.01,
            eta_ub = 0.99,
            rp_ub = 10
        )

        self.prob_id = prob_id
        self.constrained = constrained

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
        m_initial = launchers.atlas501(x[3] / 1000., declination)
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.
        g0 = 9.80665
        DV = super(_tandem_udp, self).fitness(x)[0]
        DV = DV + 0.165  # losses for 3 swgbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        # Numerical guard for the exponential
        if m_final == 0:
            m_final = 1e-320
        if self.constrained:
            retval = [-log(m_final), sum(T) - 3652.5]
        else:
            retval = [-log(m_final), ]
        return retval

    def get_nic(self):
        return int(self.constrained)

    def get_name(self):
        return "TandEM problem, id: " + str(self.prob_id)

    def get_extra_info(self):
        retval = "\t Sequence: " + \
            [pl.name for pl in self._seq].__repr__() + "\n\t Constrained: " + \
            str(self.constrained)
        return retval
    
    def pretty(self, x):
        """
        prob.plot(x)

        - x: encoded trajectory

        Prints human readable information on the trajectory represented by the decision vector x

        Example::

          print(prob.pretty(x))
        """
        super(_tandem_udp, self).pretty(x)
        T, Vinfx, Vinfy, Vinfz = self._decode_times_and_vinf(x)
        # We transform it (only the needed component) to an equatorial system rotating along x 
        # (this is an approximation, assuming vernal equinox is roughly x and the ecliptic plane is roughly xy)
        earth_axis_inclination = 0.409072975
        # This is different from the GTOP tanmEM problem, I think it was bugged there as the rotation was in the wrong direction.
        Vinfz = - Vinfy * sin(earth_axis_inclination) + Vinfz * cos(earth_axis_inclination)
        # And we find the vinf declination (in degrees)
        sindelta = Vinfz / x[3]
        declination = asin(sindelta) / np.pi * 180.
        m_initial = launchers.atlas501(x[3] / 1000., declination)
        # And we can evaluate the final mass via Tsiolkowsky
        Isp = 312.
        g0 = 9.80665
        DV = super(_tandem_udp, self).fitness(x)[0]
        DV = DV + 0.165  # losses for 3 swgbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        print("\nInitial mass:", m_initial)
        print("Final mass:", m_final)
        print("Declination:", declination)


