from pykep.trajopt import mga_1dsm, launchers
from pykep.planet import jpl_lp

import numpy as np
from numpy.linalg import norm
from math import log, acos, cos, sin, asin, exp
from copy import deepcopy

class _tandem_udp(mga_1dsm):
    def __init__(self, prob_id=1, constrained=True):
        # Redefining the planets as to change their safe radius
        earth = jpl_lp('earth')
        earth.safe_radius = 1.05
        # We need the Earth eph in the fitnes
        venus = jpl_lp('venus')
        venus.safe_radius = 1.05
        mars = jpl_lp('mars')
        mars.safe_radius = 1.05
        jupiter = jpl_lp('jupiter')
        # This safe radius is clearly too low, but
        # its set like this for consistency to the old GTOP database
        jupiter.safe_radius = 1.05
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

        self.prob_id = prob_id
        self.constrained = constrained
        super(_tandem_udp, self).__init__(
            seq=seq_tandem[prob_id - 1],
            t0=[5475, 9132],
            tof=[[20, 2500], [20, 2500], [20, 2500], [20, 2500]],
            vinf=[2.5, 4.9],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding='direct',
            multi_objective=False,
            eta_lb=0.01,
            eta_ub=0.99,
            rp_ub=10
        )

    def fitness(self, x):
        # The Earth position and velocity at x[0]
        rE, vE = self._seq[0].eph(x[0])
        # We define a reference system on the ecliptic plane
        vtemp = np.cross(rE, vE)
        iP1 = np.array(vE) / norm(vE)
        zP1 = vtemp / norm(vtemp)
        jP1 = np.cross(zP1, iP1)
        # And compute the outgoing hyperbolic velocity in it
        theta = 2. * np.pi * x[1]
        phi = acos(2. * x[2] - 1.) - np.pi / 2.
        vinf = x[3] * (cos(theta) * cos(phi) * iP1 + sin(theta)
                       * cos(phi) * jP1 + sin(phi) * zP1)
        # We transform it (only the needed component) to an equatorial system
        incl = 0.409072975
        vinf[2] = vinf[1] * sin(incl) + vinf[2] * cos(incl)
        # And we find the vinf declination (in degrees)
        sindelta = vinf[2] / x[3]
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
            T = self._decode_times_and_vinf(x)[0]
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
