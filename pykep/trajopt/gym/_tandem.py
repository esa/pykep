from pykep.trajopt import mga_1dsm, launchers
from pykep.planet import jpl_lp
from pykep import epoch_from_string

import numpy as np
from numpy.linalg import norm
from math import log, acos, cos, sin, asin, exp
from copy import deepcopy

class _tandem_udp(mga_1dsm):
    """ 
    This class represents a rendezvous mission to Saturn modelled as an MGA-1DSM transfer. Mission parameters are
    inspired to the TandEM mission. A launcher model (i.e. Atlas 501) is also used, so that the final mass delivered
    at Saturn is the main objective of this optimization problem.

    The problem draws inspiration from the work performed in April 2008 by the
    European Space Agency working group on mission analysis on the mission named TandEM. TandEM is an interplanetary
    mission aimed at reaching Titan and Enceladus (two moons of Saturn). 

    .. note::

       The Titan and Enceladus Mission (TandEM), an ambitious scientific mission to study the Saturnian system
       with particular emphasis on the moons Titan and Enceladus, was selected in October 2007 as a candidate mission
       within the ESA Cosmic Vision plan. In February 2009, TandEM exited the Cosmic Vision programme when ESA
       and NASA chose EJSM-Laplace as the L-class outer Solar System mission candidate.

    .. note::

       A significantly similar version of this problem was part of the no longer maintained GTOP database, 
       https://www.esa.int/gsp/ACT/projects/gtop/gtop.html. The exact definition is, though, different and results
       cannot thus not be compared to those posted in GTOP.
    """
    def __init__(self, prob_id = 1, constrained = True):
        """ 
        The TandEM problem of the trajectory gym consists in 48 different instances varying in fly-by sequence and
        the presence of a time constraint.

        Args:
            - prob_id (int): The problem id defines the fly-by sequence.
            - constrained (bool): Activates the constraint on the time of flight 
              (fitness will thus return two numbers, the objectove function and the inequality constraint violation)
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

        super().__init__(
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
        DV = super().fitness(x)[0]
        DV = DV + 165.  # losses for 3 swgbys + insertion
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
        super().pretty(x)
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
        DV = super().fitness(x)[0]
        DV = DV + 165.  # losses for 3 swgbys + insertion
        m_final = m_initial * exp(-DV / Isp / g0)
        print("\nInitial mass:", m_initial)
        print("Final mass:", m_final)
        print("Declination:", declination)

    def __repr__(self):
        return "TandEM (Trajectory Optimisation Gym P12, multiple instances)"

# Problem P12: TandEM mission MGA1DSM, single objective, direct encoding, time constrained
tandem = _tandem_udp
