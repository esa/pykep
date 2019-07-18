from pykep.trajopt import mga_1dsm, launchers
from pykep.planet import jpl_lp

import numpy as np
from numpy.linalg import norm
from math import log, acos, cos, sin, asin, exp
from copy import deepcopy

class _juice_udp(mga_1dsm):
    def __init__(self):
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
            tof=[[100, 400], [100, 400], [100, 400], [100, 400], [500, 800], [500, 1500]],
            vinf=[2.5, 3.475],
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
