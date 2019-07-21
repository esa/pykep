from pykep.trajopt import mga_1dsm
from pykep.planet import jpl_lp, keplerian
from pykep import AU, DEG2RAD, MU_SUN, epoch

# CASSINI2 (we need to modify the safe radius of the planets to match the wanted problem)
_earth_cassini2 = jpl_lp('earth')
_earth_cassini2.safe_radius = 1.15
_venus_cassini2 = jpl_lp('venus')
_venus_cassini2.safe_radius = 1.05
_jupiter_cassini2 = jpl_lp('jupiter')
_jupiter_cassini2.safe_radius = 1.7
_seq_cassini2 = [_earth_cassini2,
                _venus_cassini2,
                _venus_cassini2,
                _earth_cassini2,
                _jupiter_cassini2,
                jpl_lp('saturn')]


class _cassini2_udp(mga_1dsm):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
            seq =_seq_cassini2,
            t0 = [-1000, 0],
            tof = [[100, 400], [100, 500], [30, 300], [400, 1600], [800, 2200]],
            vinf = [3., 5.],
            add_vinf_dep = False,
            add_vinf_arr = True,
            tof_encoding = "direct",
            multi_objective = False,
            orbit_insertion=True,
            e_target=0.98,
            rp_target=108950000,
            eta_lb = 0.01,
            eta_ub = 0.9,
            rp_ub = 70.
        )

# Problem P11: Cassini mission MGA1DSM, single objective, direct encoding
cassini2 = _cassini2_udp()