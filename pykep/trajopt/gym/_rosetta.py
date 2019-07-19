from pykep.trajopt import mga_1dsm
from pykep.planet import jpl_lp, keplerian
from pykep import AU, DEG2RAD, MU_SUN, epoch

# ROSETTA (we need to modify the safe radius of the planets to match the wanted problem)
_churyumov = keplerian(epoch(52504.23754000012, "mjd"),
                                 [3.50294972836275 * AU,
                                  0.6319356,
                                  7.12723 * DEG2RAD,
                                  50.92302 * DEG2RAD,
                                  11.36788 * DEG2RAD,
                                  0. * DEG2RAD],
                                 MU_SUN, 0., 0., 0., "Churyumov-Gerasimenko")

_mars_rosetta = jpl_lp('mars')
_mars_rosetta.safe_radius = 1.05
_seq_rosetta = [jpl_lp('earth'), 
                jpl_lp('earth'), 
                _mars_rosetta, 
                jpl_lp('earth'), 
                jpl_lp('earth'), 
                _churyumov]

class _rosetta_udp(mga_1dsm):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        super(_rosetta_udp, self).__init__(
            seq =_seq_rosetta,
            t0 = [1460, 1825],
            tof = [[300, 500], [150, 800], [150, 800], [300, 800], [700, 1850]],
            vinf = [3., 5.],
            add_vinf_dep = False,
            add_vinf_arr = True,
            tof_encoding = "direct",
            multi_objective = False,
            eta_lb = 0.01,
            eta_ub = 0.9,
            rp_ub = 9.
        )

# Problem P10: Rosetta mission MGA1DSM, single objective, direct encoding
rosetta = _rosetta_udp()