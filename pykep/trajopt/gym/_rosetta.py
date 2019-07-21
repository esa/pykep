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
    This class represents a rendezvous mission to the comet 67P/Churyumov-Gerasimenko modelled as an MGA-1DSM transfer.
    The fly-by sequence selected (i.e. E-EMEE-C) is similar to the one planned for the spacecraft Rosetta. 
    The objective function considered is the total mission delta V. No launcher model is employed and a final randezvous
    with the comet is included in the delta V computations.

    .. note::

       A significantly similar version of this problem was part of the no longer maintained GTOP database, 
       https://www.esa.int/gsp/ACT/projects/gtop/gtop.html. The exact definition is, though, different and results
       should thus not be compared to those posted in GTOP.
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
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