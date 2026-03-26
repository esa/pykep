import pykep as pk
from pykep.trajopt import mga_1dsm as _mga_1dsm

# ROSETTA (we need to modify the safe radius of the planets to match the wanted problem)
_churyumov = pk.udpla.keplerian(
        when=pk.epoch(52504.23754000012, pk.epoch.julian_type.MJD),
        elem=[
            3.50294972836275 * pk.AU,
            0.6319356,
            7.12723 * pk.DEG2RAD,
            50.92302 * pk.DEG2RAD,
            11.36788 * pk.DEG2RAD,
            0.0 * pk.DEG2RAD,
        ],
        mu_central_body=pk.MU_SUN,
        name="Churyumov-Gerasimenko",
    )

_mars_rosetta = pk.udpla.jpl_lp("mars")
_mars_rosetta.safe_radius = 1.05 * _mars_rosetta.radius

_seq_rosetta = [
    pk.planet(pk.udpla.jpl_lp("earth")),
    pk.planet(pk.udpla.jpl_lp("earth")),
    pk.planet(_mars_rosetta),
    pk.planet(pk.udpla.jpl_lp("earth")),
    pk.planet(pk.udpla.jpl_lp("earth")),
    pk.planet(_churyumov),
]


class _rosetta_udp(_mga_1dsm):
    """
    This class represents a rendezvous mission to the comet 67P/Churyumov-Gerasimenko modelled as an MGA-1DSM transfer.
    The fly-by sequence selected (i.e. E-EMEE-C) is similar to the one planned for the spacecraft Rosetta.
    The objective function considered is the total mission delta V. No launcher model is employed and a final rendezvous
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
            seq=_seq_rosetta,
            t0=[1460, 1825],
            tof=[[300, 500], [150, 800], [150, 800], [300, 800], [700, 1850]],
            vinf=[3.0, 5.0],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding="direct",
            multi_objective=False,
            eta_bounds=[0.01, 0.9],
            rp_ub=9.0,
        )

    def get_name(self):
        return "Rosetta (Trajectory Optimisation Gym P10)"

    def __repr__(self):
        return self.get_name()


# Problem P10: Rosetta mission MGA1DSM, single objective, direct encoding
rosetta = _rosetta_udp()
