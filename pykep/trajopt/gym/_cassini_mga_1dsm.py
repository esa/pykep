import pykep as pk
from pykep.trajopt import mga_1dsm as _mga_1dsm

# CASSINI2 (we need to modify the safe radius of the planets to match the wanted problem)
_earth_cassini2 = pk.udpla.jpl_lp("earth")
_earth_cassini2.safe_radius = 1.15 * _earth_cassini2.radius
_venus_cassini2 = pk.udpla.jpl_lp("venus")
_venus_cassini2.safe_radius = 1.05 * _venus_cassini2.radius
_jupiter_cassini2 = pk.udpla.jpl_lp("jupiter")
_jupiter_cassini2.safe_radius = 1.7 * _jupiter_cassini2.radius

_seq_cassini2 = [
    pk.planet(_earth_cassini2),
    pk.planet(_venus_cassini2),
    pk.planet(_venus_cassini2),
    pk.planet(_earth_cassini2),
    pk.planet(_jupiter_cassini2),
    pk.planet(pk.udpla.jpl_lp("saturn")),
]


class _cassini2_udp(_mga_1dsm):
    def __init__(self):
        super().__init__(
            seq=_seq_cassini2,
            t0=[-1000, 0],
            tof=[[100, 400], [100, 500], [30, 300], [400, 1600], [800, 2200]],
            vinf=[3.0, 5.0],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding="direct",
            multi_objective=False,
            orbit_insertion=True,
            e_target=0.98,
            rp_target=108950000,
            eta_bounds=[0.01, 0.9],
            rp_ub=70.0,
        )

    def get_name(self):
        return "Cassini MGA-1DSM direct tof encoding (Trajectory Optimisation Gym P11)"

    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P11): Cassini MGA-1DSM, single objective, direct encoding\n"
        retval += "\tPlanetary sequence" + str([pl.get_name() for pl in _seq_cassini2])
        return retval

    def __repr__(self):
        return self.get_name()


# Problem P11: Cassini mission MGA1DSM, single objective, direct encoding
cassini2 = _cassini2_udp()
