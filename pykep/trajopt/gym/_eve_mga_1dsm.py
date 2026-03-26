import pykep as pk
from pykep.trajopt import mga_1dsm as _mga_1dsm

# Earth-Venus_Earth
_seq_eve = [
    pk.planet(pk.udpla.jpl_lp("earth")),
    pk.planet(pk.udpla.jpl_lp("venus")),
    pk.planet(pk.udpla.jpl_lp("earth")),
]


class _eve_mga_1dsm_udp(_mga_1dsm):
    def __init__(
        self,
        tof_encoding="direct",
        t0=[pk.epoch(0), pk.epoch(3000)],
        tof=[[10, 500], [10, 500]],
    ):
        super().__init__(
            seq=_seq_eve,
            t0=t0,
            tof=tof,
            vinf=[0.5, 2.5],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding=tof_encoding,
            multi_objective=False,
        )

    def get_name(self):
        return "Earth-Venus-Earth mga-1dsm (Trajectory Optimisation Gym P7-9)"

    def get_extra_info(self):
        retval = "\ttof_encoding: " + self._tof_encoding
        return retval

    def __repr__(self):
        return self.get_name()


# Problem P7: Earth-Venus-Earth MGA1DSM, single objective , direct encoding
eve_mga1dsm = _eve_mga_1dsm_udp(
    tof_encoding="direct",
    t0=[pk.epoch(0), pk.epoch(3000)],
    tof=[[10, 500], [10, 500]],
)
# Problem P8: Earth-Venus-Earth MGA1DSM, single objective , alpha encoding
eve_mga1dsm_a = _eve_mga_1dsm_udp(
    tof_encoding="alpha", t0=[pk.epoch(0), pk.epoch(3000)], tof=[300, 700]
)
# Problem P9: Earth-Venus-Earth MGA1DSM, single objective , eta encoding
eve_mga1dsm_n = _eve_mga_1dsm_udp(
    tof_encoding="eta", t0=[pk.epoch(0), pk.epoch(3000)], tof=700
)
