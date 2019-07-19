from pykep.trajopt import mga_1dsm
from pykep import epoch
from pykep.planet import jpl_lp


class _eve_mga1dsm_udp(mga_1dsm):
    """
    Write Me
    """
    def __init__(self, tof_encoding='direct', t0=[epoch(0), epoch(3000)], tof=[[10, 500], [10, 500]]):
        """
        Write Me
        """
        super(_eve_mga1dsm_udp, self).__init__(
            seq=[jpl_lp('earth'), jpl_lp('venus'), jpl_lp('earth')],
            t0=t0,
            tof=tof,
            vinf=[0.5, 2.5],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding=tof_encoding,
            multi_objective=False
        )

    def get_name(self):
        return "Earth-Venus-Earth mga-1dsm"

    def get_extra_info(self):
        retval = "\ttof_encoding: " + self._tof_encoding
        return retval


# Problem P7: Earth-Venus-Earth MGA1DSM, single objective , direct encoding
eve_mga1dsm = _eve_mga1dsm_udp(tof_encoding='direct', t0=[
                                epoch(0), epoch(3000)], tof=[[10, 500], [10, 500]])
# Problem P8: Earth-Venus-Earth MGA1DSM, single objective , alpha encoding
eve_mga1dsm_a = _eve_mga1dsm_udp(tof_encoding='alpha', t0=[
                                epoch(0), epoch(3000)], tof=[300, 700])
# Problem P9: Earth-Venus-Earth MGA1DSM, single objective , eta encoding
eve_mga1dsm_n = _eve_mga1dsm_udp(
    tof_encoding='eta', t0=[epoch(0), epoch(3000)], tof=700)


