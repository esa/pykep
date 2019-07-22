from pykep.trajopt import mga
from pykep.planet import jpl_lp
from pykep import AU, DEG2RAD, MU_SUN, epoch

# CASSINI1
_seq_cassini1 = [jpl_lp('earth'), 
                jpl_lp('venus'), 
                jpl_lp('venus'), 
                jpl_lp('earth'), 
                jpl_lp('jupiter'), 
                jpl_lp('saturn')]

class _cassini1_udp(mga):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
            seq=_seq_cassini1,
            t0=[-1000., 0.],
            tof=[[30, 400], [100, 470], [30, 400], [400, 2000], [1000, 6000]],
            vinf=3.,
            tof_encoding='direct',
            orbit_insertion=True,
            e_target=0.98,
            rp_target=108950000)

    def get_name(self):
        return "Cassini MGA direct tof encoding (Trajectory Optimisation Gym P1)"

    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P1): Cassini MGA, single objective, direct encoding\n"
        retval += "\tPlanetary sequence" + \
            str([pl.name for pl in _seq_cassini1])
        return retval

    def __repr__(self):
        return self.get_name()


class _cassini1a_udp(mga):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
            seq=_seq_cassini1,
            t0=[-1000., 0.],
            tof=[4000., 7000.],
            vinf=3.,
            tof_encoding='alpha',
            orbit_insertion=True,
            e_target=0.98,
            rp_target=108950000)

    def get_name(self):
        return "Cassini MGA alpha tof encoding (Trajectory Optimisation Gym P2)"

    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P2): Cassini MGA, single objective, alpha encoding\n"
        retval += "\tPlanetary sequence" + \
            str([pl.name for pl in _seq_cassini1])
        return retval

    def __repr__(self):
        return self.get_name()


class _cassini1n_udp(mga):
    """
    Write Me
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
            seq=_seq_cassini1,
            t0=[-1000., 0.],
            tof=7000.,
            vinf=3.,
            tof_encoding='eta',
            orbit_insertion=True,
            e_target=0.98,
            rp_target=108950000)

    def get_name(self):
        return "Cassini1 MGA eta tof encoding (Trajectory Optimisation Gym P3)"

    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P3): Cassini MGA, single objective, eta encoding\n"
        retval += "\tPlanetary sequence" + \
            str([pl.name for pl in _seq_cassini1])
        return retval

    def __repr__(self):
        return self.get_name()

# Problem P1: Cassini MGA, single objective, direct encoding
cassini1 = _cassini1_udp()
# Problem P2: Cassini MGA, single objective, alpha encoding
cassini1_a = _cassini1a_udp()
# Problem P3: Cassini MGA, single objective, eta encoding
cassini1_n = _cassini1n_udp()