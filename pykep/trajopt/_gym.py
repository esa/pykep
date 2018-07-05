from pykep.trajopt import mga, pl2pl_N_impulses
from pykep import planet

# Some hidden variables
_seq_cassini = [planet.jpl_lp('earth'), planet.jpl_lp('venus'), planet.jpl_lp(
'venus'), planet.jpl_lp('earth'), planet.jpl_lp('jupiter'), planet.jpl_lp('saturn')]
class _cassini1_udp(mga):
    def __init__(self):
        super(_cassini1_udp, self).__init__(
            seq=_seq_cassini, 
            t0=[-1000., 0.], 
            tof=[[30,400], [100,470], [30, 400], [400, 2000], [1000, 6000]], 
            vinf=3., 
            tof_encoding='direct',
            orbit_insertion=True,
            e_target = 0.98,
            rp_target = 108950000) 
    def get_name(self):
        return "Cassini MGA direct tof encoding (Trajectory Optimisation Gym P1)"
    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P1): Cassini MGA, single objective, direct encoding\n"
        retval+= "\tPlanetary sequence" + str([pl.name for pl in _seq_cassini])
        return retval
    def __repr__(self):
        return self.get_name()

class _cassini1a_udp(mga):
    def __init__(self):
        super(_cassini1a_udp, self).__init__(
            seq=_seq_cassini,
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
        retval+= "\tPlanetary sequence" + str([pl.name for pl in _seq_cassini])
        return retval
    def __repr__(self):
        return self.get_name()

class _cassini1n_udp(mga):
    def __init__(self):
        super(_cassini1n_udp, self).__init__(
            seq=_seq_cassini,
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
        retval+= "\tPlanetary sequence" + str([pl.name for pl in _seq_cassini])
        return retval
    def __repr__(self):
        return self.get_name()

class _emNimp_udp(pl2pl_N_impulses):
    def __init__(self, N = 3):
        super(_emNimp_udp, self).__init__(
                start=planet.jpl_lp('earth'),
                target=planet.jpl_lp('mars'),
                N_max=N,
                tof=[200., 700.],
                vinf=[0., 4.],
                phase_free=False,
                multi_objective=False,
                t0=[10000, 11000]
                )
        self.N = N
    def get_name(self):
        return "Earth-Mars " + str(self.N) + " impulses (Trajectory Optimisation Gym P" + str(int(4 + (self.N - 3.)/2.)) + ")"
    def get_extra_info(self):
        retval = "\tTrajectory Optimisation Gym problem (P3): Earth-Mars, " + str(int(4 + (self.N - 3.)/2.)) + " impulses, single objective\n"
        return retval
    def __repr__(self):
        return self.get_name()

class gym:
    def __init__(self):
        pass
    # Problem P1: Cassini MGA, single objective, direct encoding
    cassini1 = _cassini1_udp()
    # Problem P2: Cassini MGA, single objective, alpha encoding
    cassini1a =_cassini1a_udp()
    # Problem P3: Cassini MGA, single objective, eta encoding
    cassini1n =_cassini1n_udp()
    # Problem P4: Earth-Mars, 3 impulses, single objective, alpha encoding
    em3imp = _emNimp_udp(N = 3)
    # Problem P5: Earth-Mars, 5 impulses, single objective, alpha encoding
    em5imp = _emNimp_udp(N = 5)
    # Problem P6: Earth-Mars, 7 impulses, single objective, alpha encoding
    em7imp = _emNimp_udp(N = 7)
