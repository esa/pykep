from pykep.trajopt import pl2pl_N_impulses
from pykep.planet import jpl_lp


class _emNimp_udp(pl2pl_N_impulses):
    """
    Write Me
    """
    def __init__(self, N=3):
        """
        Write Me
        """
        super().__init__(
            start=jpl_lp('earth'),
            target=jpl_lp('mars'),
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
        retval = "\tTrajectory Optimisation Gym problem (P"+ str(int(4 + (self.N - 3.)/2.)) +"): Earth-Mars, " + str(
            int(4 + (self.N - 3.)/2.)) + " impulses, single objective\n"
        return retval

    def __repr__(self):
        return self.get_name()


# Problem P4: Earth-Mars, 3 impulses, single objective, alpha encoding
em3imp = _emNimp_udp(N=3)
# Problem P5: Earth-Mars, 5 impulses, single objective, alpha encoding
em5imp = _emNimp_udp(N=5)
# Problem P6: Earth-Mars, 7 impulses, single objective, alpha encoding
em7imp = _emNimp_udp(N=7)
