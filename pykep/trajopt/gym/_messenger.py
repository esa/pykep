from pykep.trajopt import mga_1dsm
from pykep.planet import jpl_lp, keplerian
from pykep import AU, DEG2RAD, MU_SUN, epoch

# MESSENGER (FULL) (we need to modify the safe radius of the planets to match the wanted problem)
_earth = jpl_lp('earth')
_venus = jpl_lp('venus')
_venus.safe_radius = 1.1
_mercury = jpl_lp('mercury')
_mercury.safe_radius = 1.05
_seq_messenger = [_earth,
                _venus,
                _venus,
                _mercury,
                _mercury,
                _mercury,
                _mercury
]

class _messenger_udp(mga_1dsm):
    """
    This class represents a rendezvous mission to Mercury modelled as an MGA-1DSM transfer. The selected fly-by sequence,
    E-VVMeMeMe-Me, and other parameters are inspired to the Messenger mission. We have only omitted the first Earth fly-by that
    was used to correct for launcher performances, since we here do not make use of a launcher model.
    As far as chemical propelled interplanetary trajectories go, this particular one is particularly complex and difficult
    to design. The time of flights among successive Mercury fly-bys allow for multiple rvolutions and resonances, making
    optimization techniques struggle to find the correct combination.
    The amount of specialistic knowledge that needs to be used to obtain a successfull design is significant.
    Finding a global optimization approach able to find a good trajectory in complete autonomy without making
    use of additional problem knowledge is possible, but limiting the number of fitness call is difficult.

    .. note::

       A significantly similar version of this problem was part of the no longer maintained GTOP database, 
       https://www.esa.int/gsp/ACT/projects/gtop/gtop.html. The exact definition is, though, different and results
       cannot thus not be compared to those posted in GTOP.
    """
    def __init__(self):
        """
        Write Me
        """
        super().__init__(
            seq =_seq_messenger,
            t0 = [1900, 2300],
            tof = [[100, 500], [100, 500], [100, 500], [100, 500], [100, 500], [100, 600]],
            vinf = [2.5, 4.05],
            add_vinf_dep = False,
            add_vinf_arr = True,
            tof_encoding = "direct",
            multi_objective = False,
            orbit_insertion = True,
            e_target = 0.704,
            rp_target = 2640000,
            eta_lb = 0.01,
            eta_ub = 0.99,
            rp_ub = 6.
        )

# Problem P15: Messenger mission MGA1DSM, single objective, direct encoding
messenger = _messenger_udp()

