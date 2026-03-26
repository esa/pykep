.. _gym:

Trajectory Optimization Gym
###########################

A number of interplanetary trajectory problems are provided in `pykep` in the form of 
User Defined Problems (UDP) compatible with the `pygmo <https://esa.github.io/pygmo2/>`_ :cite:p:`pagmo` python package.
All of the problems are instantiated upon import of the `pykep` module and can be used directly. The collection of all
problems is called the "pykep gym" and hopes to become an established benchmark set to test the performances of evolutionary
and traditional optimisation techniques on trajectory design problems.

.. currentmodule:: pykep.trajopt.gym
 
MGA problems
************************

.. autoattribute:: pykep.trajopt.gym.cassini1

This is an MGA problem inspired to the Cassini spacecraft interplanetary transfer to Saturn. 
The objective of this mission is to reach Saturn and to be captured by its gravity into an orbit having pericenter
radius :math:`r_p=108950` km, and eccentricity :math:`e=0.98`. The planetary fly-by sequence considered
is E-VVEJ-S (as the one used by Cassini spacecraft). As objective function we use the total :math:`\Delta V`
accumulated during the mission, including the launch :math:`\Delta V` (exceeding a value provided for free by the launcher)
and the various :math:`\Delta V` one needs to give at the planets and upon arrival to perform the final orbit injection. 

.. autoattribute:: pykep.trajopt.gym.cassini1_a

This is the same MGA problem as :class:`~pykep.trajopt.gym.cassini1`, but with using
a different encoding of the decision variables. In this case, the decision variables are no longer directly
the times of flight but an equivalent set called the alpha encoding. This encoding is particularly useful
when bounds on the times of flight are not known a priori. It creates a less evolvable problem, albeit a less biased one,
or at least one that is not using any apriori knowledge on the final solution. The decision vector is also augmented with one
more variable: the total time of flight which is thus controlled via its bounds, making the encoding particularly suitable
for multi-objective optimization problems where the total time of flight is a second objective for which we want to search solutions
on the pareto front.

.. autoattribute:: pykep.trajopt.gym.cassini1_n

This is the same MGA problem as :class:`~pykep.trajopt.gym.cassini1`, but with using
a different encoding of the decision variables. In this case, the decision variables are no longer directly
the times of flight but an equivalent set called the eta encoding. This encoding is particularly useful
when bounds on the times of flight are not known a priori. It creates a less evolvable problem, albeit a less biased one.

MGA-1DSM problems
************************

.. autoattribute:: pykep.trajopt.gym.cassini2

This is an MGA-1DSM problem inspired to the Cassini spacecraft interplanetary transfer to Saturn. 
It is fundamentally similar to :class:`~pykep.trajopt.gym.cassini1`, but it allows for one deep space maneuver (DSM) at any point
of each trajectory leg, and does not account for the launch :math:`\Delta V` which is bounded instead to a maximum value.

.. autoattribute:: pykep.trajopt.gym.rosetta

This is an MGA-1DSM problem inspired to the Rosetta spacecraft interplanetary transfer to comet 67P/Churyumov-Gerasimenko.
The objective of this mission is to rendezvous with the comet matching its heliocentric position and velocity exactly.
The planetary fly-by sequence considered is E-EMEE-comet (as the one used by Rosetta spacecraft). 
As objective function we use the total :math:`\Delta V` computed as the sum of all the :math:`\Delta V` needed to perform the various
maneuvers and fly-bys. The launch :math:`\Delta V` is not accounted for, as it is bounded to a maximum value.

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm

This is an MGA-1DSM problem useful to test the performance of evolutionary algorithms on a simpe trajectory problem.
The planetary encounter sequence is E-V-E, and the objective is to rendezvous with the Earth matching its heliocentric
position and velocity exactly. As objective function we use the total :math:`\Delta V` computed as the sum of all the :math:`\Delta V` needed to perform the various
maneuvers and fly-bys. The launch :math:`\Delta V` is not accounted for, as it is bounded to a maximum value.

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_a

This is the same MGA-1DSM problem as :class:`~pykep.trajopt.gym.eve_mga_1dsm`, but with using
a different encoding of the decision variables. In this case, the decision variables are no longer directly
the times of flight but an equivalent set called the alpha encoding (see :class:`~pykep.trajopt.gym.cassini1_a`).

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_n

This is the same MGA-1DSM problem as :class:`~pykep.trajopt.gym.eve_mga_1dsm`, but with using
a different encoding of the decision variables. In this case, the decision variables are no longer directly
the times of flight but an equivalent set called the eta encoding (see :class:`~pykep.trajopt.gym.cassini1_n`).

.. autoattribute:: pykep.trajopt.gym.juice

This is an MGA-1DSM problem inspired to the JUICE spacecraft interplanetary transfer to Jupiter. 
The selected fly-by sequence, E-EVEME-J, and other parameters are inspired by the ESA Juice mission. 
A launcher model is included, namely an Ariane5 launch from Kourou. 
JUICE - JUpiter ICy moons Explorer - is the first large-class mission in ESA's Cosmic Vision 2015-2025 programme.
Launched on the 14th of April 2023, ESA’s Jupiter Icy Moons Explorer, Juice, will make detailed observations
of the giant gas planet and its three large ocean-bearing moons – Ganymede, Callisto and Europa – with a suite of
remote sensing, geophysical and in situ instruments. The mission will characterise these moons as both
planetary objects and possible habitats, explore Jupiter’s complex environment in depth, and study the wider
Jupiter system as an archetype for gas giants across the Universe.

.. autoattribute:: pykep.trajopt.gym.juice_mo

This is the same MGA-1DSM problem as :class:`~pykep.trajopt.gym.juice`, but with using
the alpha encoding for the time of flight variables (see :class:`~pykep.trajopt.gym.cassini1_a`).
The problem is also multi-objective, with the total time of flight as a second objective, leveraging the 
alpha encoding to control the total time of flight (at the cost of evolvability).

.. autoattribute:: pykep.trajopt.gym.messenger

This an MGA-1DSM problem inspired by the messenger mission that orbited the planet Mercury between 2011 and 2015, studying Mercury's chemical composition,
geology, and magnetic field. The name is a backronym for Mercury Surface, Space Environment, Geochemistry, and Ranging, and a reference to the messenger
god Mercury from Roman mythology. The selected fly-by sequence, E-VVMeMeMe-Me, and other parameters are coincident to the actual Messenger mission.
We have only omitted the first Earth fly-by that was used to correct for launcher performances, since we here do not make use of a launcher model.
As far as chemical propelled interplanetary trajectories go, this particular one is particularly complex and difficult to design. 
The time of flights among successive Mercury fly-bys allow for multiple rvolutions and resonances, making optimization techniques struggle to find the correct combination.
The amount of specialistic knowledge that needs to be used to obtain a successful design is significant.
Finding a global optimization approach able to find a good trajectory in complete autonomy without making
use of additional problem knowledge is challenging but possible, but limiting the number of fitness call is difficult.

Multiple impulse problem
************************

.. autoattribute:: pykep.trajopt.gym.em3imp

.. autoattribute:: pykep.trajopt.gym.em5imp

.. autoattribute:: pykep.trajopt.gym.em7imp

