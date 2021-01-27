.. _trajopt:

===================================
The trajopt module
===================================

This module contains classes that help instantiating interplanetary trajectory optimization problems. While
the resulting problems are compatible with the optimization package pygmo, which can thus be used to assemble a 
solution strategy, they can be also interfaced with any other optimization package / method.

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`pykep.trajopt.mga`                      class           A Multiple Gravity Assist Trajectory with no deep space maneuvers
:class:`pykep.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space maneuver at each leg
:class:`pykep.trajopt.pl2pl_N_impulses`         class           A single leg transfer with N impulses
:class:`pykep.trajopt.lt_margo`                 class           A cubesat mission to near Earth asteroids. Solar Electric Propulsion and Earth gravity are modelled.
:class:`pykep.trajopt.mr_lt_nep`                class           A multiple rendezvous low-thrust optimization problem (e.g. for asteroids in the main belt)
:class:`pykep.trajopt.direct_pl2pl`             class           A low-thrust transfer between planets using a direct transcription.
:class:`pykep.trajopt.indirect_pt2pt`           class           A low-thrust transfer between Cartesian states using an indirect transcription.
:class:`pykep.trajopt.indirect_or2or`           class           A low-thrust transfer between orbits using an indirect transcription.
:class:`pykep.trajopt._launchers`               class           Contains functors delivering the performance of launchers in terms of mass delivered
=========================================       =========       ================================================

The space trajectory gym
^^^^^^^^^^^^^^^^^^^^^^^^

Some instances of the classes above are provided in an Interplanetary Trajectory Optimization Gym. A problem test set 
we use to test and tune optimization approaches able to tackle a wide variety of problems with little / no tuning.

.. note::

   Most problem have a `pretty` and a `plot` method that can be used to visualize decision vectors (chromosomes).

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`pykep.trajopt.gym.cassini1`             instance        Cassini inspired MGA problem. Time of flights are encoded directly in the chromosome.
:class:`pykep.trajopt.gym.cassini1_a`           instance        Cassini inspired MGA problem. Time of flights use alpha encoding.
:class:`pykep.trajopt.gym.cassini1_n`           instance        Cassini inspired MGA problem. Time of flights use eta encoding.
:class:`pykep.trajopt.gym.em3imp`               instance        Earth-Mars transfer with 3 impulses
:class:`pykep.trajopt.gym.em5imp`               instance        Earth-Mars transfer with 5 impulses
:class:`pykep.trajopt.gym.em7imp`               instance        Earth-Mars transfer with 7 impulses
:class:`pykep.trajopt.gym.eve_mga1dsm`          instance        Earth-Venus-Earth transfer. Time of flights are encoded directly in the chromosome.
:class:`pykep.trajopt.gym.eve_mga1dsm_a`        instance        Earth-Venus-Earth transfer. Time of flights use alpha encoding.
:class:`pykep.trajopt.gym.eve_mga1dsm_n`        instance        Earth-Venus-Earth transfer. Time of flights use eta encoding.
:class:`pykep.trajopt.gym.cassini2`             instance        Cassini inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.rosetta`              instance        Rosetta inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.tandem`               class           TandEM inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.juice`                instance        JUICE inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.juice_mo`             instance        JUICE inspired MGA-1DSM multiobjective problem.
:class:`pykep.trajopt.gym.messenger`            instance        Messenger inspired MGA-1DSM problem.
=========================================       =========       ================================================

---------------------------------------------------------------------------------

Detailed Documentation
======================

.. autoclass:: pykep.trajopt.mga(*args)

   .. automethod:: pykep.trajopt.mga.__init__(*args)

   .. automethod:: pykep.trajopt.mga.alpha2direct(*args)

   .. automethod:: pykep.trajopt.mga.direct2alpha(*args)

   .. automethod:: pykep.trajopt.mga.eta2direct(*args)

   .. automethod:: pykep.trajopt.mga.direct2eta(*args)

   .. automethod:: pykep.trajopt.mga.pretty(*args)

   .. automethod:: pykep.trajopt.mga.plot(*args)

   .. automethod:: pykep.trajopt.mga.get_eph_function(*args)

   .. automethod:: pykep.trajopt.mga.plot_distance_and_flybys(*args)

------------

.. autoclass:: pykep.trajopt.mga_1dsm(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.__init__(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.pretty(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.plot(*args)

   .. automethod:: pykep.trajopt.mga.get_eph_function(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.plot_distance_and_flybys(*args)

------------

.. autoclass:: pykep.trajopt.pl2pl_N_impulses(*args)

   .. automethod:: pykep.trajopt.pl2pl_N_impulses.__init__(*args)

   .. automethod:: pykep.trajopt.pl2pl_N_impulses.pretty(*args)

   .. automethod:: pykep.trajopt.pl2pl_N_impulses.plot(*args)

-------------

.. autoclass:: pykep.trajopt.lt_margo(*args)

   .. automethod:: pykep.trajopt.lt_margo.__init__(*args)

   .. automethod:: pykep.trajopt.lt_margo.pretty(*args)

   .. automethod:: pykep.trajopt.lt_margo.plot_traj(*args)

   .. automethod:: pykep.trajopt.lt_margo.plot_dists_thrust(*args)

   .. automethod:: pykep.trajopt.lt_margo.double_segments(*args)

------------

.. autoclass:: pykep.trajopt.direct_pl2pl(*args)

   .. automethod:: pykep.trajopt.direct_pl2pl.__init__(*args)

   .. automethod:: pykep.trajopt.direct_pl2pl.pretty(*args)

   .. automethod:: pykep.trajopt.direct_pl2pl.plot_traj(*args)

   .. automethod:: pykep.trajopt.direct_pl2pl.plot_control(*args)

   .. automethod:: pykep.trajopt.direct_pl2pl.get_traj(*args)

------------

.. autoclass:: pykep.trajopt.mr_lt_nep(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.__init__(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.plot(*args)

------------

.. autoclass:: pykep.trajopt.indirect_pt2pt(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2pt.__init__(*args)

------------

.. autoclass:: pykep.trajopt.indirect_or2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_or2or.__init__(*args)


------------

.. autoclass:: pykep.trajopt.indirect_pt2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2or.__init__(*args)

------------

.. autoclass:: pykep.trajopt._launchers

  .. automethod:: pykep.trajopt._launchers.atlas501(*args)

  .. automethod:: pykep.trajopt._launchers.soyuzf(*args)

  .. automethod:: pykep.trajopt._launchers.ariane5(*args)

------------

.. autoattribute:: pykep.trajopt.gym.cassini1

This is an MGA problem inspired to the Cassini spacecraft interplanetary transfer to Saturn. 
The objective of this mission is to reach Saturn and to be captured by its gravity into an orbit having
pericenter radius r_p=108950 km, and eccentricity e=0.98. The planetary fly-by sequence considered is
E-VVEJ-S (as the one used by Cassini spacecraft). 
As objective function we use the total deltaV accumulated during the mission, including the launch deltaV
and the various deltaV one needs to give at the planets and upon arrival to perform the final orbit injection.

.. note::

   A more complex representation of this interplanetary trajectory transfer can be found
   in :class:`pykep.trajopt.gym.cassini2`

------------

.. autoattribute:: pykep.trajopt.gym.cassini1_a

This is the same as :class:`pykep.trajopt.gym.cassini1`, 
but the time of flights are encoded using the alpha encoding technique.

------------

.. autoattribute:: pykep.trajopt.gym.cassini1_n

This is the same as :class:`pykep.trajopt.gym.cassini1`, 
but the time of flights are encoded using the eta encoding technique.

------------

.. autoattribute:: pykep.trajopt.gym.em3imp

This is a simple Earth Mars transfer that has the possibility to use multiple impulses. While the actual solution
may not need more than a few, the problem is used to study the alpha encoding for maneuver timings and its evolvability.

.. note::

   The problem is rather simple and is used to study the multiple impulse transcription (using alpha encoding)
   comparing it to the folowing problems :class:`pykep.trajopt.gym.em5imp` and :class:`pykep.trajopt.gym.em7imp`
   having increasing number of impulses.

------------

.. autoattribute:: pykep.trajopt.gym.em5imp

This is the same as :class:`pykep.trajopt.gym.em3imp`, 
but the number of impulses is increased.

------------

.. autoattribute:: pykep.trajopt.gym.em7imp

This is the same as :class:`pykep.trajopt.gym.em3imp`, 
but the number of impulses is increased.

------------

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm

This is an Earth - Venus - Earth transfer where 1 deep space maneuver is allowed at each leg. 

.. note::

   Together with the problems :class:`pykep.trajopt.gym.eve_mga1dsm_a` and :class:`pykep.trajopt.gym.eve_mga1dsm_n`
   this simple problem is used to study the alpha, eta and direct encoding for the time of flights of an mga-1DSM
   type of problem.

------------

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_a

This is the same as :class:`pykep.trajopt.gym.eve_mga1dsm`, 
but the time of flights are encoded using the alpha encoding technique.

------------

.. autoattribute:: pykep.trajopt.gym.eve_mga1dsm_n

This is the same as :class:`pykep.trajopt.gym.eve_mga1dsm`, 
but the time of flights are encoded using the eta encoding technique.

------------

.. autoattribute:: pykep.trajopt.gym.cassini2

This is an MGA-1DSM problem inspired to the `Cassini <https://www.nasa.gov/mission_pages/cassini/main/index.html>`_ 
spacecraft interplanetary transfer to Saturn. The objective of this mission is to reach Saturn and to be
captured by its gravity into an orbit having pericenter radius r_p=108950 km, and eccentricity e=0.98. 
The planetary fly-by sequence considered is E-VVEJ-S (as the one used by Cassini spacecraft). 
As objective function we use the total deltaV accumulated during the mission by the spacecraft,
including the one necessary for the final orbit injection.

.. note::

   A similar problem was also part of the ESA's `GTOP database <https://www.esa.int/gsp/ACT/projects/gtop/gtop.html>`_ 
   with the same name, but different implementation details and mission definition. They should not be compared.

------------

.. autoattribute:: pykep.trajopt.gym.rosetta

    This represents a rendezvous mission to the comet 67P/Churyumov-Gerasimenko modelled as an MGA-1DSM transfer.
    The fly-by sequence selected (i.e. E-EMEE-C) is similar to the one planned for the spacecraft 
    `Rosetta <http://www.esa.int/Our_Activities/Space_Science/Rosetta>`_. The objective function considered is the
    total mission delta V. No launcher model is employed and a final rendezvous with the comet is included
    in the delta V computations.

.. note::

   A similar problem was also part of the ESA's `GTOP database <https://www.esa.int/gsp/ACT/projects/gtop/gtop.html>`_ 
   with the same name, but different implementation details and mission definition. They should not be compared.

------------

.. autoclass:: pykep.trajopt.gym.tandem

This class represents a rendezvous mission to Saturn modelled as an MGA-1DSM transfer. Mission parameters are
inspired to the TandEM mission. A launcher model (i.e. Atlas 501) is also used, so that the final mass delivered
at Saturn is the main objective of this optimization problem.

The problem draws inspiration from the work performed in April 2008 by the
European Space Agency working group on mission analysis on the mission named TandEM. TandEM is an interplanetary
mission aimed at reaching Titan and Enceladus (two moons of Saturn). 

.. note::

   The Titan and Enceladus Mission (TandEM), an ambitious scientific mission to study the Saturnian system
   with particular emphasis on the moons Titan and Enceladus, was selected in October 2007 as a candidate mission
   within the ESA Cosmic Vision plan. In February 2009, TandEM exited the Cosmic Vision programme when ESA
   and NASA chose EJSM-Laplace as the L-class outer Solar System mission candidate.

.. note::

   A significantly similar version of this problem was part of the no longer maintained GTOP database, 
   https://www.esa.int/gsp/ACT/projects/gtop/gtop.html. The exact definition is, though, different and results
   cannot thus not be compared to those posted in GTOP.

.. automethod:: pykep.trajopt.gym._tandem._tandem_udp.__init__(prob_id = 1, constrained = True)

-----------

.. autoattribute:: pykep.trajopt.gym.juice

This class represents a rendezvous mission to Jupiter modelled as an MGA-1DSM transfer. The selected fly-by sequence,
E-EVEME-J, and other parameters are inspired to the ESA `JUICE <http://sci.esa.int/juice/>`_ mission. 
A launcher model (i.e. Ariane 5) is also used, so that the final mass delivered at Saturn is the main objective
of this optimization problem.

.. note::

   JUICE - JUpiter ICy moons Explorer - is the first large-class mission in ESA's Cosmic Vision 2015-2025 programme. 
   Planned for launch in 2022 and arrival at Jupiter in 2029, it will spend at least three years making detailed
   observations of the giant gaseous planet Jupiter and three of its largest moons, Ganymede, Callisto and Europa.

-----------

.. autoattribute:: pykep.trajopt.gym.juice_mo

This is the multiobjective version of :class:`pykep.trajopt.gym.juice`. 
Time of flights are encoded using the alpha encoding.

------------

.. autoattribute:: pykep.trajopt.gym.messenger

    This class represents a rendezvous mission to Mercury modelled as an MGA-1DSM transfer. The selected fly-by sequence,
    E-VVMeMeMe-Me, and other parameters are inspired to the Messenger mission. We have only omitted the first Earth fly-by that
    was used to correct for launcher performances, since we here do not make use of a launcher model.
    As far as chemical propelled interplanetary trajectories go, this particular one is particularly complex and difficult
    to design. The time of flights among successive Mercury fly-bys allow for multiple revolutions and resonances, making
    optimization techniques struggle to find the correct combination.
    The amount of specialistic knowledge that needs to be used to obtain a successful design is significant.
    Finding a global optimization approach able to find a good trajectory in complete autonomy without making
    use of additional problem knowledge is possible, but limiting the number of fitness calls is difficult.

.. note::

   A similar problem was also part of the ESA's `GTOP database <https://www.esa.int/gsp/ACT/projects/gtop/gtop.html>`_ 
   with the same name, but different implementation details and mission definition. They should not be compared.
