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
:class:`pykep.trajopt.mga`                      class           A Multiple Gravity Assist Trajectory with no deep space manouvres
:class:`pykep.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space manouvre at each leg
:class:`pykep.trajopt.pl2pl_N_impulses`         class           A single leg transfer with N impulses
:class:`pykep.trajopt.lt_margo`                 class           A cubesat mission to near Earth asteroids. Solar Electric Propulsion and Earth gravity are modelled.
:class:`pykep.trajopt.mr_lt_nep`                class           A multiple randezvous low-thrust optimization problem (e.g. for asteroids in the main belt)
:class:`pykep.trajopt.direct_pl2pl`             class           A low-thrust transfer between planets using a direct transcription.
:class:`pykep.trajopt.indirect_pt2pt`           class           A low-thrust transfer between Cartesian states using an indirect transcription.
:class:`pykep.trajopt.indirect_or2or`           class           A low-thrust transfer between orbits using an indirect transcription.
:class:`pykep.trajopt._launchers`               class           Contains functors delivering the performance of launchers in terms of mass delivered
=========================================       =========       ================================================

The space trajectory gym
^^^^^^^^^^^^^^^^^^^^^^^^

Some instances of the classes above are provided in an Interplanetary Trajectory Optimization Gym. A problem test set 
we use to test and tune optimization approaches able to tackle a wide variety of problems with little / no tuning.

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
:class:`pykep.trajopt.gym.eve_mga1dsm_n`        instance        Earth-Venus-Earth transfer. Time of flights use eta encoding.
:class:`pykep.trajopt.gym.cassini2`             instance        Cassini inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.rosetta`              instance        Rosetta inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.tandem`               class           TandEM inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.juice`                instance        JUICE inspired MGA-1DSM problem.
:class:`pykep.trajopt.gym.juice_mo`             instance        JUICE inspired MGA-1DSM multiobjective problem.
:class:`pykep.trajopt.gym.messenger`            instance        Messenger inspired MGA-1DSM problem.
=========================================       =========       ================================================


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

------------

.. autoclass:: pykep.trajopt.mga_1dsm(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.__init__(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.pretty(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.plot(*args)

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