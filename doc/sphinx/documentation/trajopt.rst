.. _trajopt:

===================================
The trajopt module (requires pygmo)
===================================

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
