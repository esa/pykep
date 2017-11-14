.. _trajopt:

===================================
The trajopt module (requires pygmo)
===================================

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`pykep.trajopt.lt_margo`                 class           A cubesat mission to near earth asteroids. Solar Electric Propulsion and Earth gravity can be modelled.
:class:`pykep.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space manouvre
:class:`pykep.trajopt.mr_lt_nep`                class           A multiple randezvous low-thrust optimization problem (e.g. for asteroids in the main belt)
:class:`pykep.trajopt.direct_pl2pl`             class           A transfer between planets using a direct transcription.
:class:`pykep.trajopt.indirect_pt2pt`           class           A transfer between Cartesian states using an indirect transcription.
:class:`pykep.trajopt.indirect_or2or`           class           A transfer between orbits using an indirect transcription.
:class:`pykep.trajopt.indirect_pt2or`           class           A transfer transfer between a Cartesian state and orbit using an indirect transcription.
=========================================       =========       ================================================

Detailed Documentation
======================

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

.. autoclass:: pykep.trajopt.mga_1dsm(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.__init__(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.pretty(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.plot(*args)

------------

.. autoclass:: pykep.trajopt.mr_lt_nep(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.__init__(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.plot(*args)

------------

.. autoclass:: pykep.trajopt.indirect_pt2pt(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2pt.__init__(x0, xf, mass, thrust, isp, mu, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, atol=1e-10, rtol=1e-10)

------------

.. autoclass:: pykep.trajopt.indirect_or2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_or2or.__init__(elem0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)


------------

.. autoclass:: pykep.trajopt.indirect_pt2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2or.__init__(x0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)
