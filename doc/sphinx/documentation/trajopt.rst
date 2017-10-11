.. _trajopt:

===================================
The trajopt module (requires PyGMO)
===================================

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`PyKEP.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space manouvre
:class:`PyKEP.trajopt.pl2pl_N_impulses`         class           A direct transfer between planets allowing for many DSMs
:class:`PyKEP.trajopt.mga_lt_nep`               class           A multiple Gravity Assist Trajectory low-thrust optimization problem
:class:`PyKEP.trajopt.mr_lt_nep`                class           A multiple randezvous low-thrust optimization problem
:class:`PyKEP.trajopt.indirect_pt2pt`           class           A direct transfer between Cartesian states in the indirect transcription.
:class:`PyKEP.trajopt.indirect_pl2pl`           class           A direct transfer between planets in the indirect transcription.
:class:`PyKEP.trajopt.indirect_or2or`           class           A direct transfer between orbits in the indirect transcription.
:class:`PyKEP.trajopt.indirect_pt2or`           class           A direct transfer between a Cartesian state and orbit in the indirect transcription.
=========================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.trajopt.mga_1dsm(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.__init__(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.set_tof(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.set_launch_window(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.set_vinf(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.pretty(*args)

   .. automethod:: PyKEP.trajopt.mga_1dsm.plot(*args)

------------

.. autoclass:: PyKEP.trajopt.mga_lt_nep(*args)

   .. automethod:: PyKEP.trajopt.mga_lt_nep.__init__(*args)

   .. automethod:: PyKEP.trajopt.mga_lt_nep.high_fidelity(*args)

   .. automethod:: PyKEP.trajopt.mga_lt_nep.ic_from_mga_1dsm(*args)

   .. automethod:: PyKEP.trajopt.mga_lt_nep.plot(*args)

------------

.. autoclass:: PyKEP.trajopt.mr_lt_nep(*args)

   .. automethod:: PyKEP.trajopt.mr_lt_nep.__init__(*args)

   .. automethod:: PyKEP.trajopt.mr_lt_nep.plot(*args)

------------

.. autoclass:: PyKEP.trajopt.pl2pl_N_impulses(*args)

   .. automethod:: PyKEP.trajopt.pl2pl_N_impulses.__init__(*args)

------------

.. autoclass:: PyKEP.trajopt.indirect_pt2pt(*args)

  .. automethod:: PyKEP.trajopt.indirect_pt2pt.__init__(x0, xf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)

  .. automethod:: PyKEP.trajopt.indirect_pt2pt.plot_traj(z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU)

------------

.. autoclass:: PyKEP.trajopt.indirect_pl2pl(*args)

  .. automethod:: PyKEP.trajopt.indirect_pl2pl.__init__(p0, pf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)

  .. automethod:: PyKEP.trajopt.indirect_pl2pl.plot_traj(z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU)

------------

.. autoclass:: PyKEP.trajopt.indirect_or2or(*args)

  .. automethod:: PyKEP.trajopt.indirect_or2or.__init__(elem0, elemf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)

  .. automethod:: PyKEP.trajopt.indirect_pl2pl.plot_traj(z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU)

------------

.. autoclass:: PyKEP.trajopt.indirect_pt2or(*args)

  .. automethod:: PyKEP.trajopt.indirect_pt2or.__init__(elemf, elemf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)

  .. automethod:: PyKEP.trajopt.indirect_pt2or.plot_traj(z, mark="k.-", atol=1e-12, rtol=1e-12, units=pk.AU)
