.. _trajopt:

===================================
The trajopt module (requires PyGMO)
===================================

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`pykep.trajopt.mga_1dsm`                 class           A multiple Gravity Assist Trajectory with one deep space manouvre
:class:`pykep.trajopt.pl2pl_N_impulses`         class           A direct transfer between planets allowing for many DSMs
:class:`pykep.trajopt.mga_lt_nep`               class           A multiple Gravity Assist Trajectory low-thrust optimization problem
:class:`pykep.trajopt.mr_lt_nep`                class           A multiple randezvous low-thrust optimization problem
:class:`pykep.trajopt._direct_base`             class           Base direct transfer between Cartesian states in the direct transcription.
:class:`pykep.trajopt.direct_or2or`             class           A direct transfer between orbits in the direct transcription.
:class:`pykep.trajopt._indirect_base`           class           Base direct transfer between Cartesian states in the indirect transcription.
:class:`pykep.trajopt.indirect_pt2pt`           class           A direct transfer between Cartesian states in the indirect transcription.
:class:`pykep.trajopt.indirect_pl2pl`           class           A direct transfer between planets in the indirect transcription.
:class:`pykep.trajopt.indirect_or2or`           class           A direct transfer between orbits in the indirect transcription.
:class:`pykep.trajopt.indirect_pt2or`           class           A direct transfer between a Cartesian state and orbit in the indirect transcription.
=========================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: pykep.trajopt.mga_1dsm(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.__init__(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.set_tof(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.set_launch_window(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.set_vinf(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.pretty(*args)

   .. automethod:: pykep.trajopt.mga_1dsm.plot(*args)

------------

.. autoclass:: pykep.trajopt.mga_lt_nep(*args)

   .. automethod:: pykep.trajopt.mga_lt_nep.__init__(*args)

   .. automethod:: pykep.trajopt.mga_lt_nep.high_fidelity(*args)

   .. automethod:: pykep.trajopt.mga_lt_nep.ic_from_mga_1dsm(*args)

   .. automethod:: pykep.trajopt.mga_lt_nep.plot(*args)

------------

.. autoclass:: pykep.trajopt.mr_lt_nep(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.__init__(*args)

   .. automethod:: pykep.trajopt.mr_lt_nep.plot(*args)

------------

.. autoclass:: pykep.trajopt.pl2pl_N_impulses(*args)

   .. automethod:: pykep.trajopt.pl2pl_N_impulses.__init__(*args)

------------

.. autoclass:: pykep.trajopt._direct_base(object)

  .. automethod:: pykep.trajopt._direct_base.plot_traj(z, units=pykep.AU, N=20)
  .. automethod:: pykep.trajopt._direct_base.plot_control(z, mark="k.-", time=True)
  .. automethod:: pykep.trajopt._direct_base.get_traj(z)

------------

.. autoclass:: pykep.trajopt.direct_or2or(_direct_base)

  .. automethod:: pykep.trajopt.direct_or2or.__init__(elem0, elemf, mass, thrust, isp, nseg, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, mu=pykep.MU_SUN, hf=True)

------------

.. autoclass:: pykep.trajopt._indirect_base(object)

  .. automethod:: pykep.trajopt._indirect_base.plot_traj(z, mark="k.-", atol=1e-12, rtol=1e-12, units=pykep.AU)
  .. automethod:: pykep.trajopt._indirect_base.plot_control(z, mark="k.-", atol=1e-12, rtol=1e-12)

------------

.. autoclass:: pykep.trajopt.indirect_pt2pt(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2pt.__init__(x0, xf, mass, thrust, isp, mu, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, atol=1e-10, rtol=1e-10)

------------

.. autoclass:: pykep.trajopt.indirect_pl2pl(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pl2pl.__init__(p0, pf, mass, thrust, isp, atol, rtol, t0lb, t0ub, Tlb, Tub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)


------------

.. autoclass:: pykep.trajopt.indirect_or2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_or2or.__init__(elem0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, M0lb, M0ub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)


------------

.. autoclass:: pykep.trajopt.indirect_pt2or(_indirect_base)

  .. automethod:: pykep.trajopt.indirect_pt2or.__init__(x0, elemf, mass, thrust, isp, atol, rtol, Tlb, Tub, Mflb, Mfub, freemass=True, freetime=True, alpha=1, bound=True, mu=pk.MU_SUN)
