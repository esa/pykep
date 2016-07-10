.. _simsflanagan:

=======================
The simsflanagan module
=======================

The sims_flanagan module contains our implementation of the low-thrust interplanetary trajectory optimization model
first described in a paper by Sims and Flanagan in 1997 and later improved to include high-fidelity propagations
and easy mesh refinement by the PyKEP development team. The key papers containing most details are:

- Sims, Jon A., and Steve N. Flanagan. "Preliminary design of low-thrust interplanetary missions." (1997).

- Izzo, Dario. "PyGMO and PyKEP: open source tools for massively parallel optimization in astrodynamics (the case of interplanetary trajectory optimization)." Proceed. Fifth International Conf. Astrodynam. Tools and Techniques, ICATT. 2012.

The list of classes and the detailed documentation follows.

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`PyKEP.sims_flanagan.spacecraft`         class           represents a nuclear electric propelled spacecraft
:class:`PyKEP.sims_flanagan.sc_state`           class           represent the spacecraft state (r,v,m)
:class:`PyKEP.sims_flanagan.leg`                class           represents one leg in the Sims-Flanagan model
=========================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.sims_flanagan.spacecraft(*args)

  .. automethod:: PyKEP.sims_flanagan.spacecraft.__init__(*args)

------------

.. autoclass:: PyKEP.sims_flanagan.sc_state(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.__init__(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.set(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.get(*args)

------------

.. autoclass:: PyKEP.sims_flanagan.leg(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.__init__(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.high_fidelity(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.mismatch_constraints(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.throttles_constraints(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set_mu(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set_spacecraft(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.get_states()