.. _planet:

==================
The planet module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`PyKEP.planet._base`                              class           The base class for all planets (cannot be instantiated)
:class:`PyKEP.planet.keplerian`                          class           A simple planet with Keplerian ephemerides
:class:`PyKEP.planet.jpl_lp`                             class           A solar system planet using jpl low-precision ephemerides
:class:`PyKEP.planet.tle`                                class           An Earth artificial satellite from its TLE (ephemerides are computed via the SGP4 propagator)
:class:`PyKEP.planet.spice`                              class           A planet with ephemerides computed using the JPL SPICE Toolbox (requires BUILD_SPICE option active when building from cmake)
:class:`PyKEP.planet.mpcorb`                             class           A planet from the MPCORB database (keplerian ephemerides)
:class:`PyKEP.planet.gtoc2`                              class           An asteroid from the GTOC2 competition (keplerian ephemerides)
:class:`PyKEP.planet.gtoc5`                              class           An asteroid from the GTOC5 competition (keplerian ephemerides)
:class:`PyKEP.planet.gtoc6`                              class           A Jupiter moon from the GTOC6 competition (keplerian ephemerides)
:class:`PyKEP.planet.gtoc7`                              class           An asteroid from the GTOC7 competition (keplerian ephemerides)
==================================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.planet._base(*args)

  .. automethod:: PyKEP.planet._base.__init__(*args)

  .. automethod:: PyKEP.planet._base.eph(*args)

  .. automethod:: PyKEP.planet._base.osculating_elements(*args)

  .. automethod:: PyKEP.planet._base.compute_period(*args)

  .. automethod:: PyKEP.planet._base.human_readable_extra(*args)

  .. autoattribute:: PyKEP.planet._base.mu_central_body

  .. autoattribute:: PyKEP.planet._base.mu_self

  .. autoattribute:: PyKEP.planet._base.radius

  .. autoattribute:: PyKEP.planet._base.safe_radius

  .. autoattribute:: PyKEP.planet._base.name


------------

.. autoclass:: PyKEP.planet.keplerian(*args)

  .. automethod:: PyKEP.planet.keplerian.__init__(*args)

------------

.. autoclass:: PyKEP.planet.jpl_lp(*args)

  .. automethod:: PyKEP.planet.jpl_lp.__init__(*args)

------------

.. autoclass:: PyKEP.planet.mpcorb(*args)

  .. automethod:: PyKEP.planet.mpcorb.__init__(*args)

  .. autoattribute:: PyKEP.planet.mpcorb.H

  .. autoattribute:: PyKEP.planet.mpcorb.n_observations

  .. autoattribute:: PyKEP.planet.mpcorb.n_oppositions

  .. autoattribute:: PyKEP.planet.mpcorb.year_of_discovery

------------

.. autoclass:: PyKEP.planet.tle(*args)

  .. automethod:: PyKEP.planet.tle.__init__(*args)

------------

.. autoclass:: PyKEP.planet.spice(*args)

  .. automethod:: PyKEP.planet.spice.__init__(*args)

------------

.. autoclass:: PyKEP.planet.gtoc2(*args)

  .. automethod:: PyKEP.planet.gtoc2.__init__(*args)

------------

.. autoclass:: PyKEP.planet.gtoc5(*args)

  .. automethod:: PyKEP.planet.gtoc5.__init__(*args)

------------

.. autoclass:: PyKEP.planet.gtoc6(*args)

  .. automethod:: PyKEP.planet.gtoc6.__init__(*args)

------------

.. autoclass:: PyKEP.planet.gtoc7(*args)

  .. automethod:: PyKEP.planet.gtoc7.__init__(*args)

