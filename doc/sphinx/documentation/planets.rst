.. _planet:

==================
The planet module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`pykep.planet._base`                              class           The base class for all planets (cannot be instantiated)
:class:`pykep.planet.keplerian`                          class           A simple planet with Keplerian ephemerides
:class:`pykep.planet.jpl_lp`                             class           A solar system planet using jpl low-precision ephemerides
:class:`pykep.planet.tle`                                class           An Earth artificial satellite from its TLE (ephemerides are computed via the SGP4 propagator)
:class:`pykep.planet.spice`                              class           A planet with ephemerides computed using the JPL SPICE Toolbox (requires BUILD_SPICE option active when building from cmake)
:class:`pykep.planet.mpcorb`                             class           A planet from the MPCORB database (keplerian ephemerides)
:class:`pykep.planet.gtoc2`                              class           An asteroid from the GTOC2 competition (keplerian ephemerides)
:class:`pykep.planet.gtoc5`                              class           An asteroid from the GTOC5 competition (keplerian ephemerides)
:class:`pykep.planet.gtoc6`                              class           A Jupiter moon from the GTOC6 competition (keplerian ephemerides)
:class:`pykep.planet.gtoc7`                              class           An asteroid from the GTOC7 competition (keplerian ephemerides)
==================================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: pykep.planet._base(*args)

  .. automethod:: pykep.planet._base.__init__(*args)

  .. automethod:: pykep.planet._base.eph(*args)

  .. automethod:: pykep.planet._base.osculating_elements(*args)

  .. automethod:: pykep.planet._base.compute_period(*args)

  .. automethod:: pykep.planet._base.human_readable_extra(*args)

  .. autoattribute:: pykep.planet._base.mu_central_body

  .. autoattribute:: pykep.planet._base.mu_self

  .. autoattribute:: pykep.planet._base.radius

  .. autoattribute:: pykep.planet._base.safe_radius

  .. autoattribute:: pykep.planet._base.name


------------

.. autoclass:: pykep.planet.keplerian(*args)

  .. automethod:: pykep.planet.keplerian.__init__(*args)

------------

.. autoclass:: pykep.planet.jpl_lp(*args)

  .. automethod:: pykep.planet.jpl_lp.__init__(*args)

------------

.. autoclass:: pykep.planet.mpcorb(*args)

  .. automethod:: pykep.planet.mpcorb.__init__(*args)

  .. autoattribute:: pykep.planet.mpcorb.H

  .. autoattribute:: pykep.planet.mpcorb.n_observations

  .. autoattribute:: pykep.planet.mpcorb.n_oppositions

  .. autoattribute:: pykep.planet.mpcorb.year_of_discovery

------------

.. autoclass:: pykep.planet.tle(*args)

  .. automethod:: pykep.planet.tle.__init__(*args)

------------

.. autoclass:: pykep.planet.spice(*args)

  .. automethod:: pykep.planet.spice.__init__(*args)

------------

.. autoclass:: pykep.planet.gtoc2(*args)

  .. automethod:: pykep.planet.gtoc2.__init__(*args)

------------

.. autoclass:: pykep.planet.gtoc5(*args)

  .. automethod:: pykep.planet.gtoc5.__init__(*args)

------------

.. autoclass:: pykep.planet.gtoc6(*args)

  .. automethod:: pykep.planet.gtoc6.__init__(*args)

------------

.. autoclass:: pykep.planet.gtoc7(*args)

  .. automethod:: pykep.planet.gtoc7.__init__(*args)

