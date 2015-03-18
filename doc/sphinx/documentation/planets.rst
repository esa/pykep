.. _planets:

==================
The planets module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`PyKEP.planets.keplerian`                         class           A simple planet with Keplerian ephemerides
:class:`PyKEP.planets.jpl_lp`                            class           A solar system planet using jpl low-precision ephemerides
:class:`PyKEP.planets.tle`                               class           An Earth artificial satellite from its TLE (ephemerides are computed via the SGP4 propagator)
:class:`PyKEP.planets.spice`                             class           A planet with ephemerides computed using the JPL SPICE Toolbox (requires BUILD_SPICE option active when building from cmake)
:class:`PyKEP.planets.mpcorb`                            class           A planet from the MPCORB database (keplerian ephemerides)
:class:`PyKEP.planets.gtoc2`                             class           An asteroid from the GTOC2 competition (keplerian ephemerides)
:class:`PyKEP.planets.gtoc5`                             class           An asteroid from the GTOC5 competition (keplerian ephemerides)
:class:`PyKEP.planets.gtoc6`                             class           A Jupiter moon from the GTOC6 competition (keplerian ephemerides)
:class:`PyKEP.planets.gtoc7`                             class           An asteroid from the GTOC7 competition (keplerian ephemerides)
==================================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.planets.keplerian(*args)

  .. automethod:: PyKEP.planets.keplerian.__init__(*args)

------------

.. autoclass:: PyKEP.planets.jpl_lp(*args)

  .. automethod:: PyKEP.planets.jpl_lp.__init__(*args)

------------

.. autoclass:: PyKEP.planets.mpcorb(*args)

  .. automethod:: PyKEP.planets.mpcorb.__init__(*args)

  .. autoattribute:: PyKEP.planets.mpcorb.H

  .. autoattribute:: PyKEP.planets.mpcorb.n_observations

  .. autoattribute:: PyKEP.planets.mpcorb.n_oppositions

  .. autoattribute:: PyKEP.planets.mpcorb.year_of_discovery

------------

.. autoclass:: PyKEP.planets.tle(*args)

  .. automethod:: PyKEP.planets.tle.__init__(*args)

------------

.. autoclass:: PyKEP.planets.spice(*args)

  .. automethod:: PyKEP.planets.spice.__init__(*args)

------------

.. autoclass:: PyKEP.planets.gtoc2(*args)

  .. automethod:: PyKEP.planets.gtoc2.__init__(*args)

------------

.. autoclass:: PyKEP.planets.gtoc5(*args)

  .. automethod:: PyKEP.planets.gtoc5.__init__(*args)

------------

.. autoclass:: PyKEP.planets.gtoc6(*args)

  .. automethod:: PyKEP.planets.gtoc6.__init__(*args)

------------

.. autoclass:: PyKEP.planets.gtoc7(*args)

  .. automethod:: PyKEP.planets.gtoc7.__init__(*args)

