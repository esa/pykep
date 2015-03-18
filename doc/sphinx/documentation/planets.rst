.. _planets:

==================
The planets module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`PyKEP.planets.keplerian`                         class           A simple planet with Keplerian ephemerides
:class:`PyKEP.planets.jpl_low_precision`                 class           A solar system planet using jpl low-precision ephemerides
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

