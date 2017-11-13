.. _util:

==================
The utility module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:func:`pykep.util.read_satcat`                           function        Reads the `SATCAT database <http://celestrak.com/satcat/search.asp>`_  into a dictionary
:func:`pykep.util.read_tle`                              function        Reads a TLE file. These can be downloaded from `Celestrack <http://www.celestrak.com/NORAD/elements/>`_ or, if you have an account, `Space Track <https://www.space-track.org/>`_
:func:`pykep.util.load_spice_kernel`                     function        Loads in memory a kernel from the JPL SPICE toolbox (requires BUILD_SPICE option active when building from cmake)
==================================================       =========       ================================================

Detailed Documentation
======================

.. autofunction:: pykep.util.read_satcat(*args)

------------

.. autofunction:: pykep.util.read_tle(*args)

------------

.. autofunction:: pykep.util.load_spice_kernel(*args)
