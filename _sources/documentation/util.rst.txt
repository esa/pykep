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
:func:`pykep.util.load_gravity_model`                    function        Loads a spherical harmonics gravity model
:func:`pykep.util.gravity_spherical_harmonic`            function        Calculates the gravitational acceleration due to the provided spherical harmonic gravity model
==================================================       =========       ================================================

Detailed Documentation
======================

.. autofunction:: pykep.util.read_satcat(*args)

------------

.. autofunction:: pykep.util.read_tle(*args)

------------

.. autofunction:: pykep.util.load_spice_kernel(*args)

------------

.. autofunction:: pykep.util.load_gravity_model(*args)

------------

.. autofunction:: pykep.util.gravity_spherical_harmonic(*args)
