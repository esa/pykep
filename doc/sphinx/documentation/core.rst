.. _core:

===============
The core module
===============

The core module, imported in the main PyKEP namespace, contains the classes and helper functions to perform basic space-flight mechanics computations. You can browse the list of classes and functions below. Find the detailed documentation later on this page.

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`PyKEP.epoch`                            class           represents an epoch (i.e. a fixed point in time) 
:func:`PyKEP.epoch_from_string`                 function        helper function to construct an epoch from a string containing a date in the format YYYY-MM-DD HH:MM:SS
:func:`PyKEP.epoch_from_iso_string`             function        helper function to construct an epoch from a string containing a date in the ISO format YYYYMMDDTHHMMSS
:class:`PyKEP.planet`                           class           a generic planet, ephemerides are keplerian
:class:`PyKEP.planet_ss`                        class           a solar system planet, ephemerides are the `JPL low-precision <http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf>`_
:class:`PyKEP.planet_mpcorb`                    class           represents an asteroid from the `MPCORB database <http://www.minorplanetcenter.org/iau/MPCORB.html>`_ 
:class:`PyKEP.planet_tle`                       class           represents a satellite defined by its `TLEs <http://www.celestrak.com/NORAD/elements/>`_, ephemerides will be computed using an SGP4 propagator
:class:`PyKEP.lambert_problem`                  class           solves the multirevolution lambert problem
:func:`PyKEP.propagate_lagrangian`              function        propagates pure keplerian motion using Lagrange coefficients and universal variables
:func:`PyKEP.propagate_taylor`                  function        propagates keplerian motion disturbed by a constant inertial thrust using Taylor integration method
:func:`PyKEP.fb_con`                            function        returns violation of velocity and angular constraint during a fly-by
:func:`PyKEP.fb_vel`                            function        returns the violation of the velocity and angular constraint during a fly-by in terms of one single DV
:func:`PyKEP.fb_prop`                           function        propoagates forward a fly-by hyperbola returning the new inetrial velocity of a spacecraft after the planetary encounter
:func:`PyKEP.barker`                            function        computes the (parabolic) time-of-flight from the Barker equation
:func:`PyKEP.ic2par`                            function        Transforms r and v into the osculating orbital elements
:func:`PyKEP.read_satcat`                       function        Reads the `SATCAT database <http://celestrak.com/satcat/search.asp>`_  into a dictionary 
:func:`PyKEP.read_tle`                          function        Reads a TLE file. These can be downloaded from `Celestrack <http://www.celestrak.com/NORAD/elements/>`_ or, if you have an account, `Space Track <https://www.space-track.org/>`_

=========================================       =========       ================================================


Detailed Documentation
======================

.. autoclass:: PyKEP.lambert_problem(*args)

  .. automethod:: PyKEP.lambert_problem.__init__(*args)

  .. automethod:: PyKEP.lambert_problem.get_v1()

  .. automethod:: PyKEP.lambert_problem.get_v2()

  .. automethod:: PyKEP.lambert_problem.get_x()

  .. automethod:: PyKEP.lambert_problem.get_Nmax()

  .. automethod:: PyKEP.lambert_problem.get_iters()

------------

.. autoclass:: PyKEP.epoch(*args)

  .. automethod:: PyKEP.epoch.__init__(*args)

  .. autoattribute:: PyKEP.epoch.jd

  .. autoattribute:: PyKEP.epoch.mjd

  .. autoattribute:: PyKEP.epoch.mjd2000

------------

.. autofunction:: PyKEP.epoch_from_string(s)

------------

.. autofunction:: PyKEP.epoch_from_iso_string(s)

------------

.. autoclass:: PyKEP.planet(*args)

  .. automethod:: PyKEP.planet.__init__(*args)

  .. automethod:: PyKEP.planet.eph(*args)

  .. autoattribute:: PyKEP.planet.orbital_elements

  .. autoattribute:: PyKEP.planet.radius

  .. autoattribute:: PyKEP.planet.ref_epoch

------------

.. autoclass:: PyKEP.planet_ss(*args)

  .. automethod:: PyKEP.planet_ss.__init__(*args)

------------

.. autoclass:: PyKEP.planet_mpcorb(*args)

  .. automethod:: PyKEP.planet_mpcorb.__init__(*args)

  .. autoattribute:: PyKEP.planet_mpcorb.H

  .. autoattribute:: PyKEP.planet_mpcorb.n_observations

  .. autoattribute:: PyKEP.planet_mpcorb.n_oppositions

  .. autoattribute:: PyKEP.planet_mpcorb.year_of_discovery

------------

.. autoclass:: PyKEP.planet_tle(*args)

  .. automethod:: PyKEP.planet_tle.__init__(*args)

------------

.. autofunction:: PyKEP.propagate_lagrangian(*args)

------------

.. autofunction:: PyKEP.propagate_taylor(*args)

------------

.. autofunction:: PyKEP.fb_con(*args)

------------

.. autofunction:: PyKEP.fb_prop(*args)

------------

.. autofunction:: PyKEP.fb_vel(*args)

------------

.. autofunction:: PyKEP.barker(*args)

------------

.. autofunction:: PyKEP.ic2par(*args)

------------

.. autofunction:: PyKEP.read_satcat(*args)

------------

.. autofunction:: PyKEP.read_tle(*args)