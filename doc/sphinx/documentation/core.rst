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
:class:`PyKEP.lambert_problem`                  class           solves the multirevolution lambert problem
:func:`PyKEP.propagate_lagrangian`              function        propagates pure keplerian motion using Lagrange coefficients and universal variables
:func:`PyKEP.propagate_taylor`                  function        propagates keplerian motion disturbed by a constant inertial thrust using Taylor integration method
:func:`PyKEP.fb_con`                            function        returns violation of velocity and angular constraint during a fly-by
:func:`PyKEP.fb_vel`                            function        returns the violation of the velocity and angular constraint during a fly-by in terms of one single DV
:func:`PyKEP.fb_prop`                           function        propoagates forward a fly-by hyperbola returning the new inetrial velocity of a spacecraft after the planetary encounter
:func:`PyKEP.barker`                            function        computes the (parabolic) time-of-flight from the Barker equation
:func:`PyKEP.ic2par`                            function        Transforms r and v into the osculating orbital elements
:func:`PyKEP.par2ic`                            function        Transforms osculating orbital elements into r and v
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

.. autofunction:: PyKEP.par2ic(*args)
