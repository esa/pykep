.. _core:

===============
The core module
===============

The core module, imported in the main pykep namespace, contains the classes and helper functions to perform basic space-flight mechanics computations. You can browse the list of classes and functions below. Find the detailed documentation later on this page.

=========================================       =========       ================================================
Name                                            Type            Description
=========================================       =========       ================================================
:class:`pykep.epoch`                            class           represents an epoch (i.e. a fixed point in time)
:func:`pykep.epoch_from_string`                 function        helper function to construct an epoch from a string containing a date in the format YYYY-MM-DD HH:MM:SS
:func:`pykep.epoch_from_iso_string`             function        helper function to construct an epoch from a string containing a date in the ISO format YYYYMMDDTHHMMSS
:class:`pykep.lambert_problem`                  class           solves the multirevolution lambert problem
:func:`pykep.propagate_lagrangian`              function        propagates pure keplerian motion using Lagrange coefficients and universal variables
:func:`pykep.propagate_taylor`                  function        propagates keplerian motion disturbed by a constant inertial thrust using Taylor integration method
:func:`pykep.fb_con`                            function        returns violation of velocity and angular constraint during a fly-by
:func:`pykep.fb_vel`                            function        returns the violation of the velocity and angular constraint during a fly-by in terms of one single DV
:func:`pykep.fb_prop`                           function        propagates forward a fly-by hyperbola returning the new inetrial velocity of a spacecraft after the planetary encounter
:func:`pykep.barker`                            function        computes the (parabolic) time-of-flight from the Barker equation
:func:`pykep.ic2par`                            function        Transforms r and v into the osculating orbital elements
:func:`pykep.par2ic`                            function        Transforms osculating orbital elements into r and v
:func:`pykep.par2eq`                            function        Transforms osculating orbital elements into modified equinoctial
:func:`pykep.eq2par`                            function        Transforms modified equinoctial into osculating orbital elements 
:func:`pykep.eq2ic`                             function        Transforms modified equinoctial elements into r and v
:func:`pykep.ic2eq`                             function        Transforms r and v into modified equinoctial elements
=========================================       =========       ================================================


Detailed Documentation
======================

.. autoclass:: pykep.lambert_problem(*args)

  .. automethod:: pykep.lambert_problem.__init__(*args)

  .. automethod:: pykep.lambert_problem.get_v1()

  .. automethod:: pykep.lambert_problem.get_v2()

  .. automethod:: pykep.lambert_problem.get_x()

  .. automethod:: pykep.lambert_problem.get_Nmax()

  .. automethod:: pykep.lambert_problem.get_iters()

------------

.. autoclass:: pykep.epoch(*args)

  .. automethod:: pykep.epoch.__init__(*args)

  .. autoattribute:: pykep.epoch.jd

  .. autoattribute:: pykep.epoch.mjd

  .. autoattribute:: pykep.epoch.mjd2000

------------

.. autofunction:: pykep.epoch_from_string(s)

------------

.. autofunction:: pykep.epoch_from_iso_string(s)

------------

.. autofunction:: pykep.propagate_lagrangian(*args)

------------

.. autofunction:: pykep.propagate_taylor(*args)

------------

.. autofunction:: pykep.fb_con(*args)

------------

.. autofunction:: pykep.fb_prop(*args)

------------

.. autofunction:: pykep.fb_vel(*args)

------------

.. autofunction:: pykep.barker(*args)

------------

.. autofunction:: pykep.ic2par(*args)

------------

.. autofunction:: pykep.par2ic(*args)

------------

.. autofunction:: pykep.par2eq(*args)

------------

.. autofunction:: pykep.eq2par(*args)

------------

.. autofunction:: pykep.eq2ic(*args)

------------

.. autofunction:: pykep.ic2eq(*args)
