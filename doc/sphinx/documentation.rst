===================
PyKEP Documentation
===================

List of main PyKEP elements in the core module
==============================================

=========================================	========= 	================================================
Name						Type      	Description
=========================================	=========	================================================
:class:`PyKEP.epoch`       			class      	represents an epoch (i.e. a fixed point in time)
:func:`PyKEP.epoch_from_string`			function	helper function to construct an epoch from a string containing a date in the format YYYY-MM-DD HH:MM:SS
:func:`PyKEP.epoch_from_iso_string`		function	helper function to construct an epoch from a string containing a date in the ISO format YYYYMMDDTHHMMSS
:class:`PyKEP.planet`      			class      	represents a generic planet in a keplerian motion around its primary 
:class:`PyKEP.planet_ss`      			class      	represents a solar system planet in keplerian motion around the sun
:class:`PyKEP.planet_mpcorb`      		class      	represents an asteroid from the `MPCORB database <http://www.minorplanetcenter.org/iau/MPCORB.html>`_ 
:class:`PyKEP.planet_gtoc5`      		class      	represents an asteroid from the `GTOC5 competition <http://gtoc5.math.msu.su/>`_ asteroid list
:class:`PyKEP.planet_gtoc2` 			class		represents an asteroid from the `GTOC2 competition <http://www.esa.int/gsp/ACT/mad/op/GTOC/index.htm>`_ asteroid list
:class:`PyKEP.lambert_problem`			class		solves the multirevolution lambert problem
:func:`PyKEP.propagate_lagrangian`		function	propagates pure keplerian motion using Lagrange coefficients and universal variables
:func:`PyKEP.propagate_taylor`			function	propagates keplerian motion disturbed by a constant inertial thrust using Taylor integration method
:func:`PyKEP.fb_con`     			function	computes the feasibility of a fly-by defined by a planet and the spacecarft relative in and out velocities
:func:`PyKEP.fb_prop`   			function	propoagates forward a fly-by hyperbola returning the new inetrial velocity of a spacecraft after the planetary encounter
=========================================  	=========	================================================

List of main PyKEP elements in the sims_flanagan module
=======================================================

=========================================	========= 	================================================
Name						Type      	Description
=========================================	=========	================================================
:class:`PyKEP.sims_flanagan.spacecraft`		class		represents a nuclear electric propelled spacecraft
:class:`PyKEP.sims_flanagan.sc_state`		class		represent the spacecraft state (r,v,m)
:class:`PyKEP.sims_flanagan.leg`		class		represents one leg in the Sims-Flanagan model
=========================================  	=========	================================================

List of constants defined within PyKEP
=========================================

=================================	=============================================================================
Name					Description
=================================	=============================================================================
PyKEP.AU				An astronomical unit in meters
PyKEP.DAY2SEC				Conversion factor from days to seconds
PyKEP.DAY2YEAR				Conversion factor from days to years
PyKEP.DEG2RAD				Conversion factor from degrees to radians
PyKEP.RAD2DEG				Conversion factor from radians to degree
PyKEP.SEC2DAY				Conversion factor from seconds to days
PyKEP.MU_SUN				Sun gravitational constyant in m^3/s^2
PyKEP.EARTH_VELOCITY			Square root of MU_SUN/AU. Average earth velocity in meters per seconds.
PyKEP.G0				The standard gravity acceleration at ground level 
=================================	=============================================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.lambert_problem(*args)

  .. automethod:: PyKEP.lambert_problem.__init__(*args)

  .. automethod:: PyKEP.lambert_problem.get_v1()

  .. automethod:: PyKEP.lambert_problem.get_v2()

  .. automethod:: PyKEP.lambert_problem.get_a()

  .. automethod:: PyKEP.lambert_problem.get_p()

  .. automethod:: PyKEP.lambert_problem.get_Nmax()

  .. automethod:: PyKEP.lambert_problem.get_iters()
 
  .. automethod:: PyKEP.lambert_problem.is_reliable()

.. autoclass:: PyKEP.epoch(*args)

  .. automethod:: PyKEP.epoch.__init__(*args)

  .. automethod:: PyKEP.epoch.jd() 

  .. automethod:: PyKEP.epoch.mjd()

  .. automethod:: PyKEP.epoch.mjd2000()

  .. autoattribute:: PyKEP.epoch.epoch_type

.. autofunction:: PyKEP.epoch_from_string(s)

.. autofunction:: PyKEP.epoch_from_iso_string(s)

.. autoclass:: PyKEP.planet(*args)

  .. automethod:: PyKEP.planet.__init__(*args)

  .. automethod:: PyKEP.planet.eph(*args)

  .. autoattribute:: PyKEP.planet.orbital_elements

  .. autoattribute:: PyKEP.planet.radius

  .. autoattribute:: PyKEP.planet.ref_epoch

.. autoclass:: PyKEP.planet_ss(*args)

  .. automethod:: PyKEP.planet_ss.__init__(*args)

.. autoclass:: PyKEP.planet_mpcorb(*args)

  .. automethod:: PyKEP.planet_mpcorb.__init__(*args)

  .. autoattribute:: PyKEP.planet_mpcorb.H

  .. autoattribute:: PyKEP.planet_mpcorb.n_observations

  .. autoattribute:: PyKEP.planet_mpcorb.n_oppositions

  .. autoattribute:: PyKEP.planet_mpcorb.year_of_discovery

.. autoclass:: PyKEP.planet_gtoc5(*args)

  .. automethod:: PyKEP.planet_gtoc5.__init__(*args)

.. autoclass:: PyKEP.planet_gtoc2(*args)

  .. automethod:: PyKEP.planet_gtoc2.__init__(*args)

.. autoclass:: PyKEP.sims_flanagan.spacecraft(*args)

  .. automethod:: PyKEP.sims_flanagan.spacecraft.__init__(*args)

.. autoclass:: PyKEP.sims_flanagan.sc_state(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.__init__(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.set(*args)

  .. automethod:: PyKEP.sims_flanagan.sc_state.get(*args)

.. autoclass:: PyKEP.sims_flanagan.leg(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.__init__(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.high_fidelity(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.mismatch_constraints(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.throttles_constraints(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set_mu(*args)

  .. automethod:: PyKEP.sims_flanagan.leg.set_spacecraft(*args)

.. autofunction:: PyKEP.propagate_lagrangian(*args)

.. autofunction:: PyKEP.propagate_taylor(*args)

.. autofunction:: PyKEP.fb_con(*args)

.. autofunction:: PyKEP.fb_prop(*args)
