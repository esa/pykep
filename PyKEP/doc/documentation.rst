===================
PyKEP Documentation
===================

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

  .. automethod:: PyKEP.planet.orbital_elements(*args)

.. autoclass:: PyKEP.planet_ss(*args)

  .. automethod:: PyKEP.planet_ss.__init__(*args)

.. autoclass:: PyKEP.planet_mpcorb(*args)

  .. automethod:: PyKEP.planet_mpcorb.__init__(*args)

.. autoclass:: PyKEP.planet_gtoc5(*args)

  .. automethod:: PyKEP.planet_gtoc5.__init__(*args)

.. autoclass:: PyKEP.planet_gtoc2(*args)

  .. automethod:: PyKEP.planet_gtoc2.__init__(*args)

.. autofunction:: PyKEP.propagate_lagrangian(*args)


