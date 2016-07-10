.. _phasing:

==================
The phasing module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`PyKEP.phasing.knn`                               class           Finds the nearest-neighbour in a large list of :py:class:`PyKEP.planet`
:class:`PyKEP.phasing.dbscan`                            class           Detects clusters in a large list of :py:class:`PyKEP.planet`
:func:`PyKEP.phasing.three_impulses_approx`              function        Computes the orbital transfer cost between two :py:class:`PyKEP.planet`
==================================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: PyKEP.phasing.knn(*args)
 
   .. automethod:: PyKEP.phasing.knn.__init__(*args)

   .. automethod:: PyKEP.phasing.knn.find_neighbours(*args)

------------
       
.. autoclass:: PyKEP.phasing.dbscan(*args)
 
   .. automethod:: PyKEP.phasing.dbscan.__init__(*args)

   .. automethod:: PyKEP.phasing.dbscan.cluster(*args)

   .. automethod:: PyKEP.phasing.dbscan.pretty(*args)

   .. automethod:: PyKEP.phasing.dbscan.plot(*args)

------------

.. autofunction:: PyKEP.phasing.three_impulses_approx(*args)
