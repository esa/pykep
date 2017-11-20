.. _phasing:

==================
The phasing module
==================

==================================================       =========       ================================================
Name                                                     Type            Description
==================================================       =========       ================================================
:class:`pykep.phasing.knn`                               class           Finds the nearest-neighbour in a large list of :py:class:`pykep.planet`
:class:`pykep.phasing.dbscan`                            class           Detects clusters in a large list of :py:class:`pykep.planet`
:func:`pykep.phasing.three_impulses_approx`              function        Computes the orbital transfer cost between two :py:class:`pykep.planet`
==================================================       =========       ================================================

Detailed Documentation
======================

.. autoclass:: pykep.phasing.knn(*args)
 
   .. automethod:: pykep.phasing.knn.__init__(*args)

   .. automethod:: pykep.phasing.knn.find_neighbours(*args)

------------
       
.. autoclass:: pykep.phasing.dbscan(*args)
 
   .. automethod:: pykep.phasing.dbscan.__init__(*args)

   .. automethod:: pykep.phasing.dbscan.cluster(*args)

   .. automethod:: pykep.phasing.dbscan.pretty(*args)

   .. automethod:: pykep.phasing.dbscan.plot(*args)

------------

.. autofunction:: pykep.phasing.three_impulses_approx(*args)
