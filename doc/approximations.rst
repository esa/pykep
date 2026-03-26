.. _approximations:

Various approximations for orbital transfers
############################################

Computing the exact optimal transfer between orbits is often more expensive than allowed in
preliminary phases of the mission design. For this reason, **pykep** offers a number of approximations
that can be used as surrogates of the actual complete computation. These can be used to quickly pre-screen 
a large number of possible transfers, to then focus on the most promising ones.
Many of these approximations have been used and developed during the various edition of the
[GTOC competition](https://sophia.estec.esa.int/gtoc_portal/) and can
be found in our papers, e.g. :cite:p:`approximations`, :cite:p:`gtoc12`.

.. currentmodule:: pykep

Basic transfers
----------------
Sometimes basic ideal transfers models can be used to get some bounds or information on the possible
difficulty of a certain orbital geometry.

.. autofunction:: hohmann

.. autofunction:: bielliptic

Maximum initial mass approximation (MIMA)
-----------------------------------------
If we are computing the low-thrust transfer between two arbitrary orbits and we know the starting and final time of arrival,
these approximations can be used to compute the maximum initial mass (MIM) that will result in a feasible transfer.
In other words if our spacecraft is fat (heavier than the MIM), the transfer will not be feasible.
In general the best way to approximate the MIM is via neural method and machine learning :cite:p:`acciarini2024computing`,
but that requires learning from a vast database of precomputed trajectories and sometimes we do not have that luxury.

.. autofunction:: mim_from_hop

.. autofunction:: mima

.. autofunction:: mima_from_hop

.. autofunction:: mima2

.. autofunction:: mima2_from_hop


Phasing indicators
----------------------

.. currentmodule:: pykep.utils

.. autoclass:: knn
    :members: find_neighbours