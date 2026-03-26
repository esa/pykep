## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
## This file is part of the kep3 library.
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pykep as pk

class knn:
    R"""
    The class finds the k-nearest neighbours to a given planet
    at some epoch from a list of planets.
    The idea being that under some definition of "closeness"
    close-by planets are good candidates for orbital transfers (i.e. resulting in a low :math:`\Delta V`).
    The use of this class is thus in preliminary mission phases of
    multiple randevous trajectories where target selection is to be performed efficiently.

    The problem of finding who is "close-by" under a given metric can be efficiently solved
    using appropriate data structures. Here a kdtree is employed bringing complexity down to O(log N).
    The k-d-tree can then be queried efficiently for all asteroid within a given distance ('ball' query)
    or for all k closest asteroids ('knn' query).

    The notion of distance used (metric) can be:

    - 'euclidean': the simple Euclidean distance over the asteroid's (r,v).
      The position and velocity vectors are scaled w.r.t. some reference values.

    - 'orbital', the distance is computed with respect to :math:`\frac{r}{T} + v`, :math:`\frac{r}{T}`,
      corresponding to the :math:`\Delta V` computed over a linear model of the orbital transfer.
      The distance returned will thus be in m/s.

    The class is initialized with a list of planets, an epoch, and the metric to be used.
    """

    from pykep.core import AU, EARTH_VELOCITY

    def _eph_normalize(self, eph):
        """
        Normalize a body's ephemerides:
            ( r / refn_r, v / refn_v )
        where `refn_r` and `refn_v` are reference values you can set in the
        arguments (defaults: refn_r == pk.AU, refn_v == pk.EARTH_VELOCITY).

        Accepts either a single ephemeride, or a matrix with one per row.

        Example:
        >>> eph_normalize( asteroid_list[0].eph(TRAJ_START_MIN) )
        array([ -1.80619704e-01,   9.66554711e-01,  -1.46653866e-05,
                -9.99416814e-01,  -1.87484002e-01,   4.12493855e-06])
        """
        import numpy as np

        if type(eph) == list:
            eph = np.hstack(eph)
        e = eph.reshape(-1, 6)

        # normalize r
        e[:, :3] /= self._ref_r
        # normalize v
        e[:, -3:] /= self._ref_v

        return e.reshape(eph.shape)

    def _orbital_metric(self, r, v):
        """
        Compute the orbital metric for given position and velocity vectors.

        Args:
            *r* (:class:`list`): position vector.

            *v* (:class:`list`): velocity vector.

        Returns:
            :class:`tuple`[:class:`list`, :class:`list`]: (DV1, DV2) where DV1 and DV2 are components of the orbital metric.
        """
        DV2 = [a / (self._tof * pk.DAY2SEC) for a in r]
        DV1 = [a + b for a, b in zip(DV2, v)]
        return (DV1, DV2)

    def _make_kdtree(self, t):
        """
        Returns a kd-tree data structure indexing the normalized ephemerides
        of planets_list, at epoch t.

        Args:
            *t* (:class:`~pykep.epoch`): epoch.

        Returns:
            :class:`scipy.spatial.cKDTree`: kd-tree data structure.

        See also:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
            https://en.wikipedia.org/wiki/K-d_tree

        Examples::

          ast_list = [pykep.planet_gtoc7(idx) for idx in range(1,16256)]
          kdt = pykep.phasing.make_kdtree(ast_list,epykep.epoch(6754.3))
        """

        from scipy.spatial import cKDTree
        import numpy as np

        if self._metric == "euclidean":
            e = np.array([a.eph(t) for a in self._asteroids])
        else:
            e = np.array([self._orbital_metric(*a.eph(t)) for a in self._asteroids])

        # reshape memory area, so each asteroid's ephemeride gets represented by
        # a single 6 dimensional vector
        e = e.reshape((e.shape[0], -1))

        # normalize the full matrix
        # (if the `normalize_f` is set to None, no normalization takes place)
        if self._metric == "euclidean":
            self._eph_normalize(e)

        return cKDTree(e)

    def __init__(
        self,
        planet_list,
        when,
        metric="orbital",
        ref_r=pk.AU,
        ref_v=pk.EARTH_VELOCITY,
        tof=180.0,
    ):
        """
        Initializes the knn class.

        Args:
            *planet_list* (:class:`list` of :class:`~pykep.planet`): list of pykep planets (typically thousands).

            *when* (:class:`~pykep.epoch`): epoch.

            *metric* (:class:`str`, optional): one of ['euclidean', 'orbital']. Defaults to 'orbital'.

            *ref_r* (:class:`float`, optional): reference radius (used as a scaling factor for r if the metric is 'euclidean'). Defaults to AU.

            *ref_v* (:class:`float`, optional): reference velocity (used as a scaling factor for v if the metric is 'euclidean'). Defaults to EARTH_VELOCITY.

            *tof* (:class:`float`, optional): average transfer time in days (used in the definition of the 'orbital' metric). Defaults to 180.0.

        Example::

            from pykep import *
            pl_list = ...... # a list of planets
            knn = phasing.knn(pl_list, epoch(t0), metric='orbital', tof=180)
            neighb, ids, dists = knn.find_neighbours(pl_list[ast_0], query_type='knn', k=10000)
            neighb, ids, _ = knn.find_neighbours(pl_list[ast_0], query_type='ball', r=5000)
        """
        import numpy as np

        self._asteroids = np.array(planet_list)
        self._ref_r = ref_r
        self._ref_v = ref_v
        self._when = when
        self._metric = metric
        self._tof = tof
        self._kdtree = self._make_kdtree(self._when)

    def find_neighbours(self, query_planet, query_type="knn", *args, **kwargs):
        """
        Finds the neighbours of a given planet at a given epoch. The user may query for the
        k-nearest neighbours or for all neighbours within a given distance.

        Args:
            *query_planet* (:class:`~pykep.planet` or :class:`int`): the planet we want to find neighbours of. Can be an integer, in which case it refers to the idx in self.asteroid_list.

            *query_type* (:class:`str`, optional): one of 'knn' or 'ball'. Defaults to 'knn'.

            *\\*args*: according to the query type (read below).

            *\\*\\*kwargs*: according to the query type (read below).

        Returns:
            :class:`tuple`: (neighb, neighb_ids, dists), where dist is only computed if 'knn' type query is made.

        The following kinds of spatial queries are currently implemented:

        query_type = 'knn':
            The kwarg 'k' determines how many k-nearest neighbours are returned.
            For arguments, see:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html

        query_type = 'ball':
            The kwarg 'r' determines the distance within which all asteroids are returned.
            For arguments, see:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
        """
        if type(query_planet) == int:
            query_planet = self._asteroids[query_planet]

        # generate the query vector
        x = query_planet.eph(self._when)
        if self._metric == "euclidean":
            x = self._eph_normalize(x)
        else:
            DV1, DV2 = self._orbital_metric(x[0], x[1])
            x = DV1 + DV2

        if query_type == "knn":
            # Query for the k nearest neighbors
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html
            dists, idxs = self._kdtree.query(x, *args, **kwargs)
        elif query_type == "ball":
            # Query for all neighbors within a sphere of given radius
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
            idxs = self._kdtree.query_ball_point(x, *args, **kwargs)
            dists = [None] * len(idxs)
        else:
            raise Exception("Unrecognized query type: %s" % str(query_type))

        neighb = [
            # (ast. object, ast. ID, distance)
            (self._asteroids[i], i, d)
            for i, d in zip(idxs, dists)
        ]

        # split into three lists, one of objects, one of IDs, and one for
        # distances
        neighb, neighb_ids, dists = list(zip(*neighb)) if neighb != [] else ([], [], [])

        return neighb, neighb_ids, dists
