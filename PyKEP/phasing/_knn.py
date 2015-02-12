from PyKEP.core import AU, EARTH_VELOCITY


class knn():

    """
    The class helps finding the k-nearest neighbours to a given planet from a list of planets.
    The problem of finding who is "close-by" can be efficiently solved using appropriate data structures.
    Here a kdtree is employed bringng complexity down to O(log N). The k-d-tree can also be queried efficiently
    for all asteroid within a given distance ('ball' query)
    """

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

        if type(eph) == tuple:
            eph = np.hstack(eph)
        e = eph.reshape(-1, 6)

        # normalize r
        e[:, :3] /= self._ref_r
        # normalize v
        e[:, -3:] /= self._ref_v

        return e.reshape(eph.shape)

    def _make_kdtree(self, t):
        """
        Returns a kd-tree data structure indexing the normalized ephemerides
        of planets_list, at epoch t.

        make_kdtree( planets_list, t, ref_r=AU, ref_v=EARTH_VELOCITY )

        - planets_list: list of PyKEP.planet objects
        - t: epoch

        The returned kd-tree can then be used for efficient nearest-neighbor queries.

        See also:
            http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
            https://en.wikipedia.org/wiki/K-d_tree

        Examples::

          ast_list = [PyKEP.planet_gtoc7(idx) for idx in range(1,16256)]
          kdt = PyKEP.phasing.make_kdtree(ast_list,ePyKEP.epoch(6754.3))
        """

        from scipy.spatial import cKDTree
        import numpy as np

        e = np.array([a.eph(t) for a in self._asteroids])

        # reshape memory area, so each asteroid's ephemeride gets represented by
        # a single 6 dimensional vector
        e = e.reshape((e.shape[0], -1))

        # normalize the full matrix
        # (if the `normalize_f` is set to None, no normalization takes place)
        self._eph_normalize(e)

        return cKDTree(e)

    def __init__(self, planet_list, t, ref_r=AU, ref_v=EARTH_VELOCITY):
        """
        USAGE: knn = knn(planet_list, t, ref_r=AU, ref_v=EARTH_VELOCITY):

        - planet_list   list of PyKEP planets (typically thousands)
        - t             epoch
        - ref_r         reference radius
        - ref_v         reference velocity
        """
        import numpy as np
        import PyKEP as pk
        self._asteroids = np.array(planet_list, dtype=np.object)
        self._ref_r = ref_r
        self._ref_v = ref_v
        self._t = t
        self._kdtree = self._make_kdtree(self._t)

    def find_neighbours(self, query_planet, query_type='knn', *args, **kwargs):
        """
        Finds the neighbours of a given planet at a given epoch. The user may query for the
        k-nearest neighbours or for all neighbours within a given distance

        knn.find_neighbours(query_planet, query_type='knn', \*args, \*\*kwargs )

        - query_planet: the planet we want to find neighbours of. Can be an integer, in which case it refers to the idx in self.asteroid_list
        - query_type: one of 'knn' or 'ball'.
        - \*args, \*\*args: according to the query type (read below)

        Returns (neighb, neighb_ids, dists), where dist is only computed if 'knn' type query is made

        The following kinds of spatial queries are currently implemented:

        query_type = 'knn':
            Obtain the `k` nearest asteroids
            For arguments, see:
            http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html

        query_type = 'ball':
            Obtain all asteroids within a radius `r` of the given query asteroid
            For arguments, see:
            http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
        """
        if type(query_planet) == int:
            query_planet = self._asteroids[query_planet]

        # generate the query vector
        x = query_planet.eph(self._t)
        x = self._eph_normalize(x)

        if query_type == 'knn':
            # Query for the k nearest neighbors
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html
            dists, idxs = self._kdtree.query(x, *args, **kwargs)
        elif query_type == 'ball':
            # Query for all neighbors within a sphere of given radius
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_point.html
            idxs = self._kdtree.query_ball_point(x, *args, **kwargs)
            dists = [None] * len(idxs)
        else:
            raise Exception('Unrecognized query type: %s' % str(query_type))

        neighb = [
            # (ast. object, ast. ID, distance)
            (self._asteroids[i], i, d)
            for i, d in zip(idxs, dists)
        ]

        # split into three lists, one of objects, one of IDs, and one for
        # distances
        neighb, neighb_ids, dists = list(
            zip(*neighb)) if neighb != [] else ([], [], [])

        return neighb, neighb_ids, dists
