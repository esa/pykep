# Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
#                          Advanced Concepts Team, European Space Agency (ESA)
#
# This file is part of the pykep library.
#
# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import pykep as _pk
import pygmo as pg
import numpy as np

import unittest as _ut


def float_rel_error(a: float, b: float):
    return abs(a - b) / abs(a)


class encoding_tests(_ut.TestCase):
    def test_alpha_direct_conversion(self):
        import pykep as _pk

        tofs = [12.34, 232.2, 23.45, 134.3]
        alphas, T = _pk.direct2alpha(tofs)
        tofs_from_alphas = _pk.alpha2direct(alphas, T)
        err = [a - b for a, b in zip(tofs, tofs_from_alphas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

        tofs = np.random.random((4,)) * 20
        alphas, T = _pk.direct2alpha(tofs)
        tofs_from_alphas = _pk.alpha2direct(alphas, T)
        err = [a - b for a, b in zip(tofs, tofs_from_alphas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

    def test_eta_direct_conversion(self):
        import pykep as _pk

        tofs = [12.34, 232.2, 23.45, 134.3]
        tmax = 300
        etas = _pk.direct2eta(tofs, tmax)
        tofs_from_etas = _pk.eta2direct(etas, tmax)
        err = [a - b for a, b in zip(tofs, tofs_from_etas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

        tofs = np.random.random((4,)) * 100
        tmax = 400
        etas = _pk.direct2eta(tofs, tmax)
        tofs_from_etas = _pk.eta2direct(etas, tmax)
        err = [a - b for a, b in zip(tofs, tofs_from_etas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)
        
    def test_uvV_cartesian_conversion(self):
        import pykep as _pk
        
        vector = np.random.random((3,)) 
        uvV = _pk.cartesian2uvV(vector)
        vector_new = _pk.uvV2cartesian(uvV)
        err = [a - b for a, b in zip(vector, vector_new)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

class knn_tests(_ut.TestCase):
    def _planets(self):
        import pykep as _pk

        return [
            _pk.planet(_pk.udpla.jpl_lp(name)) for name in ("earth", "mars", "venus")
        ]

    def test_single_neighbour_query(self):
        import pykep as _pk

        for metric in ("orbital", "euclidean"):
            knn = _pk.utils.knn(self._planets(), _pk.epoch(0.0), metric=metric)
            neighb, ids, dists = knn.find_neighbours(0)
            self.assertEqual(len(neighb), 1)
            self.assertEqual(len(ids), 1)
            self.assertEqual(len(dists), 1)
            self.assertEqual(ids[0], 0)
            self.assertEqual(dists[0], 0.0)

    def test_more_neighbours_requested_than_available(self):
        import pykep as _pk

        planets = self._planets()
        for metric in ("orbital", "euclidean"):
            knn = _pk.utils.knn(planets, _pk.epoch(0.0), metric=metric)
            neighb, ids, dists = knn.find_neighbours(0, k=5)
            self.assertEqual(len(neighb), len(planets))
            self.assertEqual(len(ids), len(planets))
            self.assertEqual(len(dists), len(planets))

    def test_distance_upper_bound_keeps_only_reachable_neighbours(self):
        import pykep as _pk

        for metric in ("orbital", "euclidean"):
            knn = _pk.utils.knn(self._planets(), _pk.epoch(0.0), metric=metric)
            neighb, ids, dists = knn.find_neighbours(
                0, k=3, distance_upper_bound=1e-10
            )
            self.assertEqual(len(neighb), 1)
            self.assertEqual(ids[0], 0)

    def test_no_neighbour_within_distance_upper_bound(self):
        import pykep as _pk

        outsider = _pk.planet(_pk.udpla.jpl_lp("jupiter"))
        for metric in ("orbital", "euclidean"):
            knn = _pk.utils.knn(self._planets(), _pk.epoch(0.0), metric=metric)
            self.assertEqual(
                knn.find_neighbours(outsider, k=3, distance_upper_bound=1e-10),
                ([], [], []),
            )

    def test_nominal_query_unchanged(self):
        import pykep as _pk

        for metric in ("orbital", "euclidean"):
            knn = _pk.utils.knn(self._planets(), _pk.epoch(0.0), metric=metric)
            neighb, ids, dists = knn.find_neighbours(0, k=2)
            self.assertEqual(len(neighb), 2)
            self.assertEqual(ids[0], 0)

            neighb, ids, dists = knn.find_neighbours(0, "ball", r=1e12)
            self.assertEqual(len(neighb), 3)
