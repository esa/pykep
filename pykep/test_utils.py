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
    def setUp(self):
        self.planets = [
            _pk.planet(_pk.udpla.jpl_lp(name)) for name in ("earth", "mars", "venus")
        ]
        self.when = _pk.epoch(0.0)

    def test_single_neighbour(self):
        for metric in ("orbital", "euclidean"):
            finder = _pk.utils.knn(self.planets, self.when, metric=metric)
            for options in ({}, {"k": 1}):
                with self.subTest(metric=metric, options=options):
                    neighbours, ids, distances = finder.find_neighbours(0, **options)
                    self.assertEqual(tuple(ids), (0,))
                    self.assertIs(neighbours[0], self.planets[0])
                    self.assertEqual(tuple(distances), (0.0,))

    def test_more_neighbours_than_planets(self):
        for metric in ("orbital", "euclidean"):
            with self.subTest(metric=metric):
                finder = _pk.utils.knn(self.planets, self.when, metric=metric)
                expected = finder.find_neighbours(0, k=len(self.planets))
                neighbours, ids, distances = finder.find_neighbours(0, k=5)
                self.assertEqual(tuple(ids), tuple(expected[1]))
                np.testing.assert_array_equal(distances, expected[2])
                self.assertEqual(len(neighbours), len(self.planets))

    def test_distance_bound_returns_only_matches(self):
        for metric in ("orbital", "euclidean"):
            with self.subTest(metric=metric):
                finder = _pk.utils.knn(self.planets, self.when, metric=metric)
                neighbours, ids, distances = finder.find_neighbours(
                    0, k=3, distance_upper_bound=1e-10
                )
                self.assertEqual(tuple(ids), (0,))
                self.assertIs(neighbours[0], self.planets[0])
                self.assertEqual(tuple(distances), (0.0,))

    def test_no_neighbours_within_bound(self):
        query = _pk.planet(_pk.udpla.jpl_lp("jupiter"))
        for metric in ("orbital", "euclidean"):
            finder = _pk.utils.knn(self.planets, self.when, metric=metric)
            for k in (1, 3):
                with self.subTest(metric=metric, k=k):
                    result = finder.find_neighbours(query, k=k, distance_upper_bound=1e-10)
                    self.assertEqual(result, ([], [], []))

    def test_requested_ranks_skip_missing_neighbours(self):
        for metric in ("orbital", "euclidean"):
            with self.subTest(metric=metric):
                finder = _pk.utils.knn(self.planets, self.when, metric=metric)
                _, all_ids, all_distances = finder.find_neighbours(0, k=3)
                neighbours, ids, distances = finder.find_neighbours(0, k=[1, 3, 5])
                self.assertEqual(tuple(ids), (all_ids[0], all_ids[2]))
                np.testing.assert_array_equal(distances, [all_distances[0], all_distances[2]])
                self.assertEqual(len(neighbours), 2)

    def test_ball_query_is_unchanged(self):
        for metric in ("orbital", "euclidean"):
            with self.subTest(metric=metric):
                finder = _pk.utils.knn(self.planets, self.when, metric=metric)
                neighbours, ids, distances = finder.find_neighbours(0, query_type="ball", r=1e-10)
                self.assertEqual(tuple(ids), (0,))
                self.assertIs(neighbours[0], self.planets[0])
                self.assertEqual(tuple(distances), (None,))
                query = _pk.planet(_pk.udpla.jpl_lp("jupiter"))
                self.assertEqual(
                    finder.find_neighbours(query, query_type="ball", r=1e-10),
                    ([], [], []),
                )
