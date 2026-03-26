## Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
## (bluescarni@gmail.com)
##
## This file is part of the pykep library.
##
## This Source Code Form is subject to the terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pykep as pk
import pygmo as pg
import numpy as np

import unittest as _ut


def float_rel_error(a: float, b: float):
    return abs(a - b) / abs(a)


class encoding_tests(_ut.TestCase):
    def test_alpha_direct_conversion(self):
        import pykep as pk

        tofs = [12.34, 232.2, 23.45, 134.3]
        alphas, T = pk.direct2alpha(tofs)
        tofs_from_alphas = pk.alpha2direct(alphas, T)
        err = [a - b for a, b in zip(tofs, tofs_from_alphas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

        tofs = np.random.random((4,)) * 20
        alphas, T = pk.direct2alpha(tofs)
        tofs_from_alphas = pk.alpha2direct(alphas, T)
        err = [a - b for a, b in zip(tofs, tofs_from_alphas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

    def test_eta_direct_conversion(self):
        import pykep as pk

        tofs = [12.34, 232.2, 23.45, 134.3]
        tmax = 300
        etas = pk.direct2eta(tofs, tmax)
        tofs_from_etas = pk.eta2direct(etas, tmax)
        err = [a - b for a, b in zip(tofs, tofs_from_etas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)

        tofs = np.random.random((4,)) * 100
        tmax = 400
        etas = pk.direct2eta(tofs, tmax)
        tofs_from_etas = pk.eta2direct(etas, tmax)
        err = [a - b for a, b in zip(tofs, tofs_from_etas)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)
        
    def test_uvV_cartesian_conversion(self):
        import pykep as pk
        
        vector = np.random.random((3,)) 
        uvV = pk.cartesian2uvV(vector)
        vector_new = pk.uvV2cartesian(uvV)
        err = [a - b for a, b in zip(vector, vector_new)]
        err = np.linalg.norm(err)
        self.assertTrue(err < 1e-13)