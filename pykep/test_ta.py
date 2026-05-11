# Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
#                          Advanced Concepts Team, European Space Agency (ESA)
#
# This file is part of the pykep library.
#
# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import math
import unittest as _ut

import numpy as np
import pykep as _pk


class ta_api_tests(_ut.TestCase):
    def _assert_finite(self, seq):
        for x in seq:
            self.assertTrue(math.isfinite(float(x)))

    def test_zoh_kep_api(self):
        dyn = _pk.ta.zoh_kep_dyn()
        self.assertEqual(len(dyn), 7)

        ta = _pk.ta.get_zoh_kep(1e-12)
        self.assertEqual(len(ta.state), 7)
        self.assertEqual(len(ta.pars), 5)
        ta.time = 0.0
        ta.state[:] = [1.1, 0.2, -0.1, 0.0, 0.9, 0.2, 1.0]
        ta.pars[:] = [0.01, 1.0, 0.0, 0.0, 0.2]
        ta.propagate_until(0.2)
        self._assert_finite(ta.state)

        ta_var = _pk.ta.get_zoh_kep_var(1e-12)
        self.assertGreater(len(ta_var.state), 7)
        self.assertEqual(len(ta_var.pars), 5)
        ta_var.time = 0.0
        ta_var.state[:7] = [1.1, 0.2, -0.1, 0.0, 0.9, 0.2, 1.0]
        ta_var.pars[:] = [0.01, 1.0, 0.0, 0.0, 0.2]
        ta_var.propagate_until(0.2)
        self._assert_finite(ta_var.state[:7])

    def test_zoh_eq_api(self):
        dyn = _pk.ta.zoh_eq_dyn()
        self.assertEqual(len(dyn), 7)

        ta = _pk.ta.get_zoh_eq(1e-12)
        self.assertEqual(len(ta.state), 7)
        self.assertEqual(len(ta.pars), 5)
        ta.time = 0.0
        ta.state[:] = [1.0, 0.01, 0.02, 0.01, -0.02, 0.3, 1.2]
        ta.pars[:] = [0.01, 0.0, 1.0, 0.0, 0.2]
        ta.propagate_until(0.2)
        self._assert_finite(ta.state)

        ta_var = _pk.ta.get_zoh_eq_var(1e-12)
        self.assertGreater(len(ta_var.state), 7)
        self.assertEqual(len(ta_var.pars), 5)
        ta_var.time = 0.0
        ta_var.state[:7] = [1.0, 0.01, 0.02, 0.01, -0.02, 0.3, 1.2]
        ta_var.pars[:] = [0.01, 0.0, 1.0, 0.0, 0.2]
        ta_var.propagate_until(0.2)
        self._assert_finite(ta_var.state[:7])

    def test_zoh_cr3bp_api(self):
        dyn = _pk.ta.zoh_cr3bp_dyn()
        self.assertEqual(len(dyn), 7)

        ta = _pk.ta.get_zoh_cr3bp(1e-12)
        self.assertEqual(len(ta.state), 7)
        self.assertEqual(len(ta.pars), 6)
        ta.time = 0.0
        ta.state[:] = [0.9, 0.0, 0.1, 0.0, 0.8, 0.0, 1.0]
        ta.pars[:] = [0.005, 1.0, 0.0, 0.0, 0.2, _pk.CR3BP_MU_EARTH_MOON]
        ta.propagate_until(0.2)
        self._assert_finite(ta.state)

        ta_var = _pk.ta.get_zoh_cr3bp_var(1e-12)
        self.assertGreater(len(ta_var.state), 7)
        self.assertEqual(len(ta_var.pars), 6)
        ta_var.time = 0.0
        ta_var.state[:7] = [0.9, 0.0, 0.1, 0.0, 0.8, 0.0, 1.0]
        ta_var.pars[:] = [0.005, 1.0, 0.0, 0.0, 0.2, _pk.CR3BP_MU_EARTH_MOON]
        ta_var.propagate_until(0.2)
        self._assert_finite(ta_var.state[:7])

    def test_zoh_ss_api(self):
        dyn = _pk.ta.zoh_ss_dyn()
        self.assertEqual(len(dyn), 6)

        ta = _pk.ta.get_zoh_ss(1e-12)
        self.assertEqual(len(ta.state), 6)
        self.assertEqual(len(ta.pars), 3)
        ta.time = 0.0
        ta.state[:] = [0.8, -0.4, 0.3, 0.2, 0.9, -0.1]
        ta.pars[:] = [0.25, -1.1, 0.04]
        ta.propagate_until(0.2)
        self._assert_finite(ta.state)

        ta_var = _pk.ta.get_zoh_ss_var(1e-12)
        self.assertGreater(len(ta_var.state), 6)
        self.assertEqual(len(ta_var.pars), 3)
        ta_var.time = 0.0
        ta_var.state[:6] = [0.8, -0.4, 0.3, 0.2, 0.9, -0.1]
        ta_var.pars[:] = [0.25, -1.1, 0.04]
        ta_var.propagate_until(0.2)
        self._assert_finite(ta_var.state[:6])

    def test_bcp_api(self):
        # dyn is time-dependent so no dyn accessor exposed, just check the integrator shape
        ta = _pk.ta.get_bcp(1e-12)
        self.assertEqual(len(ta.state), 6)
        self.assertEqual(len(ta.pars), 4)
        ta.time = 0.0
        ta.state[:] = [0.9, 0.0, 0.1, 0.0, 0.8, 0.0]
        ta.pars[:] = [_pk.CR3BP_MU_EARTH_MOON, 3.28900541e-5, 388.81, 0.925195985520347]
        ta.propagate_until(0.2)
        self._assert_finite(ta.state)

        ta_var = _pk.ta.get_bcp_var(1e-12)
        self.assertGreater(len(ta_var.state), 6)
        self.assertEqual(len(ta_var.pars), 4)
        ta_var.time = 0.0
        ta_var.state[:6] = [0.9, 0.0, 0.1, 0.0, 0.8, 0.0]
        ta_var.pars[:] = [_pk.CR3BP_MU_EARTH_MOON, 3.28900541e-5, 388.81, 0.925195985520347]
        ta_var.propagate_until(0.2)
        self._assert_finite(ta_var.state[:6])


class ta_regression_tests(_ut.TestCase):
    def test_zoh_ss_physical(self):
        # Physical regression test with GTOC13 reference trajectory (solar sail, ALTAIRA star units)
        ta = _pk.ta.get_zoh_ss(tol=1e-16)

        MU = 139348062043.343e9  # gravitational parameter of ALTAIRA star [m^3/s^2]
        L = 149597870.691e3     # reference radius (1 AU) [m]
        V = np.sqrt(MU / L)
        TIME = L / V
        ACC = V / TIME

        SAIL_C = 5.4026e-6   # sail radiation pressure coefficient [N/m^2]
        SAIL_A = 15000       # sail area [m^2]
        SAIL_MASS = 500      # spacecraft mass [kg]

        c = (2 * SAIL_C * SAIL_A / SAIL_MASS) / ACC
        ta.pars[2] = c
        ta.pars[:2] = [np.arctan(1 / np.sqrt(2)), 0.0]

        ic = [
            1.4959787069100000e+11 / L,
            0.0,
            0.0,
            0.0,
            3.0520227081462453e+04 / V,
            0.0,
        ]
        tof = 6.1595291581690468e+07 / TIME
        ta.time = 0.0
        ta.state[:] = ic
        ta.propagate_until(tof)

        gt = [
            1.0991008964199088e+11 / L,
            -1.0318907971048611e+11 / L,
            1.2103937089378202e+09 / L,
            2.0295741060391072e+04 / V,
            2.2479218998905460e+04 / V,
            -5.0463494709976425e+02 / V,
        ]
        self.assertTrue(np.allclose(gt, ta.state, rtol=1e-13))
