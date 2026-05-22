import unittest as _ut
import numpy as _np

import pykep as _pk

class cr3bp_udpla_tests(_ut.TestCase):
    def test_construction_and_eph_acc(self):
        # Build a simple non-dimensional reference state and construct udpla
        when = 6500.0  # MJD2000
        state_nd = _np.array([1.0, 0.0, 0.0, 0.0, 0.5, 0.0])
        mu = 0.0121505856  # example CR3BP mu (Earth-Moon like)
        TIME = 375000
        L = 384400000

        udpla = _pk.udpla.cr3bp(when, state_nd, mu, TIME, L, name="test_cr3bp", tol=1e-12)

        # Basic construction checks
        self.assertIsNotNone(udpla)
        # verify methods exist
        self.assertTrue(callable(getattr(udpla, "eph", None)))
        self.assertTrue(callable(getattr(udpla, "acc", None)))

    def test_eph_and_acc(self):
        # Build a simple non-dimensional reference state
        when = 6500.0  # MJD2000
        state_nd = _np.array([1.0, 0.0, 0.0, 0.0, 0.5, 0.0])
        mu = 0.0121505856  # example CR3BP mu (Earth-Moon like)
        TIME = 375000
        L = 384400000

        udpla = _pk.udpla.cr3bp(when, state_nd, mu, TIME, L, name="test_cr3bp", tol=1e-12)

        # eph at reference epoch should return the reference state scaled to SI units
        r_si, v_si = udpla.eph(when)
        # positions scaled by L, velocities scaled by V = L / TIME
        self.assertTrue(_np.allclose(r_si, state_nd[:3] * L))
        self.assertTrue(_np.allclose(v_si, state_nd[3:] * (L / TIME)))

        # acc should return a 3-element array in SI units
        acc = _np.asarray(udpla.acc(when))
        self.assertEqual(acc.shape, (3,))

        # regression: values should be stable across runs (within tolerance)
        expected_r = _np.array([1.0, 0.0, 0.0]) * L
        expected_v = _np.array([0.0, 0.5, 0.0]) * (L / TIME)

        self.assertTrue(_np.allclose(r_si, expected_r, atol=1e-14))
        self.assertTrue(_np.allclose(v_si, expected_v, atol=1e-14))

    def test_eph_acc_reproducibility(self):
        # Ensure eph/acc calls are reproducible: call at t=1, t=10, t=1 and
        # check that the results at t=1 are equal (within tolerance).
        when = 6500.0
        state_nd = _np.array([1.0809931218390707, 0.0, -0.20235953267405354, 0.0, -0.19895001215078018, 0.0])
        mu = 0.0121505856
        TIME = 375000
        L = 384400000

        udpla = _pk.udpla.cr3bp(when, state_nd, mu, TIME, L, name="test_cr3bp", tol=1e-12)

        r1, v1 = udpla.eph(0.1)
        a1 = _np.asarray(udpla.acc(0.1))

        # We now call the methods at 10. this causes the internal integrator state to go there and the next integratiopn to
        # start from there. We want to check that the results at 1 are the same as before, which verifies that the internal
        # state is not affecting the results.
        _ = udpla.eph(1.0)
        _ = udpla.acc(1.0)

        r1b, v1b = udpla.eph(0.1)
        a1b = _np.asarray(udpla.acc(0.1))

        self.assertTrue(_np.allclose(r1, r1b, atol=1e-14))
        self.assertTrue(_np.allclose(v1, v1b, atol=1e-14))
        self.assertTrue(_np.allclose(a1, a1b, atol=1e-14))
