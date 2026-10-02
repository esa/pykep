import pykep as _pk
import numpy as np

import unittest as _ut

class leg_zoh_ms_test(_ut.TestCase):
    def test_zoh_ms_set_initial_guess(self):
        """We test that the initial guess actually sets the interior nodes (closing the defects when
        not ballistic) without touching endpoints, tgrid and controls, and that malformed input raises."""
        ta = _pk.ta.get_zoh_kep(1e-12)
        state0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        controls = [0.0] * 8
        tgrid = [0.0, 0.5, 1.0]

        ta.state[:] = state0
        ta.pars[:4] = controls[:4]
        ta.propagate_until(tgrid[-1])
        state1 = ta.state.tolist()
        states = state0 + [10.0] * 7 + state1

        leg = _pk.leg.zoh_ms_py(states, controls, tgrid, 0.5, [ta, None])
        endpoint_states = states[:7] + states[-7:]
        original_tgrid = leg.tgrid.copy()

        leg.set_initial_guess()

        self.assertEqual(leg.states[:7] + leg.states[-7:], endpoint_states)
        self.assertEqual(leg.tgrid, original_tgrid)
        self.assertTrue(np.allclose(leg.compute_defects(), 0.0, atol=1e-10))

        ballistic_states = leg.states.copy()
        nonzero_controls = [0.1, 1.0, 0.0, 0.0] * 2
        leg.controls = nonzero_controls
        leg.states = state0 + [10.0] * 7 + state1
        leg.set_initial_guess(ballistic=True)

        self.assertEqual(leg.controls, nonzero_controls)
        self.assertTrue(np.allclose(leg.states[7:14], ballistic_states[7:14], atol=1e-12))

        leg.states = leg.states[:-1]
        with self.assertRaises(ValueError):
            leg.set_initial_guess()

        leg.states = states
        leg.tgrid = leg.tgrid[:-1]
        with self.assertRaises(ValueError):
            leg.compute_defects()

    def test_zoh_ms_defects_grad(self):
        """We test that the analytical sparse gradients agree with finite differences for any cut,
        and that the sparsity patterns are consistent (no duplicates, no missing nonzeros)."""
        ta = _pk.ta.get_zoh_kep(1e-14)
        ta_var = _pk.ta.get_zoh_kep_var(1e-14)
        ta.pars[4] = ta_var.pars[4] = 10.0
        rng = np.random.default_rng(42)
        nseg, d, h = 4, 7, 1e-7

        for cut in [0.0, 0.25, 0.5, 0.75, 1.0]:
            # Random leg with non-zero defects around a circular orbit
            node = [1.0, 0.1, 0.0, 0.0, 1.0, 0.1, 1.0]
            states = (np.tile(node, nseg + 1) + rng.uniform(-0.05, 0.05, d * (nseg + 1))).tolist()
            controls = rng.uniform(0, 0.01, 4 * nseg).tolist()
            tgrid = np.sort(rng.uniform(0, 2, nseg + 1)).tolist()
            leg = _pk.leg.zoh_ms_py(states, controls, tgrid, cut, [ta, ta_var])

            sparsities = leg.defects_grad_sparsity()
            grads = leg.compute_defects_grad()
            for name, sp, g in zip(["states", "controls", "tgrid"], sparsities, grads):
                self.assertEqual(len(np.unique(sp, axis=0)), len(sp))
                self.assertEqual(len(sp), len(g))

                # Dense finite-difference Jacobian of the defects w.r.t. this block
                x = list(getattr(leg, name))
                J_fd = np.zeros((d * nseg, len(x)))
                for j in range(len(x)):
                    xp, xm = list(x), list(x)
                    xp[j] += h
                    xm[j] -= h
                    setattr(leg, name, xp)
                    dp = np.array(leg.compute_defects())
                    setattr(leg, name, xm)
                    dm = np.array(leg.compute_defects())
                    J_fd[:, j] = (dp - dm) / 2 / h
                setattr(leg, name, x)

                J = np.zeros_like(J_fd)
                J[sp[:, 0], sp[:, 1]] = g
                self.assertTrue(np.allclose(J, J_fd, atol=1e-7))
