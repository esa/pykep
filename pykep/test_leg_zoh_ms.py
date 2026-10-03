import pykep as _pk
import numpy as np
import heyoka as _hy

import copy as _copy
import gc as _gc
import pickle as _pickle
import unittest as _ut


def _make_leg(leg_type, cut=0.5):
    """Construct a one-state leg with a non-control parameter on a nonuniform mesh."""
    # A small model checks that the code works outside the usual orbital dimensions.
    x = _hy.make_vars("x")
    ta = _hy.taylor_adaptive([(x, x + _hy.par[0] + _hy.par[1])], tol=1e-12)
    ta.time = 3.0
    ta.state[:] = [0.4]
    ta.pars[:] = [0.0, 0.7]
    args = ([1.0, 1.5, 2.0, 2.5], [0.1, -0.2, 0.3], [0.2, 0.4, 0.9, 1.3], cut, (ta, None))
    return leg_type(*args, dim_dynamics=1, dim_controls=1)


class leg_zoh_ms_test(_ut.TestCase):
    """Tests specific to the C++-backed zoh_ms class."""

    def test_properties(self):
        """We test that the public binding has matching defaults and validated properties."""
        # Mix arrays and tuples to check that the constructor accepts common input types.
        ta = _pk.ta.get_zoh_kep(1e-12)
        states = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0] * 3
        leg = _pk.leg.zoh_ms(
            states=np.array(states), controls=np.zeros(8), tgrid=(0.0, 0.5, 1.0), cut=0.5, tas=(ta, None)
        )
        # Check the public class name, stored inputs and default settings.
        self.assertIs(_pk.leg.zoh_ms, _pk.core._zoh_ms_cpp)
        self.assertEqual(_pk.leg.zoh_ms.__name__, "zoh_ms")
        self.assertEqual(_pk.leg.zoh_ms.__module__, "pykep.leg")
        self.assertEqual(leg.states, states)
        self.assertEqual(leg.controls, [0.0] * 8)
        self.assertEqual(leg.tgrid, [0.0, 0.5, 1.0])
        self.assertEqual((leg.dim_dynamics, leg.dim_controls), (7, 4))
        self.assertEqual((leg.nseg, leg.nseg_fwd, leg.nseg_bck), (2, 1, 1))
        self.assertIsNone(leg.max_steps)
        self.assertIsNone(leg.ta_var)
        # The C++ binding returns cached coordinate arrays without requiring variations.
        self.assertFalse(hasattr(leg, "compute_defects_grad"))
        sparsities = leg.defects_grad_sparsity()
        self.assertEqual([pattern.shape for pattern in sparsities], [(112, 2), (56, 2), (28, 2)])
        self.assertTrue(all(pattern.dtype == np.int64 for pattern in sparsities))
        np.testing.assert_array_equal(sparsities[0][:8], [[0, j] for j in range(7)] + [[0, 7]])
        np.testing.assert_array_equal(sparsities[1][:4], [[0, j] for j in range(4)])
        np.testing.assert_array_equal(sparsities[2][:2], [[0, 0], [0, 1]])
        leg.cut = 0.0
        np.testing.assert_array_equal(leg.defects_grad_sparsity()[0][:2], [[0, 0], [0, 7]])

        # Set each writable property, then read it back to check that it changed.
        leg = _make_leg(_pk.leg.zoh_ms)
        for name, value in (
            ("states", [2.0, 3.0, 4.0, 5.0]),
            ("controls", [0.2, 0.3, 0.4]),
            ("tgrid", [0.0, 0.1, 0.3, 0.6]),
            ("cut", 1.0),
            ("max_steps", 12),
        ):
            setattr(leg, name, value)
            self.assertEqual(getattr(leg, name), value)
        # Moving the cut to either endpoint should put all segments in one direction.
        self.assertEqual((leg.nseg_fwd, leg.nseg_bck), (3, 0))
        leg.cut = 0.0
        self.assertEqual((leg.nseg_fwd, leg.nseg_bck), (0, 3))
        leg.max_steps = None
        self.assertIsNone(leg.max_steps)
        # Editing a returned list must not silently change the leg's stored states.
        detached_states = leg.states
        detached_states[0] = 99.0
        self.assertEqual(leg.states[0], 2.0)

        # Readonly properties must reject assignment, even when the value is unchanged.
        for name in ("nseg", "nseg_fwd", "nseg_bck", "dim_dynamics", "dim_controls", "ta", "ta_var"):
            with self.subTest(readonly=name), self.assertRaises(AttributeError):
                setattr(leg, name, getattr(leg, name))
        # Invalid updates must raise an error and leave the previous data intact.
        for name, value in (("states", [1.0]), ("controls", [1.0]), ("tgrid", [0.0]), ("cut", 1.1)):
            original = getattr(leg, name)
            with self.subTest(invalid=name), self.assertRaises(RuntimeError):
                setattr(leg, name, value)
            self.assertEqual(getattr(leg, name), original)
        # A missing state component must also be rejected during construction.
        with self.assertRaises(RuntimeError):
            _pk.leg.zoh_ms(states[:-1], [0.0] * 8, [0.0, 0.5, 1.0], 0.5, (ta, None))

    def test_sampling_ownership_and_failure(self):
        """We test that samples outlive the leg and failed propagation reports incomplete histories."""
        leg = _make_leg(_pk.leg.zoh_ms)
        # Asking for zero samples is invalid.
        with self.assertRaises(RuntimeError):
            leg.get_state_info(N=0)
        # Save the samples, then delete the leg to check that the arrays keep their data.
        forward, backward, success = leg.get_state_info()
        snapshots = [segment.copy() for segment in forward + backward]
        self.assertTrue(success)
        del leg
        _gc.collect()
        for actual, expected in zip(forward + backward, snapshots):
            np.testing.assert_array_equal(actual, expected)
            actual[0, 0] += 1.0
            self.assertEqual(actual[0, 0], expected[0, 0] + 1.0)

        # An infinite starting state forces a failure after the first segment succeeds.
        leg = _make_leg(_pk.leg.zoh_ms, cut=1.0)
        leg.states = [1.0, float("inf"), 2.0, 2.5]
        forward, backward, success = leg.get_state_info(N=3)
        self.assertFalse(success)
        self.assertEqual(len(forward), 2)
        self.assertEqual(forward[0].shape, (3, 1))
        self.assertEqual(forward[1].shape, (1, 1))
        self.assertFalse(np.isfinite(forward[1][0, 0]))
        self.assertEqual(backward, [])

    def test_integrators_and_pickle(self):
        """We test that integrators are copied and copy/pickle preserve independent usable legs."""
        # Save the supplied integrators so we can check that the leg never changes them.
        ta = _pk.ta.get_zoh_kep(1e-12)
        ta_var = _pk.ta.get_zoh_kep_var(1e-12)
        state = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        ta.time = ta_var.time = 0.2
        ta.state[:] = ta_var.state[:7] = state
        ta.pars[:] = ta_var.pars[:] = [0.01, 1.0, 0.0, 0.0, 0.2]
        original = [(integrator.time, integrator.state.copy(), integrator.pars.copy()) for integrator in (ta, ta_var)]
        # Repeat with and without the optional integrator used for derivatives.
        for variational in (None, ta_var):
            with self.subTest(variational=variational is not None):
                leg = _pk.leg.zoh_ms(state * 3, [0.01, 1.0, 0.0, 0.0] * 2, [0.0, 0.5, 1.0], 0.5,
                                     (ta, variational), max_steps=100)
                self.assertIsNot(leg.ta, ta)
                if variational is not None:
                    self.assertIsNot(leg.ta_var, variational)
                    np.testing.assert_array_equal(leg.ta_var.state, variational.state)
                # Use the leg first, so we also test copies of a used integrator.
                leg.set_initial_guess()
                defects = leg.compute_defects()
                leg.get_state_info()
                self.assertTrue(repr(leg))
                # Normal copies, deep copies and saving/loading must preserve data and results.
                for clone in (_copy.copy(leg), _copy.deepcopy(leg), _pickle.loads(_pickle.dumps(leg))):
                    self.assertIsInstance(clone, _pk.leg.zoh_ms)
                    for name in ("states", "controls", "tgrid", "cut", "max_steps", "dim_dynamics", "dim_controls"):
                        self.assertEqual(getattr(clone, name), getattr(leg, name))
                    self.assertEqual(clone.ta_var is None, variational is None)
                    np.testing.assert_allclose(clone.compute_defects(), defects, atol=1e-12)
                    # Changing a copy must not change the original leg.
                    clone.states = [2.0] * len(clone.states)
                    self.assertNotEqual(clone.states, leg.states)
        # Compare against the saved inputs to detect any unwanted changes.
        for integrator, (time, state, pars) in zip((ta, ta_var), original):
            self.assertEqual(integrator.time, time)
            np.testing.assert_array_equal(integrator.state, state)
            np.testing.assert_array_equal(integrator.pars, pars)


class leg_zoh_ms_py_test(_ut.TestCase):
    """Tests specific to the pure-Python zoh_ms_py class, including gradients."""

    def test_copies_integrators(self):
        """We test that propagating a leg does not mutate supplied integrators."""
        # Build the integrators and save their time, states and parameters before use.
        ta = _pk.ta.get_zoh_kep(1e-12)
        ta_var = _pk.ta.get_zoh_kep_var(1e-12)
        state = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        states = state * 3
        controls = [0.01, 1.0, 0.0, 0.0] * 2
        tgrid = [0.0, 0.5, 1.0]
        pars = [*controls[:4], 0.2]

        ta.time = ta_var.time = 0.0
        ta.state[:] = state
        ta_var.state[:7] = state
        ta.pars[:] = ta_var.pars[:] = pars
        original = [
            (integrator.time, integrator.state.copy(), integrator.pars.copy())
            for integrator in (ta, ta_var)
        ]

        # Exercise propagation and derivatives, which use different integrator copies.
        leg = _pk.leg.zoh_ms_py(states, controls, tgrid, 0.5, [ta, ta_var])
        leg.compute_defects()
        leg.compute_defects_grad()

        # The supplied integrators must still match the saved versions exactly.
        self.assertIsNot(leg.ta, ta)
        self.assertIsNot(leg.ta_var, ta_var)
        for integrator, (time, state, parameters) in zip((ta, ta_var), original):
            self.assertEqual(integrator.time, time)
            self.assertTrue(np.array_equal(integrator.state, state))
            self.assertTrue(np.array_equal(integrator.pars, parameters))

    def test_set_initial_guess(self):
        """We test that the initial guess actually sets the interior nodes (closing the defects when
        not ballistic) without touching endpoints, tgrid and controls, and that malformed input raises."""
        ta = _pk.ta.get_zoh_kep(1e-12)
        state0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        controls = [0.0] * 8
        tgrid = [0.0, 0.5, 1.0]

        # Generate a reachable final state, then give the middle node a bad guess.
        ta.state[:] = state0
        ta.pars[:4] = controls[:4]
        ta.propagate_until(tgrid[-1])
        state1 = ta.state.tolist()
        states = state0 + [10.0] * 7 + state1

        leg = _pk.leg.zoh_ms_py(states, controls, tgrid, 0.5, [ta, None])
        endpoint_states = states[:7] + states[-7:]
        original_tgrid = leg.tgrid.copy()

        # Filling the middle node should close the defects without moving the endpoints.
        leg.set_initial_guess()

        self.assertEqual(leg.states[:7] + leg.states[-7:], endpoint_states)
        self.assertEqual(leg.tgrid, original_tgrid)
        self.assertTrue(np.allclose(leg.compute_defects(), 0.0, atol=1e-10))

        # Ballistic guessing ignores thrust but must keep the stored controls unchanged.
        ballistic_states = leg.states.copy()
        nonzero_controls = [0.1, 1.0, 0.0, 0.0] * 2
        leg.controls = nonzero_controls
        leg.states = state0 + [10.0] * 7 + state1
        leg.set_initial_guess(ballistic=True)

        self.assertEqual(leg.controls, nonzero_controls)
        self.assertTrue(np.allclose(leg.states[7:14], ballistic_states[7:14], atol=1e-12))

        # Python accepts these assignments, but must reject bad sizes before computing.
        leg.states = leg.states[:-1]
        with self.assertRaises(ValueError):
            leg.set_initial_guess()

        leg.states = states
        leg.tgrid = leg.tgrid[:-1]
        with self.assertRaises(ValueError):
            leg.compute_defects()

    def test_defects_grad(self):
        """We test that the analytical sparse gradients agree with finite differences for any cut,
        and that the sparsity patterns are consistent (no duplicates, no missing nonzeros)."""
        ta = _pk.ta.get_zoh_kep(1e-14)
        ta_var = _pk.ta.get_zoh_kep_var(1e-14)
        ta.pars[4] = ta_var.pars[4] = 10.0
        # A fixed random seed gives varied inputs while keeping the test repeatable.
        rng = np.random.default_rng(42)
        nseg, d, h = 4, 7, 1e-7

        # Try several forward/backward splits, including both endpoint cases.
        for cut in [0.0, 0.25, 0.5, 0.75, 1.0]:
            # Small random offsets from a circular orbit give nonzero defects.
            node = [1.0, 0.1, 0.0, 0.0, 1.0, 0.1, 1.0]
            states = (np.tile(node, nseg + 1) + rng.uniform(-0.05, 0.05, d * (nseg + 1))).tolist()
            controls = rng.uniform(0, 0.01, 4 * nseg).tolist()
            tgrid = np.sort(rng.uniform(0, 2, nseg + 1)).tolist()
            leg = _pk.leg.zoh_ms_py(states, controls, tgrid, cut, [ta, ta_var])

            sparsities = leg.defects_grad_sparsity()
            grads = leg.compute_defects_grad()
            # Each listed matrix entry must be unique and have one derivative value.
            for name, sp, g in zip(["states", "controls", "tgrid"], sparsities, grads):
                self.assertEqual(len(np.unique(sp, axis=0)), len(sp))
                self.assertEqual(len(sp), len(g))

                # Estimate each derivative using small positive and negative input changes.
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

                # Put the sparse values into a full matrix and compare with our estimates.
                J = np.zeros_like(J_fd)
                J[sp[:, 0], sp[:, 1]] = g
                self.assertTrue(np.allclose(J, J_fd, atol=1e-7))


class leg_zoh_ms_api_test(_ut.TestCase):
    """Tests comparing the shared non-gradient APIs and results of both implementations."""

    def test_shared_api(self):
        """We test that both classes expose the same shared properties and callable methods."""
        # Give both implementations the same inputs before comparing their interfaces.
        cpp = _make_leg(_pk.leg.zoh_ms)
        python = _make_leg(_pk.leg.zoh_ms_py)
        # Shared properties should have the same Python types and initial values.
        for name in ("states", "controls", "tgrid", "cut", "max_steps", "dim_dynamics", "dim_controls",
                     "nseg", "nseg_fwd", "nseg_bck"):
            with self.subTest(property=name):
                self.assertIs(type(getattr(cpp, name)), type(getattr(python, name)))
                self.assertEqual(getattr(cpp, name), getattr(python, name))
            # Check the integrator data too, not just the leg's mesh and settings.
        self.assertIs(type(cpp.ta), type(python.ta))
        np.testing.assert_array_equal(cpp.ta.state, python.ta.state)
        np.testing.assert_array_equal(cpp.ta.pars, python.ta.pars)
        self.assertIsNone(cpp.ta_var)
        self.assertIsNone(python.ta_var)
        # Both classes must provide the shared methods; gradients are tested separately.
        for name in ("compute_defects", "set_initial_guess", "get_state_info"):
            self.assertTrue(callable(getattr(cpp, name)))
            self.assertTrue(callable(getattr(python, name)))
        # Apply the same valid updates and check that the results still agree.
        for name, value in (("states", [2.0, 3.0, 4.0, 5.0]), ("controls", [0.2, 0.3, 0.4]),
                            ("tgrid", [0.0, 0.1, 0.3, 0.6]), ("max_steps", 100)):
            setattr(cpp, name, value)
            setattr(python, name, value)
            self.assertEqual(getattr(cpp, name), getattr(python, name))
        np.testing.assert_allclose(cpp.compute_defects(), python.compute_defects(), rtol=1e-11, atol=1e-12)

    def test_results(self):
        """We test that defects, initial guesses and sampled trajectories agree across cuts."""
        # Test backward-only, mixed and forward-only propagation.
        for cut in (0.0, 0.5, 1.0):
            with self.subTest(cut=cut):
                cpp = _make_leg(_pk.leg.zoh_ms, cut)
                python = _make_leg(_pk.leg.zoh_ms_py, cut)
                # Compare both the defect list layout and its numerical values.
                self.assertIsInstance(cpp.compute_defects(), list)
                self.assertIsInstance(python.compute_defects(), list)
                self.assertEqual(len(cpp.compute_defects()), 3)
                self.assertEqual(len(python.compute_defects()), 3)
                np.testing.assert_allclose(cpp.compute_defects(), python.compute_defects(), rtol=1e-11, atol=1e-12)
                # Check the default sample count, a custom count and a single starting sample.
                for sample_count in (None, 3, 1):
                    cpp_info = cpp.get_state_info() if sample_count is None else cpp.get_state_info(N=sample_count)
                    py_info = python.get_state_info() if sample_count is None else python.get_state_info(N=sample_count)
                    self.assertIsInstance(cpp_info, tuple)
                    self.assertIsInstance(py_info, tuple)
                    self.assertEqual(len(cpp_info), 3)
                    self.assertEqual(len(py_info), 3)
                    self.assertIs(cpp_info[2], True)
                    self.assertIs(py_info[2], True)
                    # Compare segment counts, array shapes and values in each direction.
                    for actual, expected, count in zip(cpp_info[:2], py_info[:2], (cpp.nseg_fwd, cpp.nseg_bck)):
                        self.assertIsInstance(actual, list)
                        self.assertIsInstance(expected, list)
                        self.assertEqual(len(actual), count)
                        self.assertEqual(len(expected), count)
                        for segment, reference in zip(actual, expected):
                            self.assertIsInstance(segment, np.ndarray)
                            self.assertIsInstance(reference, np.ndarray)
                            self.assertEqual(segment.dtype, np.dtype("float64"))
                            self.assertEqual(segment.dtype, reference.dtype)
                            self.assertEqual(segment.shape, (5 if sample_count is None else sample_count, 1))
                            self.assertEqual(segment.shape, reference.shape)
                            np.testing.assert_allclose(segment, reference, rtol=1e-11, atol=1e-12)

                # Compare initial guesses made with the stored controls and with zero controls.
                for ballistic in (False, True):
                    cpp = _make_leg(_pk.leg.zoh_ms, cut)
                    python = _make_leg(_pk.leg.zoh_ms_py, cut)
                    original = (cpp.states, cpp.controls, cpp.tgrid)
                    if ballistic:
                        self.assertIsNone(cpp.set_initial_guess(ballistic=True))
                        self.assertIsNone(python.set_initial_guess(ballistic=True))
                    else:
                        self.assertIsNone(cpp.set_initial_guess())
                        self.assertIsNone(python.set_initial_guess())
                    np.testing.assert_allclose(cpp.states, python.states, rtol=1e-11, atol=1e-12)
                    # Only interior nodes may change; endpoints, controls and times must stay fixed.
                    for leg in (cpp, python):
                        self.assertEqual(leg.states[0], original[0][0])
                        self.assertEqual(leg.states[-1], original[0][-1])
                        self.assertEqual(leg.controls, original[1])
                        self.assertEqual(leg.tgrid, original[2])
                    np.testing.assert_allclose(cpp.compute_defects(), python.compute_defects(), rtol=1e-11, atol=1e-12)
