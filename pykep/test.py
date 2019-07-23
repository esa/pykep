# -*- coding: utf-8 -*-

# Copyright 2017-2018 PaGMO development team
#
# This file is part of the PaGMO library.
#
# The PaGMO library is free software; you can redistribute it and/or modify
# it under the terms of either:
#
#   * the GNU Lesser General Public License as published by the Free
#     Software Foundation; either version 3 of the License, or (at your
#     option) any later version.
#
# or
#
#   * the GNU General Public License as published by the Free Software
#     Foundation; either version 3 of the License, or (at your option) any
#     later version.
#
# or both in parallel, as here.
#
# The PaGMO library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received copies of the GNU General Public License and the
# GNU Lesser General Public License along with the PaGMO library.  If not,
# see https://www.gnu.org/licenses/.

from __future__ import absolute_import as _ai

import unittest as _ut
import numpy as np


class core_functions_test_case(_ut.TestCase):
    """Test case for the core functions

    """

    def runTest(self):
        self.run_epoch_test()

    def run_epoch_test(self):
        from .core import epoch
        epoch(julian_date=0, julian_date_type='mjd2000')


class lambert_test_case(_ut.TestCase):
    """Test case for the lambert problem class

    """

    def runTest(self):
        from .core import lambert_problem
        lambert_problem(r1=[1, 1, 0], r2=[0, 1, 0],
                        mu=1., cw=False, max_revs=0, tof=0.3)


class mga_1dsm_test_case(_ut.TestCase):
    """Test case for the mga1_dsm class

    """

    def runTest(self):
        self.run_construction_test()
        self.run_decode_times_and_vinf_test()

    def run_construction_test(self):
        from .trajopt import mga_1dsm
        # Correct use (nothrow)
        udp = mga_1dsm()
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='direct', tof=[[20, 400], [20, 400]])
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='eta', tof=500)
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='alpha', tof=[20, 500])
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='direct', tof=[
                       [20, 400], [20, 400]], multi_objective=True)
        self.assertEqual(udp.get_nobj(), 2)

        # Incorrect use (raise)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=34)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=[[400], [20, 400]])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=[20, 400])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='eta', tof=[20, 400])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='eta', tof=[[20, 400], [20, 400]])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='alpha', tof=4)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='alpha', tof=[[20, 400], [20, 400]])

    def run_decode_times_and_vinf_test(self):
        from .trajopt import mga_1dsm
        # Alpha
        udp = mga_1dsm(tof_encoding='alpha', tof=[20, 500])
        x = [10] + [0., 0., 1., 1., 0.5] + [0.4, 1.3, 0.5, 0.5] + [50]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 25.)
        self.assertAlmostEqual(retval[0][1], 25.)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)
        # Eta
        udp = mga_1dsm(tof_encoding='eta', tof=50)
        x = [10] + [0., 0., 1., 1., 0.5] + [0.4, 1.3, 0.5, 0.5]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 25.)
        self.assertAlmostEqual(retval[0][1], 12.5)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)
        # Direct
        udp = mga_1dsm(tof_encoding='direct', tof=[[10, 400], [10, 400]])
        x = [10] + [0., 0., 1., 1., 123] + [0.4, 1.3, 0.5, 321]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 123)
        self.assertAlmostEqual(retval[0][1], 321)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)


class gym_test_case(_ut.TestCase):
    """Test case for the gym

    """

    def runTest(self):
        self.run_rosetta_test()
        self.run_cassini2_test()
        self.run_tandem_test()
        self.run_juice_test()
        self.run_messenger_test()

    def run_rosetta_test(self):
        from .trajopt import gym
        udp = gym.rosetta
        x = [1.53488329e+03, 4.56388378e-01, 9.51717655e-01, 4.18212047e+03,
             4.32159299e-01, 3.65256539e+02, 5.03363275e+00, 2.38949977e+00,
             4.55746823e-01, 7.09999954e+02, 1.79894273e+00, 1.05000003e+00,
             6.09083347e-01, 2.60816142e+02, 4.95158968e+00, 3.16049580e+00,
             6.89049263e-01, 7.29775762e+02, 4.30823655e+00, 1.10842692e+00,
             4.16075410e-01, 1.84999995e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], 1371.4992595125018)

    def run_cassini2_test(self):
        from .trajopt import gym
        udp = gym.cassini2
        x = [-7.75699976e+02,  9.15777367e-01,  4.06442043e-01,  3.21309562e+03,
        6.81118341e-01,  1.62660490e+02, -1.58051063e+00,  1.28479507e+00,
        4.72699902e-01,  4.24319550e+02,  4.30475919e+00,  1.15739933e+00,
        2.55718252e-01,  5.44489098e+01, -1.54332794e+00,  1.27160729e+00,
        9.00000000e-01,  5.88481599e+02,  4.76774269e+00,  7.00000000e+01,
        1.00000000e-02,  2.20000000e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], 1511.731793212)

    def run_tandem_test(self):
        from .trajopt import gym
        udp = gym.tandem(prob_id = 6, constrained=False)
        x = [ 7.98593791e+03,  9.02569236e-01,  3.59510662e-01,  3.31011777e+03,
        1.00000000e-02,  1.68854894e+02, -1.15321363e+00,  1.05006538e+00,
        4.21911761e-01,  3.38579352e+02, -1.68143357e+00,  1.10121434e+00,
        6.39679220e-01,  1.17925136e+03, -2.02838607e+00,  1.05000000e+00,
        2.31554532e-01,  1.80000707e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], -7.302117213470749)

    def run_juice_test(self):
        from .trajopt import gym
        udp = gym.juice
        x = [ 8.16283083e+03,  6.41922787e-01,  6.51202691e-01,  2.51009414e+03,
        2.97841478e-01,  3.81541370e+02,  9.58572190e-01,  1.53007674e+00,
        3.06125365e-01,  1.49264351e+02,  4.10103546e+00,  2.39297670e+00,
        4.34424957e-01,  3.16066418e+02,  4.33225338e+00,  1.30946367e+00,
        4.52048883e-01,  1.63208108e+02,  5.93850330e-01,  1.34871269e+00,
        2.03288502e-01,  6.52494606e+02, -1.37902374e+00,  1.55482534e+00,
        1.96917559e-01,  1.08471126e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], -7.987614956531155)

    def run_messenger_test(self):
        from .trajopt import gym
        udp = gym.messenger
        x = [ 2.03241398e+03,  6.40762059e-01,  6.63357785e-01,  4.04989271e+03,
        6.63732323e-01,  4.50068524e+02, -3.86553343e+00,  3.52631372e+00,
        5.57888828e-01,  2.24619580e+02, -4.45910441e+00,  1.22736521e+00,
        7.08063036e-01,  2.17965497e+02, -2.47894274e+00,  1.43586128e+00,
        5.88391838e-01,  2.62423586e+02, -2.40594385e-02,  2.45470457e+00,
        7.25370468e-01,  3.58067954e+02,  1.47192632e+00,  1.05000000e+00,
        9.02984391e-01,  5.38436770e+02]
        self.assertAlmostEqual(udp.fitness(x)[0], 5855.8143406005165)


class spherical_harmonics_loader_test_case(_ut.TestCase):
    """Test case for the spherical harmonic gravity file loader function.

    """

    def runTest(self):
        self.run_egm96_test()
        self.run_ggm02c_test()
        self.run_ggm02s_test()
        self.run_jgmro120d_test()
        self.run_glgm3150_test()
        self.run_jgl150q1_test()
        self.run_lpe200_test()
        
    def run_egm96_test(self):
        from .util import load_gravity_model

        req_r = 6378137
        req_mu = 3.986004418e+14

        c20_20 = 4.01448327968e-09
        c100_100 = 1.10930637955e-09

        s100_50 = -1.06362863541e-09
        s250_125 = -8.34356616626e-11

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Earth/egm96.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[250, 125], s250_125)

    def run_ggm02c_test(self):
        from .util import load_gravity_model

        req_r = 6378136.3
        req_mu = 3.986004415e+14

        c20_20 = 3.734698472368e-09
        c100_100 = 8.8043569591782e-10

        s100_50 = -8.2691933615552e-10
        s200_100 = -1.3266576400742e-12

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Earth/ggm02c.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[200, 100], s200_100)

    def run_ggm02s_test(self):
        from .util import load_gravity_model

        req_r = 6378136.3
        req_mu = 3.986004415e+14

        c20_20 = 3.7323168217059e-09
        c100_100 = 1.0356187365035e-09

        s100_50 = -8.3581312522994e-10
        s160_80 = -2.327147027496e-09

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Earth/ggm02s.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[160, 80], s160_80)

    def run_jgmro120d_test(self):
        from .util import load_gravity_model

        req_r = 3396000.0
        req_mu = 4.28283758157561e+13

        c20_20 = -4.767064931643e-07
        c100_100 = -6.238892373683e-09

        s100_50 = -1.191278972881e-08
        s120_60 = 5.598027946471e-10

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Mars/jgmro_120d.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[120, 60], s120_60)

    def run_glgm3150_test(self):
        from .util import load_gravity_model

        req_r = 1738000.0
        req_mu = 4.9002800238e+12

        c20_20 = -1.4913227206107e-08
        c100_100 = -3.0367607155944e-08

        s100_50 = 1.7233711523911e-10
        s120_60 = -5.3052488170131e-09

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Moon/glgm3_150.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[120, 60], s120_60)

    def run_jgl150q1_test(self):
        from .util import load_gravity_model

        req_r = 1738000.0
        req_mu = 4.902801076e+12

        c20_20 = 1.17705963027e-07
        c100_100 = -1.5618538762e-08

        s100_50 = 5.1145625569e-09
        s120_60 = -1.05725130139e-09

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Moon/jgl_150q1.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[120, 60], s120_60)

    def run_lpe200_test(self):
        from .util import load_gravity_model

        req_r = 1738000.0
        req_mu = 4.902800238e+12

        c20_20 = -1.5489032082686e-08
        c100_100 = -3.0841142975112e-08

        s100_50 = -1.53194061215089e-09
        s200_100 = 4.25825545837e-09

        r, mu, c, s, n, m = load_gravity_model("pykep/util/gravity_models/Moon/lpe_200.txt")

        self.assertEqual(r, req_r)
        self.assertEqual(mu, req_mu)
        
        self.assertEqual(c[20, 20], c20_20)
        self.assertEqual(c[100, 100], c100_100)

        self.assertEqual(s[100, 50], s100_50)
        self.assertEqual(s[200, 100], s200_100)

class gravity_spherical_harmonic_test_case(_ut.TestCase):
    """Test case for the spherical harmonic gravity calculator function.

    """

    def runTest(self):
        from .util import load_gravity_model, gravity_spherical_harmonic

        r, mu, c, s, n_max, m_max = load_gravity_model("pykep/util/gravity_models/Earth/egm96.txt")

        x = np.array([[6.07303362e+06, -1.63535914e-9, -3.22908926e+06],
                      [-5874145.34596, 1745831.60905, 3123338.4834],
                      [5290507.45841, -3377313.30177, -2813012.71391],
                      [-4360347.55688, 4787584.94679, 2318437.92131],
                      [3144590.02052, -5884275.41062, -1672008.17262]])

        acc = gravity_spherical_harmonic(x, r, mu, c, s, n_max, m_max)

        req_acc = np.array([[-7.438268885207450, 4.174587578722027e-05, 3.966055730251360],
                            [7.195259383674340, -2.138389507954604, -3.836583575161607],
                            [-6.482132055071062, 4.138055191303929, 3.456233485081509],
                            [5.344543324379746, -5.868302782861193, -2.849822685569829],
                            [-3.855999728216701, 7.215225587047862, 2.055872629061049]])

        for a, req_a in zip(acc, req_acc):
            for ax, req_ax in zip(a, req_a):
                self.assertAlmostEqual(ax, req_ax, 13)


def run_test_suite(level=0):
    """Run the full test suite.

    This function will raise an exception if at least one test fails.

    Args:
        level(``int``): the test level (higher values run longer tests)

    """

    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(core_functions_test_case)
    suite.addTest(lambert_test_case())
    suite.addTest(mga_1dsm_test_case())
    suite.addTest(gym_test_case())
    suite.addTest(spherical_harmonics_loader_test_case())
    suite.addTest(gravity_spherical_harmonic_test_case())


    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')
