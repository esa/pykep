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



class trajopt_mga_tests(_ut.TestCase):
    def test_construction(self):
        import pykep as pk

        earth = pk.planet(pk.udpla.jpl_lp("earth"))
        venus = pk.planet(pk.udpla.jpl_lp("venus"))
        udp = pk.trajopt.mga(
            seq=[earth, venus, earth, venus, earth],
            tof_encoding="direct",
            t0=[0, 1000],
            tof=[[30, 200], [30, 300], [30, 300], [30, 300]],
            vinf=2.5,
        )
        prob = pg.problem(udp)
        pop = pg.population(prob, 100)

    def test_encoding_to_encoding(self):
        import pykep as pk

        udp_direct = pk.trajopt.mga(tof_encoding="direct", tof=[[30, 200], [200, 300]])
        udp_alpha = pk.trajopt.mga(tof_encoding="alpha", tof=[230, 500])
        udp_eta = pk.trajopt.mga(tof_encoding="eta", tof=500)
        prob = pg.problem(udp_direct)
        pop = pg.population(prob, 100)
        x_direct = pop.champion_x
        gt = udp_direct.fitness(x_direct)[0]
        x_alpha = udp_alpha.direct2alpha(x_direct)
        x_eta = udp_eta.direct2eta(x_direct)
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-12)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-12)
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0])
            < 1e-12
        )
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta))[0])
            < 1e-12
        )


class gym_tests(_ut.TestCase):
    def test_cassini1(self):
        import pykep as pk

        udp = pk.trajopt.gym.cassini1
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            -554.5189290104555,
            103.27184879471751,
            335.41655259663474,
            80.50258543604521,
            862.0950563689543,
            2865.018040480413,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 80400.08898184073) < 1e-12)

        x = [
            -932.0532394941108,
            37.534681289972674,
            162.5093144821548,
            336.970139545233,
            1743.2915882004586,
            2527.8785180526465,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 217105.875031613573) < 1e-12)

        x = [
            -583.0776694390058,
            388.65047998036107,
            334.9959782156864,
            65.57508619540917,
            1520.2982946551908,
            2132.7771932619144,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 107218.08496509642) < 1e-12)

    def test_cassini2(self):
        import pykep as pk

        udp = pk.trajopt.gym.cassini2
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            -7.75699976e02,
            9.15777367e-01,
            4.06442043e-01,
            3.21309562e03,
            6.81118341e-01,
            1.62660490e02,
            -1.58051063e00,
            1.28479507e00,
            4.72699902e-01,
            4.24319550e02,
            4.30475919e00,
            1.15739933e00,
            2.55718252e-01,
            5.44489098e01,
            -1.54332794e00,
            1.27160729e00,
            9.00000000e-01,
            5.88481599e02,
            4.76774269e00,
            7.00000000e01,
            1.00000000e-02,
            2.20000000e03,
        ]

        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1511.7317645968126) < 1e-12)

    def test_rosetta(self):
        import pykep as pk

        udp = pk.trajopt.gym.rosetta
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            1.53488329e03,
            4.56388378e-01,
            9.51717655e-01,
            4.18212047e03,
            4.32159299e-01,
            3.65256539e02,
            5.03363275e00,
            2.38949977e00,
            4.55746823e-01,
            7.09999954e02,
            1.79894273e00,
            1.05000003e00,
            6.09083347e-01,
            2.60816142e02,
            4.95158968e00,
            3.16049580e00,
            6.89049263e-01,
            7.29775762e02,
            4.30823655e00,
            1.10842692e00,
            4.16075410e-01,
            1.84999995e03,
        ]

        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1371.4992633334382) < 1e-12)

    def test_eve_mga1dsm(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            7.31864730e02,
            6.62420829e-01,
            3.46714249e-01,
            1.60589872e03,
            4.69582854e-01,
            2.98596210e02,
            -1.90774304e00,
            2.05710242e01,
            3.63127164e-01,
            9.95555392e01,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 47456.939061940415) < 1e-12)

    def test_eve_mga1dsm_a(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_a
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            1.12861301e02,
            5.49233732e-01,
            3.04597487e-02,
            1.92554472e03,
            5.22619135e-01,
            8.46696560e-01,
            -2.64317289e00,
            2.16924824e01,
            6.62441172e-01,
            8.89339812e-02,
            5.15383086e02,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1101622.7179572878) < 1e-12)

    def test_eve_mga1dsm_n(self):
        import pykep as pk

        udp = pk.trajopt.gym.eve_mga1dsm_n
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            5.86500918e02,
            6.32855532e-01,
            9.59033298e-01,
            1.93800759e03,
            8.02901287e-01,
            9.88679911e-01,
            5.18276555e00,
            1.04655908e01,
            5.84524787e-01,
            9.68549775e-01,
        ]
        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 1917650.9004062244) < 1e-12)

    def test_juice(self):
        import pykep as pk

        udp = pk.trajopt.gym.juice
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            8.25587848945082e03,
            3.87927696743959e-01,
            5.27687480720340e-01,
            2.06516101304215e03,
            7.74309396597230e-01,
            4.66381738991143e02,
            -1.02751095385811e00,
            8.25250602985045e00,
            1.85241667091363e-01,
            1.10839235700056e02,
            3.89240251780202e00,
            2.50359555867057e00,
            6.46178864367641e-02,
            2.15407037751815e02,
            4.78686007014207e00,
            8.41382633987924e00,
            3.38707892618604e-01,
            3.19871173483077e01,
            2.06974341215216e00,
            4.36373629930523e00,
            4.66997711732296e-01,
            7.22372043742636e02,
            4.77835192401833e00,
            7.65391290702327e00,
            4.89771110016942e-01,
            1.03158288517734e03,
        ]

        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, 204.38446462495546) < 1e-12)

        # Solution previously in the old pykep gym code
        x = [
            8.16283083e03,
            6.41922787e-01,
            6.51202691e-01,
            2.51009414e03,
            2.97841478e-01,
            3.81541370e02,
            9.58572190e-01,
            1.53007674e00,
            3.06125365e-01,
            1.49264351e02,
            4.10103546e00,
            2.39297670e00,
            4.34424957e-01,
            3.16066418e02,
            4.33225338e00,
            1.30946367e00,
            4.52048883e-01,
            1.63208108e02,
            5.93850330e-01,
            1.34871269e00,
            2.03288502e-01,
            6.52494606e02,
            -1.37902374e00,
            1.55482534e00,
            1.96917559e-01,
            1.08471126e03,
        ]

        f = udp.fitness(x)[0]
        self.assertTrue(float_rel_error(f, -7.987614927397559) < 1e-12)

    def test_juice_mo(self):
        import pykep as pk

        udp = pk.trajopt.gym.juice_mo
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            8.08458777709711e03,
            3.98862191419023e-02,
            1.13145851845895e-01,
            2.54109146808688e03,
            7.47321603643624e-01,
            2.27778122138468e-02,
            -4.41107070349122e00,
            6.27932399248486e00,
            7.64033667254947e-01,
            2.94097311318708e-01,
            2.72683135802239e00,
            3.48208267655823e00,
            9.73171666093481e-02,
            5.02201106930738e-01,
            3.76552737505492e00,
            5.97748097683895e00,
            5.76935286152452e-01,
            1.52899955890002e-01,
            -1.54719587734854e00,
            8.37088373080571e00,
            7.83924298935545e-01,
            5.32290337033144e-01,
            -5.62028545284311e00,
            3.37199051531738e00,
            9.13407577042211e-01,
            2.61007823689634e-01,
            2.27004720198764e03,
        ]

        f = udp.fitness(x)
        self.assertTrue(float_rel_error(f[0], 427.31998557200325) < 1e-12)
        self.assertTrue(float_rel_error(f[1], 2270.0472019876433) < 1e-12)

    def test_messenger(self):
        import pykep as pk

        udp = pk.trajopt.gym.messenger
        # Ground truth checked by the old pykep code (up to 8 digits only as per differences with constants and all)
        x = [
            2.03241398e03,
            6.40762059e-01,
            6.63357785e-01,
            4.04989271e03,
            6.63732323e-01,
            4.50068524e02,
            -3.86553343e00,
            3.52631372e00,
            5.57888828e-01,
            2.24619580e02,
            -4.45910441e00,
            1.22736521e00,
            7.08063036e-01,
            2.17965497e02,
            -2.47894274e00,
            1.43586128e00,
            5.88391838e-01,
            2.62423586e02,
            -2.40594385e-02,
            2.45470457e00,
            7.25370468e-01,
            3.58067954e02,
            1.47192632e00,
            1.05000000e00,
            9.02984391e-01,
            5.38436770e02,
        ]

        f = udp.fitness(x)
        self.assertTrue(float_rel_error(f[0], 5855.81434335236) < 1e-12)


class trajopt_mga1dsm_tests(_ut.TestCase):
    def test_construction(self):
        import pykep as pk

        earth = pk.planet(pk.udpla.jpl_lp("earth"))
        venus = pk.planet(pk.udpla.jpl_lp("venus"))
        udp = pk.trajopt.mga_1dsm(
            seq=[
                earth,
                venus,
                earth,
            ],
            t0=[0, 1000],
            tof=[[30, 200], [200, 300]],
            vinf=[0.5, 2.5],
            add_vinf_dep=False,
            add_vinf_arr=True,
            tof_encoding="direct",
            multi_objective=False,
            orbit_insertion=False,
            e_target=None,
            rp_target=None,
            eta_bounds=[0.1, 0.9],
            rp_ub=30,
        )
        prob = pg.problem(udp)
        pop = pg.population(prob, 100)

    def test_encoding_to_encoding(self):
        import pykep as pk

        udp_direct = pk.trajopt.mga_1dsm(
            tof_encoding="direct", tof=[[30, 200], [200, 300]]
        )
        udp_alpha = pk.trajopt.mga_1dsm(tof_encoding="alpha", tof=[230, 500])
        udp_eta = pk.trajopt.mga_1dsm(tof_encoding="eta", tof=500)
        prob = pg.problem(udp_direct)
        pop = pg.population(prob, 100)
        x_direct = pop.champion_x
        gt = udp_direct.fitness(x_direct)[0]
        x_alpha = udp_alpha.direct2alpha(x_direct)
        x_eta = udp_eta.direct2eta(x_direct, 500)
        self.assertTrue(float_rel_error(gt, udp_alpha.fitness(x_alpha)[0]) < 1e-12)
        self.assertTrue(float_rel_error(gt, udp_eta.fitness(x_eta)[0]) < 1e-12)
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_direct.alpha2direct(x_alpha))[0])
            < 1e-12
        )
        self.assertTrue(
            float_rel_error(gt, udp_direct.fitness(udp_eta.eta2direct(x_eta, 500))[0])
            < 1e-12
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.alpha2direct(x_alpha) - x_direct) < 1e-12
        )
        self.assertTrue(
            np.linalg.norm(udp_direct.eta2direct(x_eta, 500) - x_direct) < 1e-12
        )


class mit_tests(_ut.TestCase):
    def test_primer_vector(self):
        import pykep as pk

        DVi = np.random.random((3,)) / 2 - 1
        DVj = np.random.random((3,)) / 2 - 1

        # Test the primer vector at k=i
        Mji = np.random.random((6, 6)) / 2 - 1
        Mjk = Mji
        p, Aik, Ajk = pk.trajopt.primer_vector(DVi, DVj, Mji, Mjk)
        self.assertTrue(float_rel_error(np.linalg.norm(p), 1) < 1e-12)
        self.assertTrue(
            float_rel_error(np.linalg.norm(Aik), np.linalg.norm(np.eye(3))) < 1e-12
        )
        self.assertTrue(abs(np.linalg.norm(Ajk)) < 1e-12)
        
        # Test the primer vector at k=j
        Mji = np.random.random((6, 6)) / 2 - 1
        Mjk = np.eye(6)
        p, Aik, Ajk = pk.trajopt.primer_vector(DVi, DVj, Mji, Mjk)
        self.assertTrue(float_rel_error(np.linalg.norm(p), 1) < 1e-12)
        self.assertTrue(
            float_rel_error(np.linalg.norm(Ajk), np.linalg.norm(np.eye(3))) < 1e-12
        )
        self.assertTrue(abs(np.linalg.norm(Aik)) < 1e-12)
        
        # Test at a generic case
        DVi = np.array([0.87066204, 0.08769792, 0.83858668])
        DVj = np.array([0.89810309, 0.69398023, 0.79896632])
        Mji = np.array([
            [ 0.72007593,  0.30196707, -0.43032973, -0.73575787,  0.5111975 , -0.5574713 ],
            [-0.89571788, -0.08538045,  0.50266506, -0.04080879,  0.34533095, 0.125654  ],
            [-0.85658008,  0.57145651,  0.64125294,  0.23238585,  0.44306253, 0.09536816],
            [-0.92263004,  0.84154659,  0.16484889,  0.80878073,  0.73393498, 0.20910764],
            [-0.31083144,  0.04295479, -0.89515654, -0.88093621,  0.56141158, 0.58555147],
            [-0.591202  , -0.97773281,  0.21925006,  0.74073955, -0.02260352, -0.55315306]
        ])
        Mjk = np.array([
            [ 0.23272111,  0.5072838 ,  0.85937303,  0.22824718,  0.97333937, 0.09703637],
            [ 0.70277786,  0.1919488 , -0.62527082, -0.01379391, -0.52237025, 0.06973349],
            [-0.57908681, -0.54215409,  0.01683106,  0.19941052, -0.54358993, -0.56340091],
            [ 0.36043904, -0.28857117, -0.52420958, -0.90249808,  0.74278826, 0.70469376],
            [ 0.07784766, -0.64352802, -0.84122326, -0.32708089,  0.27241266, 0.86394226],
            [ 0.44508538, -0.63388016,  0.84706968, -0.71098094,  0.5239476 , 0.19557921]
        ])
        p, Aik, Ajk = pk.trajopt.primer_vector(DVi, DVj, Mji, Mjk)
        self.assertTrue(float_rel_error(np.linalg.norm(p), 2.6981957244221193) < 1e-12)
        self.assertTrue(float_rel_error(np.linalg.norm(Aik), 3.7165788775706137) < 1e-12)
        self.assertTrue(float_rel_error(np.linalg.norm(Ajk), 5.1268002110392725) < 1e-12)
class trajopt_zoh_point2point_tests(_ut.TestCase):
    def test_gradient_uniform(self):
        ta = pk.ta.get_zoh_kep(1e-14)
        ta_var = pk.ta.get_zoh_kep_var(1e-10)

        udp = pk.trajopt.zoh_point2point(
            nseg=5,
            tas=(ta, ta_var),
            time_encoding='uniform',
        )

        lb, ub = udp.get_bounds()
        x = np.array([(l + u) / 2.0 for l, u in zip(lb, ub)])

        # Compute analytical gradient (dense, flattened nf×nx)
        ag = udp.gradient(x)
        nf = 1 + udp.get_nec()
        nx = len(x)
        J_analytical = np.array(ag).reshape((nf, nx), order="C")

        # Compute numerical gradient
        J_numerical_flat = pg.estimate_gradient_h(callable=udp.fitness, x=x, dx=1e-6)
        J_numerical = np.array(J_numerical_flat).reshape((nf, nx), order="C")

        self.assertTrue(np.allclose(J_analytical, J_numerical, atol=1e-5, rtol=1e-5))

    def test_gradient_softmax(self):
        ta = pk.ta.get_zoh_kep(1e-14)
        ta_var = pk.ta.get_zoh_kep_var(1e-10)

        udp = pk.trajopt.zoh_point2point(
            nseg=5,
            tas=(ta, ta_var),
            time_encoding='softmax',
        )

        lb, ub = udp.get_bounds()
        x = np.array([(l + u) / 2.0 for l, u in zip(lb, ub)])

        # Compute analytical gradient (dense, flattened nf×nx)
        ag = udp.gradient(x)
        nf = 1 + udp.get_nec()
        nx = len(x)
        J_analytical = np.array(ag).reshape((nf, nx), order="C")

        # Compute numerical gradient
        J_numerical_flat = pg.estimate_gradient_h(callable=udp.fitness, x=x, dx=1e-6)
        J_numerical = np.array(J_numerical_flat).reshape((nf, nx), order="C")

        self.assertTrue(np.allclose(J_analytical, J_numerical, atol=1e-5, rtol=1e-5))
class trajopt_zoh_pl2pl_tests(_ut.TestCase):
    def test_gradient_uniform(self):
        # Build the variational integrators needed for gradients
        ta = pk.ta.get_zoh_kep(1e-14)
        ta_var = pk.ta.get_zoh_kep_var(1e-10)

        # Create a small instance for testing
        udp = pk.trajopt.zoh_pl2pl(
            nseg=5,
            tas=(ta, ta_var),
            time_encoding='uniform',
            t0_bounds=[6700.0, 6800.0],
            tof_bounds=[3.4, 8.6],
            vinf_dep_bounds=[0.0, 0.2],
            vinf_arr_bounds=[0.0, 0.2],
        )

        lb, ub = udp.get_bounds()
        # Use midpoint of bounds as test point; set direction vectors to non-zero
        x = np.array([(l + u) / 2.0 for l, u in zip(lb, ub)])
        # Ensure departure and arrival direction vectors are non-zero (they default to midpoint of [-1, 1] = 0)
        x[3], x[4], x[5] = 1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3)
        x[7], x[8], x[9] = 1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3)

        # Compute analytical gradient (returned as dense, flattened nf×nx matrix)
        ag = udp.gradient(x)
        nf = 1 + udp.get_nec()
        nx = len(x)
        J_analytical = np.array(ag).reshape((nf, nx), order="C")

        # Compute numerical gradient using central differences
        J_numerical_flat = pg.estimate_gradient_h(callable=udp.fitness, x=x, dx=1e-6)
        J_numerical = np.array(J_numerical_flat).reshape((nf, nx), order="C")

        # Compare analytical and numerical gradients
        self.assertTrue(np.allclose(J_analytical, J_numerical, atol=1e-4, rtol=1e-4))
    def test_gradient_softmax(self):
        # Build the variational integrators needed for gradients
        ta = pk.ta.get_zoh_kep(1e-14)
        ta_var = pk.ta.get_zoh_kep_var(1e-10)
        
        #now for the softmax:
        udp_softmax = pk.trajopt.zoh_pl2pl(
            nseg=5,
            tas=(ta, ta_var),
            time_encoding='softmax',
            t0_bounds=[6700.0, 6800.0],
            tof_bounds=[3.4, 8.6],
            vinf_dep_bounds=[0.0, 0.2],
            vinf_arr_bounds=[0.0, 0.2],
        )
        lb, ub = udp_softmax.get_bounds()
        x = np.array([(l + u) / 2.0 for l, u in zip(lb, ub)])
        x[3], x[4], x[5] = 1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3)
        x[7], x[8], x[9] = 1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3)
        nf = 1 + udp_softmax.get_nec()
        nx = len(x)
        ag = udp_softmax.gradient(x)
        J_analytical = np.array(ag).reshape((nf, nx), order="C")
        J_numerical_flat = pg.estimate_gradient_h(callable=udp_softmax.fitness, x=x, dx=1e-6)
        J_numerical = np.array(J_numerical_flat).reshape((nf, nx), order="C")
        self.assertTrue(np.allclose(J_analytical, J_numerical, atol=1e-4, rtol=1e-4))