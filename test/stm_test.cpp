// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/stm.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using xt::linalg::dot;

TEST_CASE("stm_reynolds")
{
    // Test that the STM is the identity if tof = 0
    {
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1., 0., 0.}, {0., 1., 0.}}};
        std::array<std::array<double, 3>, 2> pos_velf = {{{1, 0., 0.}, {0., 1., 0.}}};
        double tof = 0;
        double mu = 1.;
        auto computed = kep3::stm_reynolds(pos_vel0, pos_velf, tof, mu);
        std::array<double, 36> real{1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
        REQUIRE(kep3_tests::L_infinity_norm(computed, real) < 1e-16);
    }
    // Test a case for ellipses (ground truth obtained from heyoka numerical propagation of the variational equations)
    {
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.0, -0.1, 0.1}, {-0.1, 1.0, -0.1}}};
        std::array<std::array<double, 3>, 2> pos_velf
            = {{{-0.8287768689824112, -0.3754546412562854, -0.041211111611466},
                {0.2280234730659919, -1.0912316240587416, 0.119932281556794}}};
        double tof = 2.5;
        double mu = 1.1;
        auto computed = kep3::stm_reynolds(pos_vel0, pos_velf, tof, mu);
        std::array<double, 36> real{
            -4.0239266584906375e+00, -2.4418402840266060e-01, -2.6406086966749842e-01, -5.5228775886036308e-01,
            -4.6974357664764570e+00, 4.1891908684062967e-01,  6.9462648764588364e+00,  1.9945075251451927e+00,
            3.7060766250619542e-01,  2.8792848784718301e+00,  6.3862600330165646e+00,  -3.6090337292505165e-01,
            -9.1773804283672533e-01, -2.8306951066303160e-01, -9.3277020303401248e-01, -2.6987369906375991e-01,
            -1.0496961588294416e+00, -3.9206899674304574e-01, -9.0598278389024092e+00, -1.1027442662181999e+00,
            -7.3428953605302905e-01, -2.2081442425288338e+00, -8.9181443163834704e+00, 7.0811105141122388e-01,
            -1.9238348584522358e+00, -8.4747014330838855e-01, -8.6933035567703432e-02, -1.1859991646656800e+00,
            -8.1877937553924118e-01, -1.3149466189047918e-01, -6.5964493675902569e-01, -1.2288436273700180e-02,
            6.1250732945897668e-02,  5.1887648914252108e-03,  -8.3441694841027836e-01, -1.0028936995504940e+00};
        REQUIRE(kep3_tests::L_infinity_norm(computed, real) < 1e-13);
    }
    // Test a case for hyperbolas (ground truth obtained from heyoka numerical propagation of the variational equations)
    {
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.0, -0.1, 0.1}, {-0.1, 4.0, -0.1}}};
        std::array<std::array<double, 3>, 2> pos_velf
            = {{{0.0620813012976075, 9.381317850166885, -0.20554057619272},
                {-0.4006410093516651, 3.7283262854102937, -0.1232579847955077}}};
        double tof = 2.5;
        double mu = 1.1;
        auto computed = kep3::stm_reynolds(pos_vel0, pos_velf, tof, mu);
        std::array<double, 36> real{
            1.8050918967262499e+00,  6.2308736854929558e-01,  1.3331791887418351e-01, 2.5422350090621637e+00,
            1.5974372140151119e-01,  1.4916539150549199e-02,  7.0165876441824671e-01, 9.2648969848674534e-01,
            5.4392238849208477e-02,  1.6835227906166478e-01,  2.6116515738430328e+00, 1.0615903135503009e-02,
            1.3154563174931996e-01,  4.6712327974799746e-02,  3.0916182382747365e-01, 1.4722361158214910e-02,
            9.7744651687210766e-03,  2.3539819547075598e+00,  3.2308728658335789e-01, 2.7652591046357017e-01,
            5.5467674302453585e-02,  1.0074022094092847e+00,  7.5520868950613509e-02, 6.4118178953852549e-03,
            3.2054236344814069e-01,  -1.3743648698706183e-03, 2.4410253666026935e-02, 8.1021085851723754e-02,
            1.0723653741971355e+00,  4.5812590638524440e-03,  5.4474821979493340e-02, 2.0107893599865902e-02,
            -3.0333231546807449e-01, 6.2877528525030705e-03,  4.0436438780296381e-03, 9.2489986970891724e-01};
        REQUIRE(kep3_tests::L_infinity_norm(computed, real) < 1e-13);
    }
}

TEST_CASE("propagate_stm_reynolds")
{
    { // We test the identity stm02 = stm12stm01 (ellipses)
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.23, -0.12, 0.12}, {-0.12, 1.23, 0.12}}};
        double mu = 1.02;
        double dt1 = 1.1;
        double dt2 = 2.3;
        auto res01 = kep3::propagate_stm_reynolds(pos_vel0, dt1, mu, true);
        auto res02 = kep3::propagate_stm_reynolds(pos_vel0, dt2, mu, true);
        // ... and compute the stm from 1 to 2. (pos_vel1 wil also be propagated to pos_vel2)
        auto res12 = kep3::propagate_stm_reynolds(res01.first, dt2 - dt1, mu, true);
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M01 = xt::adapt(res01.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M02 = xt::adapt(res02.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M12 = xt::adapt(res12.second.value(), {6, 6});
        REQUIRE(xt::linalg::norm(M02 - dot(M12, M01)) < 1e-13);
    }
    { // We test the identity stm02 = stm12stm01 (hyperbolas)
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.23, -0.12, 0.12}, {-4.12, 1.23, 0.12}}};
        double mu = 1.02;
        double dt1 = 1.6;
        double dt2 = 3.4;
        auto res01 = kep3::propagate_stm_reynolds(pos_vel0, dt1, mu, true);
        auto res02 = kep3::propagate_stm_reynolds(pos_vel0, dt2, mu, true);
        // ... and compute the stm from 1 to 2. (pos_vel1 wil also be propagated to pos_vel2)
        auto res12 = kep3::propagate_stm_reynolds(res01.first, dt2 - dt1, mu, true);
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M01 = xt::adapt(res01.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M02 = xt::adapt(res02.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M12 = xt::adapt(res12.second.value(), {6, 6});
        REQUIRE(xt::linalg::norm(M02 - dot(M12, M01)) < 1e-8);
    }
}

TEST_CASE("propagate_lagrangian(stm)")
{
    { // We test the identity stm02 = stm12stm01 (ellipses)
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.23, -0.12, 0.12}, {-0.12, 1.23, 0.12}}};
        double mu = 1.02;
        double dt1 = 1.1;
        double dt2 = 2.3;
        auto res01 = kep3::propagate_lagrangian(pos_vel0, dt1, mu, true);
        auto res02 = kep3::propagate_lagrangian(pos_vel0, dt2, mu, true);
        // ... and compute the stm from 1 to 2. (pos_vel1 wil also be propagated to pos_vel2)
        auto res12 = kep3::propagate_lagrangian(res01.first, dt2 - dt1, mu, true);
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M01 = xt::adapt(res01.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M02 = xt::adapt(res02.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M12 = xt::adapt(res12.second.value(), {6, 6});
        REQUIRE(xt::linalg::norm(M02 - dot(M12, M01)) < 1e-13);
    }
    { // We test the identity stm02 = stm12stm01 (hyperbolas)
        std::array<std::array<double, 3>, 2> pos_vel0 = {{{1.23, -0.12, 0.12}, {-4.12, 1.23, 0.12}}};
        double mu = 1.02;
        double dt1 = 1.6;
        double dt2 = 3.4;
        auto res01 = kep3::propagate_lagrangian(pos_vel0, dt1, mu, true);
        auto res02 = kep3::propagate_lagrangian(pos_vel0, dt2, mu, true);
        // ... and compute the stm from 1 to 2. (pos_vel1 wil also be propagated to pos_vel2)
        auto res12 = kep3::propagate_lagrangian(res01.first, dt2 - dt1, mu, true);
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M01 = xt::adapt(res01.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M02 = xt::adapt(res02.second.value(), {6, 6});
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        auto M12 = xt::adapt(res12.second.value(), {6, 6});
        REQUIRE(xt::linalg::norm(M02 - dot(M12, M01)) < 1e-12);
    }
}

TEST_CASE("reynolds_vs_lagrange")
{
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(1220202343u);
    { // Targeting Ellipses
        std::uniform_real_distribution<double> sma_d(1.1, 10.);
        std::uniform_real_distribution<double> ecc_d(0, 0.9);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(-2. * kep3::pi, 2. * kep3::pi);
        std::uniform_real_distribution<double> mu_d(0.8, 10);

        auto N = 1000u;
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);
            double tof = time_d(rng_engine);
            double mu = mu_d(rng_engine);

            std::array<double, 6> par = {sma, ecc, incl, Omega, omega, f};
            auto pos_vel = kep3::par2ic(par, 1.);
            auto pos_velr = pos_vel;
            auto res = kep3::propagate_lagrangian(pos_vel, tof, mu, true);
            auto res_r = kep3::propagate_stm_reynolds(pos_velr, tof, mu, true);
            REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[0], pos_velr[0]) < 1e-14);
            REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[1], pos_velr[1]) < 1e-14);
            // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
            REQUIRE(kep3_tests::L_infinity_norm(res.second.value(), res_r.second.value()) < 1e-10);
        }
    }

    { // Targeting Hyperbolas
        std::uniform_real_distribution<double> sma_d(-1.1, -10.);
        std::uniform_real_distribution<double> ecc_d(1.1, 10.);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(-2. * kep3::pi, 2. * kep3::pi);
        std::uniform_real_distribution<double> mu_d(0.8, 10);

        auto N = 1000u;
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);
            double tof = time_d(rng_engine);
            double mu = mu_d(rng_engine);
            if (std::cos(f) > 0.) {
                std::array<double, 6> par = {sma, ecc, incl, Omega, omega, f};
                auto pos_vel = kep3::par2ic(par, 1.);
                auto pos_velr = pos_vel;
                auto res = kep3::propagate_lagrangian(pos_vel, tof, mu, true);
                auto res_r = kep3::propagate_stm_reynolds(pos_velr, tof, mu, true);
                REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[0], pos_velr[0]) < 1e-14);
                REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[1], pos_velr[1]) < 1e-14);
                // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
                REQUIRE(kep3_tests::L_infinity_norm(res.second.value(), res_r.second.value()) < 1e-8);
            }
        }
    }
}
