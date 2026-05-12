// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <heyoka/expression.hpp>
#include <heyoka/kw.hpp>
#include <heyoka/math/cos.hpp>
#include <heyoka/math/sin.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using kep3::ic2mee;
using kep3::mee2ic;
using kep3::pi;

TEST_CASE("fb_con")
{
    // Zero inclination and eccentricity
    {
        auto par = ic2mee({{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}}, 1.0);
        REQUIRE(par[0] == 1.); // p is 1
        REQUIRE(par[1] == 0.); // f is zero
        REQUIRE(par[2] == 0.); // g is zero
        REQUIRE(par[0] == 1.); // h is 1
        REQUIRE(par[1] == 0.); // k is zero
        REQUIRE(par[2] == 0.); // L is zero
    }
    // Orbit at 90 degrees inclination
    {
        auto par = ic2mee({{{1.0, 0.0, 0.0}, {0.0, 0.0, 1.1}}}, 1.0);
        REQUIRE_THAT(par[0], WithinRel(1.2658227848101269 * (1. - 0.21 * 0.21), 1e-14));
        REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
        REQUIRE(par[2] == 0.); // f is zero
        REQUIRE(par[3] == 1.); // h is 1
        REQUIRE(par[4] == 0.); // k is zero
        REQUIRE(par[5] == 0.); // L is zero
    }
    // Orbit at 90 degrees inclination
    {
        auto par = ic2mee({{{1.0, 0.0, 0.0}, {0.0, 0.0, -1.1}}}, 1.0);
        REQUIRE_THAT(par[0], WithinRel(1.2658227848101269 * (1. - 0.21 * 0.21), 1e-14));
        REQUIRE_THAT(par[1], WithinRel(0.21, 1e-14));
        REQUIRE(par[2] == 0.);  // f is zero
        REQUIRE(par[3] == -1.); // h is 1
        REQUIRE(par[4] == 0.);  // k is zero
        REQUIRE(par[5] == 0.);  // L is zero
    }
}

TEST_CASE("ic2mee2ic")
{
    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    // Distributions for the elements
    std::uniform_real_distribution<double> sma_d(1.1, 100.);
    std::uniform_real_distribution<double> ecc_d(0, 0.99);
    std::uniform_real_distribution<double> incl_d(0., pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
    std::uniform_real_distribution<double> omega_d(0., pi);
    std::uniform_real_distribution<double> ni_d(0, 2 * pi);

    {
        // Testing on N random calls on ellipses
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Compute the initial r,v
            auto pos_vel = kep3::par2ic({sma, ecc, incl, Omega, omega, ni}, 1.);
            // Test ic2mee2ic
            auto mee = ic2mee(pos_vel, 1.);
            auto pos_vel_new = mee2ic(mee, 1.0);
            double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1]
                                 + pos_vel[0][2] * pos_vel[0][2]);
            double V = std::sqrt(pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1]
                                 + pos_vel[1][2] * pos_vel[1][2]);
            double R_new = std::sqrt(pos_vel_new[0][0] * pos_vel_new[0][0] + pos_vel_new[0][1] * pos_vel_new[0][1]
                                     + pos_vel_new[0][2] * pos_vel_new[0][2]);
            double V_new = std::sqrt(pos_vel_new[1][0] * pos_vel_new[1][0] + pos_vel_new[1][1] * pos_vel_new[1][1]
                                     + pos_vel_new[1][2] * pos_vel_new[1][2]);
            // Here we do not use catch matchers to test floating point as for small numbers (<1) we care about absolute
            // while for large (>1) we care for relative error.
            double rel_err_V = std::abs(V_new - V) / std::max(1., std::max(V_new, V));
            double rel_err_R = std::abs(R_new - R) / std::max(1., std::max(R_new, R));
            REQUIRE(rel_err_V < 1e-13);
            REQUIRE(rel_err_R < 1e-13);
        }
    }
    {
        // Testing on N random calls on ellipses
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = -sma_d(rng_engine);
            double ecc = ecc_d(rng_engine) + 1.1;
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Skipping if true anomaly is way out of asymptotes
            if (std::cos(ni) < -1 / ecc + 0.1) {
                continue;
            }
            // Compute the initial r,v
            auto pos_vel = kep3::par2ic({sma, ecc, incl, Omega, omega, ni}, 1.);
            // Test ic2mee2ic
            auto mee = ic2mee(pos_vel, 1.);
            auto pos_vel_new = mee2ic(mee, 1.0);

            double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1]
                                 + pos_vel[0][2] * pos_vel[0][2]);
            double V = std::sqrt(pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1]
                                 + pos_vel[1][2] * pos_vel[1][2]);
            double R_new = std::sqrt(pos_vel_new[0][0] * pos_vel_new[0][0] + pos_vel_new[0][1] * pos_vel_new[0][1]
                                     + pos_vel_new[0][2] * pos_vel_new[0][2]);
            double V_new = std::sqrt(pos_vel_new[1][0] * pos_vel_new[1][0] + pos_vel_new[1][1] * pos_vel_new[1][1]
                                     + pos_vel_new[1][2] * pos_vel_new[1][2]);
            // Here we do not use catch matchers to test floating point as for small numbers (<1) we care about absolute
            // while for large (>1) we care for relative error.
            double rel_err_V = std::abs(V_new - V) / std::max(1., std::max(V_new, V));
            double rel_err_R = std::abs(R_new - R) / std::max(1., std::max(R_new, R));
            REQUIRE(rel_err_V < 1e-13);
            REQUIRE(rel_err_R < 1e-13);
        }
    }
}

TEST_CASE("ic2mee2ic_retrogade")
{
    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    // Distributions for the elements
    std::uniform_real_distribution<double> sma_d(1.1, 100.);
    std::uniform_real_distribution<double> ecc_d(0, 0.99);
    std::uniform_real_distribution<double> incl_d(0., pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
    std::uniform_real_distribution<double> omega_d(0., pi);
    std::uniform_real_distribution<double> ni_d(0, 2 * pi);

    {
        // Testing on N random calls on ellipses
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Compute the initial r,v
            auto pos_vel = kep3::par2ic({sma, ecc, incl, Omega, omega, ni}, 1.);
            // Test ic2mee2ic
            auto mee = ic2mee(pos_vel, 1., true);
            auto pos_vel_new = mee2ic(mee, 1.0, true);
            double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1]
                                 + pos_vel[0][2] * pos_vel[0][2]);
            double V = std::sqrt(pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1]
                                 + pos_vel[1][2] * pos_vel[1][2]);
            double R_new = std::sqrt(pos_vel_new[0][0] * pos_vel_new[0][0] + pos_vel_new[0][1] * pos_vel_new[0][1]
                                     + pos_vel_new[0][2] * pos_vel_new[0][2]);
            double V_new = std::sqrt(pos_vel_new[1][0] * pos_vel_new[1][0] + pos_vel_new[1][1] * pos_vel_new[1][1]
                                     + pos_vel_new[1][2] * pos_vel_new[1][2]);
            // Here we do not use catch matchers to test floating point as for small numbers (<1) we care about absolute
            // while for large (>1) we care for relative error.
            REQUIRE(kep3_tests::floating_point_error(R, R_new) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(V, V_new) < 1e-13);
        }
    }
    {
        // Testing on N random calls on ellipses
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = -sma_d(rng_engine);
            double ecc = ecc_d(rng_engine) + 1.1;
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Skipping if true anomaly is way out of asymptotes
            if (std::cos(ni) < -1 / ecc + 0.1) {
                continue;
            }
            // Compute the initial r,v
            auto pos_vel = kep3::par2ic({sma, ecc, incl, Omega, omega, ni}, 1.);
            // Test ic2mee2ic
            auto mee = ic2mee(pos_vel, 1., true);
            auto pos_vel_new = mee2ic(mee, 1.0, true);

            double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1]
                                 + pos_vel[0][2] * pos_vel[0][2]);
            double V = std::sqrt(pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1]
                                 + pos_vel[1][2] * pos_vel[1][2]);
            double R_new = std::sqrt(pos_vel_new[0][0] * pos_vel_new[0][0] + pos_vel_new[0][1] * pos_vel_new[0][1]
                                     + pos_vel_new[0][2] * pos_vel_new[0][2]);
            double V_new = std::sqrt(pos_vel_new[1][0] * pos_vel_new[1][0] + pos_vel_new[1][1] * pos_vel_new[1][1]
                                     + pos_vel_new[1][2] * pos_vel_new[1][2]);
            // Here we do not use catch matchers to test floating point as for small numbers (<1) we care about absolute
            // while for large (>1) we care for relative error.
            REQUIRE(kep3_tests::floating_point_error(R, R_new) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(V, V_new) < 1e-13);
        }
    }
}

TEST_CASE("ic2mee2ic_symbolic")
{
    using namespace heyoka;

    // Get symbolic expressions and Jacobians
    auto [ic2mee_exprs, ic2mee_jac_opt] = kep3::ic2mee(true);
    auto [mee2ic_exprs, mee2ic_jac_opt] = kep3::mee2ic(true);

    REQUIRE(ic2mee_jac_opt.has_value());
    REQUIRE(mee2ic_jac_opt.has_value());

    // Symbolic variables for ic2mee transformation
    auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");
    // Symbolic variables for mee2ic transformation
    auto [p, f, g, h, k, L] = make_vars("p", "f", "g", "h", "k", "L");

    // Compile Jacobians to callable functions
    auto J_ic2mee = cfunc<double>(ic2mee_jac_opt.value(), {x, y, z, vx, vy, vz});
    auto J_mee2ic = cfunc<double>(mee2ic_jac_opt.value(), {p, f, g, h, k, L});

    // Compile forward/inverse transformations
    auto ic2mee_cf = cfunc<double>(ic2mee_exprs, {x, y, z, vx, vy, vz});
    auto mee2ic_cf = cfunc<double>(mee2ic_exprs, {p, f, g, h, k, L});

    // Helper: reshape flat Jacobian (36 elements) to 6x6 matrix
    auto reshape_jac = [](const std::vector<double> &jac) {
        std::array<std::array<double, 6>, 6> mat;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                mat[i][j] = jac[i * 6 + j];
            }
        }
        return mat;
    };

    // Helper: multiply two 6x6 matrices
    auto matmul = [](const std::array<std::array<double, 6>, 6> &A, const std::array<std::array<double, 6>, 6> &B) {
        std::array<std::array<double, 6>, 6> C = {};
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                for (std::size_t m = 0; m < 6; ++m) {
                    C[i][j] += A[i][m] * B[m][j];
                }
            }
        }
        return C;
    };

    // Random test engine
    std::mt19937 rng_engine(122012203u);
    std::uniform_real_distribution<double> sma_d(1.1, 10.);
    std::uniform_real_distribution<double> ecc_d(0, 0.5);
    std::uniform_real_distribution<double> incl_d(0., pi/2);
    std::uniform_real_distribution<double> angle_d(0, 2 * pi);

    for (auto t = 0u; t < 100u; ++t) {
        // Generate random Keplerian state
        auto sma = sma_d(rng_engine);
        auto ecc = ecc_d(rng_engine);
        auto incl = incl_d(rng_engine);
        auto Omega = angle_d(rng_engine);
        auto omega = angle_d(rng_engine);
        auto nu = angle_d(rng_engine);

        auto pos_vel = kep3::par2ic({sma, ecc, incl, Omega, omega, nu}, 1.);
        std::vector<double> ic_vec{pos_vel[0][0], pos_vel[0][1], pos_vel[0][2],
                                   pos_vel[1][0], pos_vel[1][1], pos_vel[1][2]};

        // TEST 1: Heyoka compiled function matches C++ version
        auto mee_cpp = kep3::ic2mee(pos_vel, 1.);
        std::vector<double> mee_heyoka(6), J1_flat(36), J2_flat(36);
        ic2mee_cf(mee_heyoka, ic_vec, kw::pars = {1., 1.});

        for (auto i = 0u; i < 5u; ++i) { // Check p, f, g, h, k (skip L for angle wrapping)
            REQUIRE_THAT(mee_heyoka[i], WithinRel(mee_cpp[i], 1e-10));
        }

        // TEST 2: J_ic2mee @ J_mee2ic = Identity
        J_ic2mee(J1_flat, ic_vec, kw::pars = {1., 1.});
        J_mee2ic(J2_flat, mee_heyoka, kw::pars = {1., 1.});

        auto J1 = reshape_jac(J1_flat);
        auto J2 = reshape_jac(J2_flat);
        auto product = matmul(J1, J2);

        // Check that product is identity (with tolerance for numerical precision)
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                REQUIRE(kep3_tests::floating_point_error(product[i][j], expected) < 1e-13);
            }
        }

        // TEST 3: Inverse transformation consistency
        auto pos_vel_roundtrip
            = mee2ic({mee_heyoka[0], mee_heyoka[1], mee_heyoka[2], mee_heyoka[3], mee_heyoka[4], mee_heyoka[5]}, 1.0);
        double R
            = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1] + pos_vel[0][2] * pos_vel[0][2]);
        double V
            = std::sqrt(pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1] + pos_vel[1][2] * pos_vel[1][2]);
        double R_rt = std::sqrt(pos_vel_roundtrip[0][0] * pos_vel_roundtrip[0][0]
                                + pos_vel_roundtrip[0][1] * pos_vel_roundtrip[0][1]
                                + pos_vel_roundtrip[0][2] * pos_vel_roundtrip[0][2]);
        double V_rt = std::sqrt(pos_vel_roundtrip[1][0] * pos_vel_roundtrip[1][0]
                                + pos_vel_roundtrip[1][1] * pos_vel_roundtrip[1][1]
                                + pos_vel_roundtrip[1][2] * pos_vel_roundtrip[1][2]);

        REQUIRE_THAT(R_rt, WithinRel(R, 1e-13));
        REQUIRE_THAT(V_rt, WithinRel(V, 1e-13));
    }
}
