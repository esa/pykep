// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cmath>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using Catch::Matchers::WithinRel;
using kep3::mee2ic;
using kep3::ic2mee;
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
            auto eq = ic2mee(pos_vel, 1.);
            auto pos_vel_new = mee2ic(eq, 1.0);
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
            auto eq = ic2mee(pos_vel, 1.);
            auto pos_vel_new = mee2ic(eq, 1.0);

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
            auto eq = ic2mee(pos_vel, 1., true);
            auto pos_vel_new = mee2ic(eq, 1.0, true);
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
            auto eq = ic2mee(pos_vel, 1., true);
            auto pos_vel_new = mee2ic(eq, 1.0, true);

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
