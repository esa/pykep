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

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/mee2par2mee.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

constexpr double pi{boost::math::constants::pi<double>()};

TEST_CASE("mee2par2mee")
{
    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    // Distributions for the elements
    std::uniform_real_distribution<double> sma_d(1.1, 100.);
    std::uniform_real_distribution<double> ecc_d(0, 0.9);
    std::uniform_real_distribution<double> incl_d(0., pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
    std::uniform_real_distribution<double> omega_d(0., 2 * pi);
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
            // Compute the initial eq
            auto eq = kep3::par2mee({sma, ecc, incl, Omega, omega, ni});
            // Test mee2par2mee
            auto par = kep3::mee2par(eq);

            // Here we do not use catch matchers to test floating point as for small
            // numbers (<1) we care about absolute while for large (>1) we care for
            // relative error.
            REQUIRE(kep3_tests::floating_point_error(sma, par[0]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ecc, par[1]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(incl, par[2]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(Omega, par[3]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(omega, par[4]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ni, par[5]) < 1e-13);
        }
    }
    {
        // Testing on N random calls on hyperbolas
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = -sma_d(rng_engine);
            double ecc = ecc_d(rng_engine) + 2.;
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Skipping if true anomaly is way out of asymptotes
            if (std::cos(ni) < -1 / ecc + 0.1) {
                continue;
            }
            // Compute the initial eq
            auto eq = kep3::par2mee({sma, ecc, incl, Omega, omega, ni});
            // Test mee2par2mee
            auto par = kep3::mee2par(eq);

            // Here we do not use catch matchers to test floating point as for small
            // numbers (<1) we care about absolute while for large (>1) we care for
            // relative error.
            REQUIRE(kep3_tests::floating_point_error(sma, par[0])
                    < 1e-13); // errors arise when e is close to 1 and a is high
            REQUIRE(kep3_tests::floating_point_error(ecc, par[1]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(incl, par[2]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(Omega, par[3]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(omega, par[4]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ni, par[5]) < 1e-13);
        }
    }
}

TEST_CASE("mee2par2mee_retrogade")
{
    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    // Distributions for the elements
    std::uniform_real_distribution<double> sma_d(1.1, 100.);
    std::uniform_real_distribution<double> ecc_d(0, 0.99);
    std::uniform_real_distribution<double> incl_d(0., pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * pi);
    std::uniform_real_distribution<double> omega_d(0., 2 * pi);
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
            // Compute the initial eq
            auto eq = kep3::par2mee({sma, ecc, incl, Omega, omega, ni}, true);
            // Test mee2par2mee
            auto par = kep3::mee2par(eq, true);

            // Here we do not use catch matchers to test floating point as for small
            // numbers (<1) we care about absolute while for large (>1) we care for
            // relative error.
            REQUIRE(kep3_tests::floating_point_error(sma, par[0]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ecc, par[1]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(incl, par[2]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(Omega, par[3]) < 1e-13);
            if (kep3_tests::floating_point_error(omega, par[4]) > 1e-13) {
                fmt::print("\n{}, {}", omega / pi * 180, par[4] / pi * 180);
            }
            REQUIRE(kep3_tests::floating_point_error(omega, par[4]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ni, par[5]) < 1e-13);
        }
    }
    {
        // Testing on N random calls on hyperbolas
        unsigned N = 10000;
        for (auto i = 0u; i < N; ++i) {
            // We sample randomly on the Keplerian space
            double sma = -sma_d(rng_engine);
            double ecc = ecc_d(rng_engine) + 1.;
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double ni = ni_d(rng_engine);
            // Skipping if true anomaly is way out of asymptotes
            if (std::cos(ni) < -1 / ecc + 0.1) {
                continue;
            }
            // Compute the initial eq
            auto eq = kep3::par2mee({sma, ecc, incl, Omega, omega, ni}, true);
            // Test mee2par2mee
            auto par = kep3::mee2par(eq, true);

            // Here we do not use catch matchers to test floating point as for small
            // numbers (<1) we care about absolute while for large (>1) we care for
            // relative error.
            REQUIRE(kep3_tests::floating_point_error(sma, par[0])
                    < 1e-10); // errors arise since p = a * (1-e^2) and a = p / (1-e^2)
                              // [when e is close to 1 and a is high]
            REQUIRE(kep3_tests::floating_point_error(ecc, par[1]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(incl, par[2]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(Omega, par[3]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(omega, par[4]) < 1e-13);
            REQUIRE(kep3_tests::floating_point_error(ni, par[5]) < 1e-13);
        }
    }
}