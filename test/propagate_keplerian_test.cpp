// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <functional>
#include <stdexcept>

#include <fmt/base.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using kep3::propagate_keplerian;

void test_propagate_keplerian(
    const std::function<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>(
        const std::array<std::array<double, 3>, 2> &, double, double, bool)> &propagate,
    unsigned int N = 1000000)
{

    // Engines
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(12523334u);



    { // Targeting Ellipses
        std::uniform_real_distribution<double> sma_d(1.1, 10.);
        std::uniform_real_distribution<double> ecc_d(0, 0.9);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., 2 * kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(-2. * kep3::pi, 2. * kep3::pi);

        // Testing on N random calls on ellipses
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);

            std::array<double, 6> par = {sma, ecc, incl, Omega, omega, f};
            auto pos_vel0 = kep3::par2ic(par, 1.1);
            double tof = time_d(rng_engine);
            auto res0 = propagate(pos_vel0, tof, 1.1, false);
            auto res1 = propagate(res0.first, -tof, 1.1, false);
            // precision loss around incl = pi and generic .. not sure why. Worstened after change of par2ic to differentiable version
            // using atan2 and arrays (corelated?)
            REQUIRE(kep3_tests::floating_point_error_vector(res1.first[0], pos_vel0[0]) < 1e-10);
            REQUIRE(kep3_tests::floating_point_error_vector(res1.first[1], pos_vel0[1]) < 1e-10);
        }
    }

    { // Targeting Hyperbolas
        std::uniform_real_distribution<double> sma_d(-100, -1.1);
        std::uniform_real_distribution<double> ecc_d(2., 20.);
        std::uniform_real_distribution<double> incl_d(0., kep3::pi);
        std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> omega_d(0., kep3::pi);
        std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
        std::uniform_real_distribution<double> time_d(0.1, 20.);
        // Testing on N random calls on hyperbolas
        for (auto i = 0u; i < N; ++i) {
            double sma = sma_d(rng_engine);
            double ecc = ecc_d(rng_engine);
            double incl = incl_d(rng_engine);
            double Omega = Omega_d(rng_engine);
            double omega = omega_d(rng_engine);
            double f = f_d(rng_engine);
            if (std::cos(f) > -1 / ecc) {
                std::array<double, 6> par = {sma, ecc, incl, Omega, omega, f};
                auto pos_vel0 = kep3::par2ic(par, 1.);
                double tof = time_d(rng_engine);
                auto res0 = propagate(pos_vel0, tof, 1., false);
                auto res1 = propagate(res0.first, -tof, 1., false);
                REQUIRE(kep3_tests::floating_point_error_vector(res1.first[0], pos_vel0[0]) < 1e-10);
                REQUIRE(kep3_tests::floating_point_error_vector(res1.first[1], pos_vel0[1]) < 1e-10);
            }
        }
    }
}

TEST_CASE("propagate_keplerian")
{
    // We test both Normal and Universal variables version nwith the same data.
    test_propagate_keplerian(&propagate_keplerian, 10000u);
}