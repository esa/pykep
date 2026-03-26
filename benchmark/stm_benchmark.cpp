// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <chrono>
#include <random>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/stm.hpp>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

// In this benchmark we test the speed and accuracy of the Lagrangian
// propagation solvers
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void perform_test_speed(double min_ecc, double max_ecc, unsigned N, bool warming = false)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> sma_d(0.5, 20.);
    std::uniform_real_distribution<double> ecc_d(min_ecc, max_ecc);
    std::uniform_real_distribution<double> incl_d(0., kep3::pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
    std::uniform_real_distribution<double> omega_d(0., 2 * kep3::pi);
    std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
    std::uniform_real_distribution<double> tof_d(10., 100.);

    // We generate the random dataset
    std::vector<std::array<std::array<double, 3>, 2>> pos_vels(N);
    std::vector<double> tofs(N);
    for (auto i = 0u; i < N; ++i) {
        auto ecc = ecc_d(rng_engine);
        auto sma = sma_d(rng_engine);
        ecc > 1. ? sma = -sma : sma;
        double f = kep3::pi;
        while (std::cos(f) < -1. / ecc && sma < 0.) {
            f = f_d(rng_engine);
        }
        pos_vels[i] = kep3::par2ic({sma, ecc, incl_d(rng_engine), Omega_d(rng_engine), omega_d(rng_engine), f}, 1.);
        tofs[i] = tof_d(rng_engine);
        // fmt::print("[{}, {},{},{},{},{},{}],", pos_vels[i][0][0], pos_vels[i][0][1], pos_vels[i][0][2],
        // pos_vels[i][1][0],
        //            pos_vels[i][1][1], pos_vels[i][1][2], tofs[i]);
    }

    if (warming) {
        for (auto i = 0u; i < N; ++i) {
            kep3::propagate_lagrangian(pos_vels[i], tofs[i], 1., true);
        }
    } else {
        // We log progress
        fmt::print("{:.2f} min_ecc, {:.2f} max_ecc, on {} data points: ", min_ecc, max_ecc, N);

        auto start = high_resolution_clock::now();
        for (auto i = 0u; i < N; ++i) {
            kep3::propagate_lagrangian(pos_vels[i], tofs[i], 1., true);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        fmt::print("{:.3f}s\n", (static_cast<double>(duration.count()) / 1e6));
    }
}

int main()
{
    // warming up
    perform_test_speed(0, 0.5, 100000, true);
    // performing tests
    fmt::print("\nComputes speed at different eccentricity ranges:\n");
    perform_test_speed(0, 0.5, 100000);
    perform_test_speed(0.5, 0.9, 100000);
    perform_test_speed(1.1, 2., 100000);
}