// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <chrono>
#include <functional>

#include <random>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

// In this benchmark we test the speed and accuracy of the Lagrangian
// propagation solvers

void perform_test_speed(
    double min_ecc, double max_ecc, unsigned N,
    const std::function<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>(
        const std::array<std::array<double, 3>, 2> &, double, double, bool)> &propagate)
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
    }

    // We log progress
    fmt::print("{:.2f} min_ecc, {:.2f} max_ecc, on {} data points: ", min_ecc, max_ecc, N);

    auto start = high_resolution_clock::now();
    for (auto i = 0u; i < N; ++i) {
        propagate(pos_vels[i], tofs[i], 1., false);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("{:.3f}s\n", (static_cast<double>(duration.count()) / 1e6));
}

void perform_test_accuracy(
    double min_ecc, double max_ecc, unsigned N,
    const std::function<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>(
        const std::array<std::array<double, 3>, 2> &, double, double, bool)> &propagate)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(53534535u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> sma_d(0.5, 10.);
    std::uniform_real_distribution<double> ecc_d(min_ecc, max_ecc);
    std::uniform_real_distribution<double> incl_d(0., kep3::pi);
    std::uniform_real_distribution<double> Omega_d(0, 2 * kep3::pi);
    std::uniform_real_distribution<double> omega_d(0., 2 * kep3::pi);
    std::uniform_real_distribution<double> f_d(0, 2 * kep3::pi);
    std::uniform_real_distribution<double> tof_d(10, 100.);

    // We generate the random dataset
    std::vector<std::array<std::array<double, 3>, 2>> pos_vels(N);
    std::vector<double> tofs(N);
    for (auto i = 0u; i < N; ++i) {
        double f = kep3::pi;
        double ecc = 10.;
        double sma = -1.;
        while (std::cos(f) < -1. / ecc && sma < 0.) {
            ecc = ecc_d(rng_engine);
            sma = sma_d(rng_engine);
            ecc > 1. ? sma = -sma : sma;
            f = f_d(rng_engine);
        }

        pos_vels[i] = kep3::par2ic({sma, ecc, incl_d(rng_engine), Omega_d(rng_engine), omega_d(rng_engine), f}, 1.);
        tofs[i] = tof_d(rng_engine);
    }
    // We log progress
    fmt::print("{:.2f} min_ecc, {:.2f} max_ecc, on {} data points: ", min_ecc, max_ecc, N);
    std::vector<double> err(N);
    for (auto i = 0u; i < N; ++i) {
        auto res1 = propagate(pos_vels[i], tofs[i], 1., false);
        auto res2 = propagate(res1.first, -tofs[i], 1., false);
        err[i] = std::abs(pos_vels[i][0][0] - res2.first[0][0])
                 + std::abs(pos_vels[i][0][1] - res2.first[0][1])
                 + std::abs(pos_vels[i][0][2] - res2.first[0][2])
                 + std::abs(pos_vels[i][1][0] - res2.first[1][0])
                 + std::abs(pos_vels[i][1][1] - res2.first[1][1])
                 + std::abs(pos_vels[i][1][2] - res2.first[1][2]);
    }
    auto max_it = max_element(std::begin(err), std::end(err));
    auto min_it = min_element(std::begin(err), std::end(err));
    auto avg = std::accumulate(err.begin(), err.end(), 0.0) / static_cast<double>(err.size()) / 6.;
    fmt::print("{:.3e} avg, {:.3e} min, {:.3e} max\n", avg, *min_it, *max_it);
}

int main()
{
    fmt::print("\nComputes speed at different eccentricity ranges:\n");
    perform_test_speed(0, 0.5, 1000000, &kep3::propagate_lagrangian);
    perform_test_speed(0.5, 0.9, 1000000, &kep3::propagate_lagrangian);
    perform_test_speed(0.9, 0.99, 1000000, &kep3::propagate_lagrangian);
    perform_test_speed(1.1, 10., 1000000, &kep3::propagate_lagrangian);

    fmt::print("\nComputes error at different eccentricity ranges:\n");
    perform_test_accuracy(0, 0.5, 100000, &kep3::propagate_lagrangian);
    perform_test_accuracy(0.5, 0.9, 100000, &kep3::propagate_lagrangian);
    perform_test_accuracy(0.9, 0.99, 100000, &kep3::propagate_lagrangian);
    //
    // fmt::print("\nComputes speed at different eccentricity ranges [Universal "
    //           "Anomaly]:\n");
    // perform_test_speed(0, 0.5, 1000000, &kep3::propagate_lagrangian_u);
    // perform_test_speed(0.5, 0.9, 1000000, &kep3::propagate_lagrangian_u);
    // perform_test_speed(0.9, 0.99, 1000000, &kep3::propagate_lagrangian_u);
    // perform_test_speed(1.1, 10., 1000000, &kep3::propagate_lagrangian_u);
    //
    // fmt::print("\nComputes error at different eccentricity ranges [Universal "
    //           "Anomaly]:\n");
    // perform_test_accuracy(0, 0.5, 100000, &kep3::propagate_lagrangian_u);
    // perform_test_accuracy(0.5, 0.9, 100000, &kep3::propagate_lagrangian_u);
    // perform_test_accuracy(0.9, 0.99, 100000, &kep3::propagate_lagrangian_u);
    //
    // fmt::print("\nComputes speed at different eccentricity ranges [keplerian "
    //           "propagation]:\n");
    // perform_test_speed(0, 0.5, 1000000, &kep3::propagate_keplerian);
    // perform_test_speed(0.5, 0.9, 1000000, &kep3::propagate_keplerian);
    // perform_test_speed(0.9, 0.99, 1000000, &kep3::propagate_keplerian);
    //
    // fmt::print("\nComputes error at different eccentricity ranges [keplerian "
    //           "propagation]:\n");
    // perform_test_accuracy(0, 0.5, 100000, &kep3::propagate_keplerian);
    // perform_test_accuracy(0.5, 0.9, 100000, &kep3::propagate_keplerian);
    // perform_test_accuracy(0.9, 0.99, 100000, &kep3::propagate_keplerian);
}