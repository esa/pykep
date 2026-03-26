// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <chrono>
#include <iostream>
#include <random>

#include <fmt/core.h>

#include <kep3/core_astro/convert_anomalies.hpp>

using kep3::e2m;
using kep3::m2e;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

// In this benchmark we test the speed and accuracy of the Kepler's equation
// solvers

void perform_test_speed(double min_ecc, double max_ecc, unsigned N)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> ecc_d(min_ecc, max_ecc);
    std::uniform_real_distribution<double> M_d(-1e8, 1e8);

    // We generate the random dataset
    std::vector<double> eccenricities(N);
    std::vector<double> mean_anomalies(N);

    for (auto i = 0u; i < N; ++i) {
        mean_anomalies[i] = M_d(rng_engine);
        eccenricities[i] = ecc_d(rng_engine);
    }

    // We log progress
    fmt::print("{:.2f} min_ecc, {:.2f} max_ecc, on {} data points: ", min_ecc, max_ecc, N);

    auto start = high_resolution_clock::now();
    for (auto i = 0u; i < N; ++i) {
        m2e(mean_anomalies[i], eccenricities[i]);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("{:.3f}s\n", (static_cast<double>(duration.count()) / 1e6));
}

void perform_test_accuracy(double min_ecc, double max_ecc, unsigned N)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> ecc_d(min_ecc, max_ecc);
    std::uniform_real_distribution<double> M_d(-1e8, 1e8);

    // We generate the random dataset
    std::vector<double> eccenricities(N);
    std::vector<double> mean_anomalies(N);

    for (auto i = 0u; i < N; ++i) {
        mean_anomalies[i] = M_d(rng_engine);
        eccenricities[i] = ecc_d(rng_engine);
    }

    // We log progress
    fmt::print("{:.2f} min_ecc, {:.2f} max_ecc, on {} data points: ", min_ecc, max_ecc, N);
    std::vector<double> err(N);
    for (auto i = 0u; i < N; ++i) {
        auto res = e2m(m2e(mean_anomalies[i], eccenricities[i]), eccenricities[i]);
        // error is arbitrarily: (|sinM-sinMtrue| +|cosM-cosMtrue|)/2
        err[i] = (std::abs(std::sin(res) - std::sin(mean_anomalies[i]))
                  + std::abs(std::cos(res) - std::cos(mean_anomalies[i])))
                 / 2.;
    }
    auto max_it = max_element(std::begin(err), std::end(err));
    auto min_it = min_element(std::begin(err), std::end(err));
    auto avg = std::accumulate(err.begin(), err.end(), 0.0) / static_cast<double>(err.size());
    fmt::print("{:.3e} avg, {:.3e} min, {:.3e} max\n", avg, *min_it, *max_it);
}

int main()
{
    fmt::print("\nComputes speed at different eccentricity ranges:\n");
    perform_test_speed(0, 0.5, 1000000);
    perform_test_speed(0.5, 0.9, 1000000);
    perform_test_speed(0.9, 0.99, 1000000);
    fmt::print("\nComputes error at different eccentricity ranges:\n");
    perform_test_accuracy(0, 0.5, 100000);
    perform_test_accuracy(0.5, 0.9, 100000);
    perform_test_accuracy(0.9, 0.99, 100000);
}