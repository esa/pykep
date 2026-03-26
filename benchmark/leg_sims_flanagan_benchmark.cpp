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
#include <fmt/ranges.h>

#include <xtensor/containers/xadapt.hpp>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "leg_sims_flanagan_udp_bench.hpp"

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void perform_convergence_benchmark(unsigned N, unsigned nseg)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> dv_pert_random(0., 0.1);
    std::uniform_real_distribution<double> mass_random(1.0, 1.2);
    std::uniform_real_distribution<double> tof_random(kep3::pi / 12, 2 * kep3::pi);
    std::uniform_real_distribution<double> ts_random(2170, 2200);

    // Create test leg for initial conditions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    // auto rvs = earth.eph(1000);
    // auto rvf = jupiter.eph(1000);
    int count_n = 0;
    int count_a = 0;
    for (decltype(N) i = 0; i < N; ++i) {
        auto rvs = earth.eph(ts_random(rng_engine));
        auto rvf = jupiter.eph(ts_random(rng_engine));
        double tof_ic = tof_random(rng_engine);
        double mu = 1;
        rvs[0][0] /= kep3::AU;
        rvs[0][1] /= kep3::AU;
        rvs[0][2] /= kep3::AU;
        rvf[0][0] /= kep3::AU;
        rvf[0][1] /= kep3::AU;
        rvf[0][2] /= kep3::AU;
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof_ic, mu};

        // Create HF legs
        std::array<std::array<double, 3>, 2> rvs_udp_ic = {{{lp.get_r0()[0], lp.get_r0()[1], lp.get_r0()[2]},
                                                            {lp.get_v0()[0][0], lp.get_v0()[0][1], lp.get_v0()[0][2]}}};
        std::array<std::array<double, 3>, 2> rvf_udp_ic
            = {{{lp.get_r1()[0], lp.get_r1()[1], lp.get_r1()[2]},
                {lp.get_v1()[0][0] + dv_pert_random(rng_engine), lp.get_v1()[0][1] + dv_pert_random(rng_engine),
                 lp.get_v1()[0][2] + dv_pert_random(rng_engine)}}};
        // double mass = 1;
        double mass = mass_random(rng_engine);
        double max_thrust = 1;
        double isp = 1;
        auto bench_udp_a = sf_bench_udp{rvs_udp_ic, mass, rvf_udp_ic, max_thrust, isp, nseg, true};
        auto bench_udp_n = sf_bench_udp{rvs_udp_ic, mass, rvf_udp_ic, max_thrust, isp, nseg, false};
        pagmo::problem prob_a{bench_udp_a};
        pagmo::problem prob_n{bench_udp_n};
        prob_a.set_c_tol(1e-8);
        prob_n.set_c_tol(1e-8);

        // We construct the same random population
        pagmo::population pop_a{prob_a, 1u};
        pagmo::population pop_n{prob_n};
        pop_n.push_back(pop_a.get_x()[0]);

        // We construct the solver
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(0);
        uda.set_xtol_rel(0);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};

        // We solve first a
        pop_a = algo.evolve(pop_a);
        if (prob_a.feasibility_f(pop_a.get_f()[0])) {
            count_a++;
            std::cout << "." << std::flush;
        } else {
            std::cout << "x" << std::flush;
        }
        // then n
        pop_n = algo.evolve(pop_n);
        if (prob_n.feasibility_f(pop_n.get_f()[0])) {
            count_n++;
            std::cout << "." << std::flush;
        } else {
            std::cout << "x" << std::flush;
        }
    }
    fmt::print("\n{} nseg - success rates: analytical {}/{} - numerical {}/{}\n", nseg, count_a, N, count_n, N);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void perform_speed_benchmark(unsigned N, unsigned nseg, unsigned pop_size)
{
    //
    // Engines
    //
    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    //
    // Distributions
    //
    std::uniform_real_distribution<double> dv_pert_random(0., 0.1);
    std::uniform_real_distribution<double> mass_random(1.0, 1.2);
    std::uniform_real_distribution<double> tof_random(kep3::pi / 12, 2 * kep3::pi);
    std::uniform_real_distribution<double> ts_random(2170, 2200);

    // Create test leg for initial conditions
    kep3::udpla::vsop2013 udpla_earth("earth_moon", 1e-2);
    kep3::udpla::vsop2013 udpla_jupiter("jupiter", 1e-2);
    kep3::planet earth{udpla_earth};
    kep3::planet jupiter{udpla_jupiter};
    double count_n = 0;
    double count_a = 0;
    for (decltype(N) i = 0; i < N; ++i) {
        auto rvs = earth.eph(ts_random(rng_engine));
        auto rvf = jupiter.eph(ts_random(rng_engine));
        double tof_ic = tof_random(rng_engine);
        double mu = 1;
        rvs[0][0] /= kep3::AU;
        rvs[0][1] /= kep3::AU;
        rvs[0][2] /= kep3::AU;
        rvf[0][0] /= kep3::AU;
        rvf[0][1] /= kep3::AU;
        rvf[0][2] /= kep3::AU;
        const kep3::lambert_problem lp{rvs[0], rvf[0], tof_ic, mu};

        // Create HF legs
        std::array<std::array<double, 3>, 2> rvs_udp_ic = {{{lp.get_r0()[0], lp.get_r0()[1], lp.get_r0()[2]},
                                                            {lp.get_v0()[0][0], lp.get_v0()[0][1], lp.get_v0()[0][2]}}};
        std::array<std::array<double, 3>, 2> rvf_udp_ic
            = {{{lp.get_r1()[0], lp.get_r1()[1], lp.get_r1()[2]},
                {lp.get_v1()[0][0] + dv_pert_random(rng_engine), lp.get_v1()[0][1] + dv_pert_random(rng_engine),
                 lp.get_v1()[0][2] + dv_pert_random(rng_engine)}}};
        double mass = mass_random(rng_engine);
        double max_thrust = 1;
        double isp = 1;
        auto bench_udp_a = sf_bench_udp{rvs_udp_ic, mass, rvf_udp_ic, max_thrust, isp, nseg, true};
        auto bench_udp_n = sf_bench_udp{rvs_udp_ic, mass, rvf_udp_ic, max_thrust, isp, nseg, false};
        pagmo::problem prob_a{bench_udp_a};
        pagmo::problem prob_n{bench_udp_n};
        prob_a.set_c_tol(1e-8);
        prob_n.set_c_tol(1e-8);

        // We construct the same random population
        pagmo::population pop{prob_a, pop_size};

        // We construct the solver
        pagmo::nlopt uda{"slsqp"};
        uda.set_xtol_abs(0);
        uda.set_xtol_rel(0);
        uda.set_ftol_abs(0);
        uda.set_maxeval(1000);
        pagmo::algorithm algo{uda};

        // First we time the analytical gradients
        auto start = high_resolution_clock::now();
        for (decltype(pop_size) j = 0u; j < pop_size; ++j) {
            prob_a.gradient(pop.get_x()[j]);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        count_a += static_cast<double>(duration.count()) / 1e6;

        // then the numerical ones
        auto start2 = high_resolution_clock::now();
        for (decltype(pop_size) j = 0u; j < pop_size; ++j) {
            prob_n.gradient(pop.get_x()[j]);
        }
        auto stop2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(stop2 - start2);
        count_n += static_cast<double>(duration2.count()) / 1e6;
    }
    fmt::print("{} nseg - timing: analytical {} - numerical {}\n", nseg, count_a, count_n);
}

int main()
{
    fmt::print("\nComputes the same analytical and numerical gradients and tests for speed:\n");
    perform_speed_benchmark(100, 5, 10);
    perform_speed_benchmark(100, 10, 10);
    perform_speed_benchmark(100, 20, 10);
    perform_speed_benchmark(100, 40, 10);

    // // performing tests
    fmt::print("\nSolves the same optimization problems with and without analytical gradients:\n");
    perform_convergence_benchmark(100, 5);
    perform_convergence_benchmark(100, 10);
    perform_convergence_benchmark(100, 15);

    fmt::print("\n");
}