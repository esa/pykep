// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <chrono>
#include <fmt/core.h>
#include <fmt/ranges.h>

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

void perform_single_nogradient_speed_test()
{
    std::array<std::array<double, 3>, 2> m_rvs{{{1, 0.1, -0.1}, {0.2, 1, -0.2}}};
    std::array<std::array<double, 3>, 2> m_rvf{{{1.2, -0.1, 0.1}, {-0.2, 1.023, -0.44}}};
    double m_ms = 1;
    double m_mf = m_ms * 13 / 15;
    double m_isp = 1;
    double m_max_thrust = 1;
    double m_cut = 0.5;
    double m_mu = 1;
    double m_tof = 1;
    std::vector<double> m_throttles = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    auto start_con = high_resolution_clock::now();
    auto sf_leg
        = kep3::leg::sims_flanagan(m_rvs, m_ms, m_throttles, m_rvf, m_mf, m_tof, m_max_thrust, m_isp*kep3::G0, m_mu, m_cut);
    auto stop_con = high_resolution_clock::now();
    auto duration_con = duration_cast<microseconds>(stop_con - start_con);
    fmt::print("\nLow-fidelity leg construction: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(duration_con.count()) / 1e6);

    auto start = high_resolution_clock::now();
    [[maybe_unused]] auto mc = sf_leg.compute_mismatch_constraints();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("\nLow-fidelity leg mc: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(duration.count()) / 1e6);

    auto two_start = high_resolution_clock::now();
    auto two_mc = sf_leg.compute_mc_grad();
    auto two_stop = high_resolution_clock::now();
    auto two_duration = duration_cast<microseconds>(two_stop - two_start);
    fmt::print("\nLow-fidelity leg mc_grad: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(two_duration.count()) / 1e6);

    fmt::print("\n\nBelow are the numerical and analytical gradient method calls from the UDPs.\n");

    // Create chromosome
    auto chromosome = m_throttles;
    chromosome.push_back(m_tof);
    chromosome.push_back(m_mf);

    // Create analytical benchmark
    auto bench_udp_a = sf_bench_udp{m_rvs, m_ms, m_rvf, 1, 1, static_cast<unsigned int>(m_throttles.size() / 3), true};

    auto lf_agrad_start = high_resolution_clock::now();
    auto lf_agrad = bench_udp_a.gradient(chromosome);
    auto lf_agrad_stop = high_resolution_clock::now();
    auto lf_agrad_duration = duration_cast<microseconds>(lf_agrad_stop - lf_agrad_start);
    fmt::print("\nLow-fidelity leg analytical gradient: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(lf_agrad_duration.count()) / 1e6);

    // Create numerical benchmark
    auto bench_udp_n = sf_bench_udp{m_rvs, m_ms, m_rvf, 1, 1, static_cast<unsigned int>(m_throttles.size() / 3), false};

    auto lf_ngrad_start = high_resolution_clock::now();
    auto lf_ngrad = bench_udp_n.gradient(chromosome);
    auto lf_ngrad_stop = high_resolution_clock::now();
    auto lf_ngrad_duration = duration_cast<microseconds>(lf_ngrad_stop - lf_ngrad_start);
    fmt::print("\nLow-fidelity leg numerical gradient: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(lf_ngrad_duration.count()) / 1e6);

    fmt::print("\n\nBelow are the numerical and analytical gradient method calls from the pagmo::problems.\n");

    // Create analytical benchmark
    auto bench_udp_a2 = sf_bench_udp{m_rvs, m_ms, m_rvf, 1, 1, static_cast<unsigned int>(m_throttles.size() / 3), true};
    pagmo::problem prob_a{bench_udp_a2};

    auto lf_agrad_start2 = high_resolution_clock::now();
    auto lf_agrad2 = prob_a.gradient(chromosome);
    auto lf_agrad_stop2 = high_resolution_clock::now();
    auto lf_agrad_duration2 = duration_cast<microseconds>(lf_agrad_stop2 - lf_agrad_start2);
    fmt::print("\nPagmo problem Low-fidelity leg analytical gradient: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(lf_agrad_duration2.count()) / 1e6);

    // Create numerical benchmark
    auto bench_udp_n2 = sf_bench_udp{m_rvs, m_ms, m_rvf, 1, 1, static_cast<unsigned int>(m_throttles.size() / 3), false};
    pagmo::problem prob_n{bench_udp_n2};

    auto lf_ngrad_start2 = high_resolution_clock::now();
    auto lf_ngrad2 = prob_n.gradient(chromosome);
    auto lf_ngrad_stop2 = high_resolution_clock::now();
    auto lf_ngrad_duration2 = duration_cast<microseconds>(lf_ngrad_stop2 - lf_ngrad_start2);
    fmt::print("\nPagmo problem Low-fidelity leg numerical gradient: {} nseg - timing: {}", m_throttles.size() / 3,
               static_cast<double>(lf_ngrad_duration2.count()) / 1e6);
}

int main()
{
    perform_single_nogradient_speed_test();
}