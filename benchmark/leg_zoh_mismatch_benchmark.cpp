// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <fmt/core.h>

#include <kep3/leg/zoh.hpp>
#include <kep3/ta/zoh_kep.hpp>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

namespace
{

struct zoh_case {
    std::vector<double> state0;
    std::vector<double> controls;
    std::vector<double> state1;
    std::vector<double> tgrid;
};

std::array<double, 3> random_unit_vector(std::mt19937 &rng)
{
    std::normal_distribution<double> n01(0.0, 1.0);
    std::array<double, 3> v = {n01(rng), n01(rng), n01(rng)};
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n <= std::numeric_limits<double>::epsilon()) {
        return {1.0, 0.0, 0.0};
    }
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    return v;
}

std::vector<zoh_case> generate_cases(unsigned ncases, unsigned nseg, std::uint32_t seed)
{
    std::mt19937 rng(seed);

    std::uniform_real_distribution<double> r_dist(-0.3, 0.3);
    std::uniform_real_distribution<double> v_dist(-0.35, 0.35);
    std::uniform_real_distribution<double> m_dist(0.85, 1.15);
    std::uniform_real_distribution<double> thrust_dist(0.0, 0.03);
    std::uniform_real_distribution<double> tof_dist(0.3, 2.0);

    std::vector<zoh_case> out;
    out.reserve(ncases);

    for (unsigned i = 0u; i < ncases; ++i) {
        zoh_case c;
        c.state0.resize(7u);
        c.state1.resize(7u);
        c.controls.resize(static_cast<std::size_t>(4u * nseg));
        c.tgrid.resize(static_cast<std::size_t>(nseg + 1u));

        // Keep states near non-dimensional heliocentric magnitudes used in tests.
        c.state0[0] = 1.0 + r_dist(rng);
        c.state0[1] = r_dist(rng);
        c.state0[2] = r_dist(rng);
        c.state0[3] = v_dist(rng);
        c.state0[4] = 1.0 + v_dist(rng);
        c.state0[5] = v_dist(rng);
        c.state0[6] = m_dist(rng);

        c.state1[0] = 1.0 + r_dist(rng);
        c.state1[1] = r_dist(rng);
        c.state1[2] = r_dist(rng);
        c.state1[3] = v_dist(rng);
        c.state1[4] = 1.0 + v_dist(rng);
        c.state1[5] = v_dist(rng);
        c.state1[6] = m_dist(rng);

        const double tof = tof_dist(rng);
        const double dt = tof / static_cast<double>(nseg);
        for (unsigned j = 0u; j <= nseg; ++j) {
            c.tgrid[j] = static_cast<double>(j) * dt;
        } 

        for (unsigned j = 0u; j < nseg; ++j) {
            const auto dir = random_unit_vector(rng);
            const double thrust = thrust_dist(rng);
            const auto idx = static_cast<std::size_t>(4u * j);
            c.controls[idx] = thrust;
            c.controls[idx + 1u] = dir[0];
            c.controls[idx + 2u] = dir[1];
            c.controls[idx + 3u] = dir[2];
        }

        out.push_back(std::move(c));
    }

    return out;
}

std::vector<kep3::leg::zoh> make_legs(const std::vector<zoh_case> &cases, double cut)
{
    auto ta = kep3::ta::get_ta_zoh_kep(1e-14);
    auto ta_var = kep3::ta::get_ta_zoh_kep_var(1e-8);

    // Parameter c = 1 / v_eff in non-dimensional units.
    *(ta.get_pars_data() + 4) = 1.0 / 3000.0;
    *(ta_var.get_pars_data() + 4) = 1.0 / 3000.0;

    std::vector<kep3::leg::zoh> legs;
    legs.reserve(cases.size());

    for (const auto &c : cases) {
        legs.emplace_back(c.state0, c.controls, c.state1, c.tgrid, cut, std::make_pair(ta, ta_var));
    }

    return legs;
}

double benchmark_mc(const std::vector<kep3::leg::zoh> &legs, unsigned repeats)
{
    volatile double checksum = 0.0;

    const auto t0 = high_resolution_clock::now();
    for (unsigned r = 0u; r < repeats; ++r) {
        for (const auto &leg : legs) {
            const auto mc = leg.compute_mismatch_constraints();
            checksum += mc[0] + mc[6];
        }
    }
    const auto t1 = high_resolution_clock::now();

    const auto elapsed_us = static_cast<double>(duration_cast<microseconds>(t1 - t0).count());
    const double ncalls = static_cast<double>(repeats) * static_cast<double>(legs.size());

    fmt::print("mismatch constraints : total {:.6f} s | {:.3f} us/call\n", elapsed_us / 1e6,
               elapsed_us / ncalls);

    return elapsed_us;
}

double benchmark_mc_grad(const std::vector<kep3::leg::zoh> &legs, unsigned repeats)
{
    volatile double checksum = 0.0;

    const auto t0 = high_resolution_clock::now();
    for (unsigned r = 0u; r < repeats; ++r) {
        for (const auto &leg : legs) {
            const auto [dmc_dx0, dmc_dx1, dmc_du, dmc_dtgrid] = leg.compute_mc_grad();
            checksum += dmc_dx0[0] + dmc_dx1[0] + dmc_du[0] + dmc_dtgrid[0];
        }
    }
    const auto t1 = high_resolution_clock::now();

    const auto elapsed_us = static_cast<double>(duration_cast<microseconds>(t1 - t0).count());
    const double ncalls = static_cast<double>(repeats) * static_cast<double>(legs.size());

    fmt::print("mismatch gradient    : total {:.6f} s | {:.3f} us/call\n", elapsed_us / 1e6,
               elapsed_us / ncalls);

    return elapsed_us;
}

void run_benchmark(unsigned ncases, unsigned nseg, unsigned repeats, std::uint32_t seed)
{
    constexpr double cut = 0.5;

    fmt::print("\nZOH mismatch benchmark (Keplerian dynamics)\n");
    fmt::print("cases: {} | nseg: {} | repeats: {} | seed: {}\n", ncases, nseg, repeats, seed);

    const auto cases = generate_cases(ncases, nseg, seed);
    const auto legs = make_legs(cases, cut);

    benchmark_mc(legs, repeats);
    benchmark_mc_grad(legs, repeats);
    fmt::print("\n");
}

} // namespace

int main()
{
    run_benchmark(200u, 20u, 5u, 424242u);
    run_benchmark(200u, 40u, 5u, 424242u);
    run_benchmark(200u, 80u, 5u, 424242u);
}
