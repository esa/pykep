/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <array>
#include <chrono>
#include <random>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/lambert_problem.hpp>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

int main()
{
    // Number of trials
    const unsigned trials = 50000u;

    std::array<std::array<double, 3>, trials> r1s{}, r2s{};
    std::array<double, trials> tof{};
    std::array<bool, trials> cw{};
    std::array<double, trials> mu{};

    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(122012203u);
    std::uniform_int_distribution<unsigned> cw_d(0, 1);
    std::uniform_real_distribution<double> r_d(-2, 2);
    std::uniform_real_distribution<double> tof_d(2., 40.);
    std::uniform_real_distribution<double> mu_d(0.9, 1.1);
    unsigned revs_max = 20u;
    unsigned long count = 0lu;

    for (auto i = 0u; i < trials; ++i) {
        // 1 - generate a random problem geometry
        r1s[i][0] = r_d(rng_engine);
        r1s[i][1] = r_d(rng_engine);
        r1s[i][2] = r_d(rng_engine);
        r2s[i][0] = r_d(rng_engine);
        r2s[i][1] = r_d(rng_engine);
        r2s[i][2] = r_d(rng_engine);
        tof[i] = tof_d(rng_engine);
        cw[i] = static_cast<bool>(cw_d(rng_engine));
        mu[i] = mu_d(rng_engine);
    }

    auto start = high_resolution_clock::now();
    for (auto i = 0u; i < trials; ++i) {
        // 2 - Solve the lambert problem
        kep3::lambert_problem lp(r1s[i], r2s[i], tof[i], mu[i], cw[i], revs_max);
        count = count + lp.get_v0().size();
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("Lambert:\n{} solutions computed in {:.3f}s\n", count, (static_cast<double>(duration.count()) / 1e6));
    fmt::print("Projected number of solutions per second: {}\n",
               static_cast<double>(count) / ((static_cast<double>(duration.count()) / 1e6)));

    count = 0; // reset counter
    start = high_resolution_clock::now();
    for (auto i = 0u; i < trials; ++i) {
        // 3 - Solve the lambert problem for the singe rev case
        kep3::lambert_problem lp(r1s[i], r2s[i], tof[i], mu[i], cw[i], 0u);
        count += 1u;
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    fmt::print("\nLambert (0 revs only):\n{} solutions computed in {:.3f}s\n", count,
               (static_cast<double>(duration.count()) / 1e6));
    fmt::print("Projected number of solutions per second: {}\n",
               static_cast<double>(count) / ((static_cast<double>(duration.count()) / 1e6)));
}