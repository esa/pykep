// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <vector>

#include <heyoka/taylor.hpp>

#include <kep3/ta/zoh_cr3bp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_zoh_cr3bp;
using kep3::ta::get_ta_zoh_cr3bp_cache_dim;
using kep3::ta::get_ta_zoh_cr3bp_var;
using kep3::ta::get_ta_zoh_cr3bp_var_cache_dim;

TEST_CASE("caches")
{
    REQUIRE(get_ta_zoh_cr3bp_cache_dim() == 0u);
    auto ta_cached = get_ta_zoh_cr3bp(1e-16);
    REQUIRE(get_ta_zoh_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_zoh_cr3bp(1e-16);
    REQUIRE(get_ta_zoh_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_zoh_cr3bp(1e-8);
    REQUIRE(get_ta_zoh_cr3bp_cache_dim() == 2u);

    REQUIRE(get_ta_zoh_cr3bp_var_cache_dim() == 0u);
    ta_cached = get_ta_zoh_cr3bp_var(1e-16);
    REQUIRE(get_ta_zoh_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_cr3bp_var(1e-16);
    REQUIRE(get_ta_zoh_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_cr3bp_var(1e-8);
    REQUIRE(get_ta_zoh_cr3bp_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_zoh_cr3bp(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 7);
    REQUIRE(ta_cached.get_pars().size() == 6);

    const std::vector<double> ic = {
        0.7505822770782307,
        -0.5633872910541109,
        -0.3617217251653817,
        -0.6180765264054857,
        -0.6418676529509515,
        -0.0933954683734568,
        1.0450135084271823,
    };
    const std::vector<double> pars = {
        0.001342869355874922,
        -0.7647203372569936,
        -0.5324337043971957,
        -0.3629285827920274,
        0.4265354685403151,
        0.4863813830217882,
    };
    constexpr double tof = 1.9473133692654372;
    const std::vector<double> ground_truth = {
        0.01822831627059477,
        -0.1520005589632799,
        -0.38639917403339696,
        -0.7216404879147444,
        0.6308233815628889,
        -0.06657963277600498,
        1.0438981235300238,
    };

    taylor_adaptive<double> ta(ta_cached);
    ta.set_time(0.);
    std::copy(ic.begin(), ic.end(), ta.get_state_data());
    std::copy(pars.begin(), pars.end(), ta.get_pars_data());

    auto out = ta.propagate_until(tof);
    REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
    REQUIRE(ta.get_state().size() == ground_truth.size());
    for (auto i = 0u; i < ground_truth.size(); ++i) {
        REQUIRE(ta.get_state()[i] == Approx(ground_truth[i]).epsilon(1e-11));
    }
}

TEST_CASE("variational_propagation")
{
    auto ta_cached = get_ta_zoh_cr3bp(1e-16);
    auto ta_var_cached = get_ta_zoh_cr3bp_var(1e-16);

    REQUIRE(ta_var_cached.is_variational() == true);
    REQUIRE(ta_var_cached.get_dim() == 84);
    REQUIRE(ta_var_cached.get_pars().size() == 6);
    REQUIRE(ta_var_cached.get_vorder() == 1);
    REQUIRE(ta_var_cached.get_vargs().size() == 11);

    const std::vector<double> ic = {
        0.7505822770782307,
        -0.5633872910541109,
        -0.3617217251653817,
        -0.6180765264054857,
        -0.6418676529509515,
        -0.0933954683734568,
        1.0450135084271823,
    };
    const std::vector<double> pars = {
        0.001342869355874922,
        -0.7647203372569936,
        -0.5324337043971957,
        -0.3629285827920274,
        0.4265354685403151,
        0.4863813830217882,
    };
    constexpr double tof = 1.9473133692654372;

    taylor_adaptive<double> ta(ta_cached);
    ta.set_time(0.);
    std::copy(ic.begin(), ic.end(), ta.get_state_data());
    std::copy(pars.begin(), pars.end(), ta.get_pars_data());

    taylor_adaptive<double> ta_var(ta_var_cached);
    ta_var.set_time(0.);
    std::copy(ic.begin(), ic.end(), ta_var.get_state_data());
    std::copy(pars.begin(), pars.end(), ta_var.get_pars_data());

    auto out = ta.propagate_until(tof);
    auto out_var = ta_var.propagate_until(tof);

    REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
    REQUIRE(std::get<0>(out_var) == taylor_outcome::time_limit);

    REQUIRE(kep3_tests::L_infinity_norm_rel(
                std::vector<double>(ta_var.get_state().begin(), ta_var.get_state().begin() + 7), ta.get_state())
            <= 1e-13);
}
