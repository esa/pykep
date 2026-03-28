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

#include <kep3/ta/zoh_ss.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_zoh_ss;
using kep3::ta::get_ta_zoh_ss_cache_dim;
using kep3::ta::get_ta_zoh_ss_var;
using kep3::ta::get_ta_zoh_ss_var_cache_dim;

TEST_CASE("caches")
{
    REQUIRE(get_ta_zoh_ss_cache_dim() == 0u);
    auto ta_cached = get_ta_zoh_ss(1e-16);
    REQUIRE(get_ta_zoh_ss_cache_dim() == 1u);
    ta_cached = get_ta_zoh_ss(1e-16);
    REQUIRE(get_ta_zoh_ss_cache_dim() == 1u);
    ta_cached = get_ta_zoh_ss(1e-8);
    REQUIRE(get_ta_zoh_ss_cache_dim() == 2u);

    REQUIRE(get_ta_zoh_ss_var_cache_dim() == 0u);
    ta_cached = get_ta_zoh_ss_var(1e-16);
    REQUIRE(get_ta_zoh_ss_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_ss_var(1e-16);
    REQUIRE(get_ta_zoh_ss_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_ss_var(1e-8);
    REQUIRE(get_ta_zoh_ss_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_zoh_ss(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 6);
    REQUIRE(ta_cached.get_pars().size() == 3);

    const std::vector<double> ic = {
        0.8,
        -0.4,
        0.3,
        0.2,
        0.9,
        -0.1,
    };
    const std::vector<double> pars = {
        0.25,
        -1.1,
        0.04,
    };
    constexpr double tof = 0.75;
    const std::vector<double> ground_truth = {
        0.6167384690293904,
        0.3246225470281653,
        0.12203176418176016,
        -0.7886148445687173,
        0.8714494740768439,
        -0.3751503867134815,
    };

    taylor_adaptive<double> ta(ta_cached);
    ta.set_time(0.);
    std::copy(ic.begin(), ic.end(), ta.get_state_data());
    std::copy(pars.begin(), pars.end(), ta.get_pars_data());

    auto out = ta.propagate_until(tof);
    REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
    REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
}

TEST_CASE("variational_propagation")
{
    auto ta_cached = get_ta_zoh_ss(1e-16);
    auto ta_var_cached = get_ta_zoh_ss_var(1e-16);

    REQUIRE(ta_var_cached.is_variational() == true);
    REQUIRE(ta_var_cached.get_dim() == 54);
    REQUIRE(ta_var_cached.get_pars().size() == 3);
    REQUIRE(ta_var_cached.get_vorder() == 1);
    REQUIRE(ta_var_cached.get_vargs().size() == 8);

    const std::vector<double> ic = {
        0.8,
        -0.4,
        0.3,
        0.2,
        0.9,
        -0.1,
    };
    const std::vector<double> pars = {
        0.25,
        -1.1,
        0.04,
    };
    constexpr double tof = 0.75;

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
                std::vector<double>(ta_var.get_state().begin(), ta_var.get_state().begin() + 6), ta.get_state())
            <= 1e-13);
}
