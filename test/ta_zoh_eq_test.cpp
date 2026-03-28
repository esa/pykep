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

#include <kep3/ta/zoh_eq.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_zoh_eq;
using kep3::ta::get_ta_zoh_eq_cache_dim;
using kep3::ta::get_ta_zoh_eq_var;
using kep3::ta::get_ta_zoh_eq_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_zoh_eq_cache_dim() == 0u);
    auto ta_cached = get_ta_zoh_eq(1e-16);
    REQUIRE(get_ta_zoh_eq_cache_dim() == 1u);
    ta_cached = get_ta_zoh_eq(1e-16);
    REQUIRE(get_ta_zoh_eq_cache_dim() == 1u);
    ta_cached = get_ta_zoh_eq(1e-8);
    REQUIRE(get_ta_zoh_eq_cache_dim() == 2u);

    // The variational integrator.
    REQUIRE(get_ta_zoh_eq_var_cache_dim() == 0u);
    ta_cached = get_ta_zoh_eq_var(1e-16);
    REQUIRE(get_ta_zoh_eq_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_eq_var(1e-16);
    REQUIRE(get_ta_zoh_eq_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_eq_var(1e-8);
    REQUIRE(get_ta_zoh_eq_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_zoh_eq(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 7);
    REQUIRE(ta_cached.get_pars().size() == 5);

    // Deterministic test vector generated from python module pykep/ta/_zoh_eq.py.
    // Script used: standalone importlib load + numpy rng(seed=424242).
    const std::vector<double> ic = {
        1.0542743379437618,
        -0.710842284896237,
        -0.38121600613133444,
        -0.36323147869252675,
        0.2968969678867671,
        3.3119188796480334,
        0.6910327780298356,
    };
    const std::vector<double> pars = {
        0.04607805345650951,
        -0.1326700998998882,
        0.8277589663109869,
        -0.5451731268911926,
        0.8978481215075556,
    };
    constexpr double tof = 1.1950891197852191;
    const std::vector<double> ground_truth = {
        1.1780477233825897,
        -0.7421098584932024,
        -0.5051860860910615,
        -0.3655690218786264,
        0.3157660329438347,
        5.412410892975656,
        0.6415906340291585,
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
    auto ta_cached = get_ta_zoh_eq(1e-16);
    auto ta_var_cached = get_ta_zoh_eq_var(1e-16);

    REQUIRE(ta_var_cached.is_variational() == true);
    REQUIRE(ta_var_cached.get_dim() == 84);
    REQUIRE(ta_var_cached.get_pars().size() == 5);
    REQUIRE(ta_var_cached.get_vorder() == 1);
    REQUIRE(ta_var_cached.get_vargs().size() == 11);

    // Same initial condition used in propagation.
    const std::vector<double> ic = {
        1.0542743379437618,
        -0.710842284896237,
        -0.38121600613133444,
        -0.36323147869252675,
        0.2968969678867671,
        3.3119188796480334,
        0.6910327780298356,
    };
    const std::vector<double> pars = {
        0.04607805345650951,
        -0.1326700998998882,
        0.8277589663109869,
        -0.5451731268911926,
        0.8978481215075556,
    };
    constexpr double tof = 1.1950891197852191;

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
