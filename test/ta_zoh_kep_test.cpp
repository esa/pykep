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

#include <kep3/ta/zoh_kep.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_zoh_kep;
using kep3::ta::get_ta_zoh_kep_cache_dim;
using kep3::ta::get_ta_zoh_kep_var;
using kep3::ta::get_ta_zoh_kep_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_zoh_kep_cache_dim() == 0u);
    auto ta_cached = get_ta_zoh_kep(1e-16);
    REQUIRE(get_ta_zoh_kep_cache_dim() == 1u);
    ta_cached = get_ta_zoh_kep(1e-16);
    REQUIRE(get_ta_zoh_kep_cache_dim() == 1u);
    ta_cached = get_ta_zoh_kep(1e-8);
    REQUIRE(get_ta_zoh_kep_cache_dim() == 2u);

    // The variational integrator.
    REQUIRE(get_ta_zoh_kep_var_cache_dim() == 0u);
    ta_cached = get_ta_zoh_kep_var(1e-16);
    REQUIRE(get_ta_zoh_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_kep_var(1e-16);
    REQUIRE(get_ta_zoh_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_zoh_kep_var(1e-8);
    REQUIRE(get_ta_zoh_kep_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_zoh_kep(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 7);
    REQUIRE(ta_cached.get_pars().size() == 5);

    // Deterministic test vector generated from python module pykep/ta/_zoh_kep.py.
    // Script used: standalone importlib load + numpy rng(seed=123456).
    const std::vector<double> ic = {
        0.27302749971786167,
        -0.23037667022284958,
        -0.9051091571460308,
        0.9105054792315637,
        0.8121018732062479,
        -0.08606089029218555,
        0.920800214245817,
    };
    const std::vector<double> pars = {
        0.03217129137127832,
        0.7368299604928924,
        0.6518130968649866,
        -0.17950291383517422,
        0.57233084801041,
    };
    constexpr double tof = 0.15999279961029858;
    const std::vector<double> ground_truth = {
        0.41469732570119566,
        -0.09762160334579481,
        -0.9066761774933562,
        0.8573492162271017,
        0.8433266702279382,
        0.0643140450442284,
        0.917854327228336,
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
    auto ta_cached = get_ta_zoh_kep(1e-16);
    auto ta_var_cached = get_ta_zoh_kep_var(1e-16);

    REQUIRE(ta_var_cached.is_variational() == true);
    REQUIRE(ta_var_cached.get_dim() == 84);
    REQUIRE(ta_var_cached.get_pars().size() == 5);
    REQUIRE(ta_var_cached.get_vorder() == 1);
    REQUIRE(ta_var_cached.get_vargs().size() == 11);

    // Same initial condition used in python_regression_propagation.
    const std::vector<double> ic = {
        0.27302749971786167,
        -0.23037667022284958,
        -0.9051091571460308,
        0.9105054792315637,
        0.8121018732062479,
        -0.08606089029218555,
        0.920800214245817,
    };
    const std::vector<double> pars = {
        0.03217129137127832,
        0.7368299604928924,
        0.6518130968649866,
        -0.17950291383517422,
        0.57233084801041,
    };
    constexpr double tof = 0.15999279961029858;

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
