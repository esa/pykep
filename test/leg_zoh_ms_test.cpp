// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <heyoka/expression.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <heyoka/kw.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/leg/zoh_ms.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include "catch.hpp"
#include "test_helpers.hpp"

namespace
{

using integrator = heyoka::taylor_adaptive<double>;
using integrator_pair = std::pair<integrator, std::optional<integrator>>;

// Builds x_i' = (i + 1) * (x_i + sum(pars)).
integrator make_integrator(unsigned dim, unsigned npars)
{
    std::vector<std::pair<heyoka::expression, heyoka::expression>> sys;
    auto control = heyoka::par[0];
    for (unsigned i = 1u; i < npars; ++i) {
        control += heyoka::par[i];
    }
    for (unsigned i = 0u; i < dim; ++i) {
        const auto var = heyoka::make_vars("x" + std::to_string(i));
        const double rate = static_cast<double>(i + 1u);
        sys.emplace_back(var, rate * (var + control));
    }
    return integrator{std::move(sys)};
}

// Uniform grid, state dimension 1 and control dimension 2.
kep3::leg::zoh_ms make_test_leg_12()
{
    const auto ta = make_integrator(1u, 2u);
    return {{0., 1., 2.}, {10., 11., 12., 13.}, {0., 1., 2.}, 0.5, {ta, std::nullopt}, std::nullopt, 1u, 2u};
}

// Nonuniform grid, state dimension 1 and control dimension 1.
kep3::leg::zoh_ms make_test_leg_11()
{
    const auto ta = make_integrator(1u, 1u);
    return {{1., 99., 4., 8.}, {0.5, 1., 1.5}, {0., 0.5, 2., 5.}, 0.5, {ta, std::nullopt}, std::nullopt, 1u, 1u};
}

// Nonuniform grid, state dimension 2 and control dimension 1.
kep3::leg::zoh_ms make_test_leg_21()
{
    const std::vector<double> states{0., 1., 2., 3., 4.,5., 6., 7.};
    const std::vector<double> controls{10., 11., 12.};
    const std::vector<double> tgrid{0., 0.5, 2., 5.};
    auto ta = make_integrator(2u, 1u);
    return {states, controls, tgrid, 0.5, {ta, std::nullopt}, std::nullopt,2u, 1u};
}

} // namespace

TEST_CASE("zoh_ms constructor")
{
    // We test constructor validation on a uniform-grid 1-state, 1-control leg and malformed integrator variants.
    const std::vector<double> states(3u, 0.0);
    const std::vector<double> controls(2u, 0.0);
    const std::vector<double> tgrid{0., 1., 2.};
    const auto ta = make_integrator(1u, 1u);
    const integrator_pair nominal_only{ta, std::nullopt};
    const integrator_pair wrong_ta_dimension{make_integrator(2u, 1u), std::nullopt};
    const integrator_pair wrong_var_dimension{ta, ta};
    const integrator_pair wrong_var_pars{ta, make_integrator(3u, 2u)};
    const auto make_leg = [](const std::vector<double> &leg_states, const std::vector<double> &leg_controls,
                             const std::vector<double> &leg_tgrid, double cut, const integrator_pair &tas,
                             unsigned dim_dynamics = 1u, unsigned dim_controls = 1u) {
        return kep3::leg::zoh_ms{leg_states, leg_controls, leg_tgrid, cut, tas, std::nullopt, dim_dynamics,
                                 dim_controls};
    };

    auto leg = make_leg(states, controls, tgrid, 0.5, nominal_only);

    REQUIRE(leg.get_states() == states);
    REQUIRE(leg.get_controls() == controls);
    REQUIRE(leg.get_tgrid() == tgrid);
    REQUIRE(leg.get_cut() == 0.5);
    REQUIRE(leg.get_dim_dynamics() == 1u);
    REQUIRE(leg.get_dim_controls() == 1u);
    REQUIRE(leg.get_nseg() == 2u);
    REQUIRE(leg.get_nseg_fwd() == 1u);
    REQUIRE(leg.get_nseg_bck() == 1u);
    REQUIRE_FALSE(leg.has_ta_var());

    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 0.5, nominal_only, 0u), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 0.5, nominal_only, 1u, 0u), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, -0.1, nominal_only), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 1.1, nominal_only), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(std::vector<double>(2u, 0.0), controls, tgrid, 0.5, nominal_only), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 0.5, wrong_ta_dimension), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(std::vector<double>(2u, 0.0), std::vector<double>(2u, 0.0),
                               std::vector<double>{0., 1.}, 0.5, nominal_only, 1u, 2u),
                      std::logic_error);
    REQUIRE_THROWS_AS(make_leg(std::vector<double>(2u, 0.0), std::vector<double>(3u, 0.0),
                               std::vector<double>{0., 1.}, 0.5, nominal_only, 1u, 2u),
                      std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, std::vector<double>{0., 2.}, 0.5, nominal_only), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 0.5, wrong_var_dimension), std::logic_error);
    REQUIRE_THROWS_AS(make_leg(states, controls, tgrid, 0.5, wrong_var_pars), std::logic_error);
}

TEST_CASE("zoh_ms setters and getters")
{
    // We test that setters update their getters and reject incompatible
    // inputs without changing stored data on a uniform-grid 1-state, 2-control leg.
    auto leg = make_test_leg_12();
    const std::vector<double> states{0., 1., 2.};
    const std::vector<double> controls{10., 11., 12., 13.};
    const std::vector<double> tgrid{0., 1., 2.};
    const std::vector<double> new_states{2., 3., 4.};
    const std::vector<double> new_controls{20., 21., 22., 23.};
    const std::vector<double> new_tgrid{1., 2., 4.};

    REQUIRE_THROWS_AS(leg.set_states({0., 1.}), std::logic_error);
    REQUIRE(leg.get_states() == states);
    leg.set_states(new_states);
    REQUIRE(leg.get_states() == new_states);

    REQUIRE_THROWS_AS(leg.set_controls({10., 11.}), std::logic_error);
    REQUIRE(leg.get_controls() == controls);
    leg.set_controls(new_controls);
    REQUIRE(leg.get_controls() == new_controls);

    REQUIRE_THROWS_AS(leg.set_tgrid({0., 1.}), std::logic_error);
    REQUIRE(leg.get_tgrid() == tgrid);
    leg.set_tgrid(new_tgrid);
    REQUIRE(leg.get_tgrid() == new_tgrid);

    REQUIRE_THROWS_AS(leg.set_cut(-0.1), std::logic_error);
    REQUIRE_THROWS_AS(leg.set_cut(1.1), std::logic_error);
    REQUIRE(leg.get_cut() == 0.5);
    leg.set_cut(0.75);
    REQUIRE(leg.get_cut() == 0.75);
    leg.set_max_steps(8u);
    REQUIRE(leg.get_max_steps() == std::optional<unsigned>{8u});

    REQUIRE(leg.get_nseg() == 2u);
    REQUIRE(leg.get_nseg_fwd() == 1u);
    REQUIRE(leg.get_nseg_bck() == 1u);

    // Unsuccesful set attempts
    REQUIRE_THROWS_AS(leg.set(states, std::vector<double>(5u, 0.0), tgrid), std::logic_error);
    REQUIRE_THROWS_AS(leg.set(states, std::vector<double>(6u, 0.0), tgrid), std::logic_error);
    REQUIRE_THROWS_AS(leg.set(states, controls, tgrid, -0.1), std::logic_error);
    REQUIRE(leg.get_states() == new_states);
    REQUIRE(leg.get_controls() == new_controls);
    REQUIRE(leg.get_tgrid() == new_tgrid);
    REQUIRE(leg.get_cut() == 0.75);
    REQUIRE(leg.get_max_steps() == std::optional<unsigned>{8u});

    // Successful set attempts
    const std::vector<double> resized_states{0., 1., 2., 3.};
    const std::vector<double> resized_controls{30., 31., 32., 33., 34., 35.};
    const std::vector<double> resized_tgrid{0., 1., 2., 3.};
    leg.set(resized_states, resized_controls, resized_tgrid, std::nullopt, 4u);
    REQUIRE(leg.get_states() == resized_states);
    REQUIRE(leg.get_controls() == resized_controls);
    REQUIRE(leg.get_tgrid() == resized_tgrid);
    REQUIRE(leg.get_cut() == 0.75);
    REQUIRE(leg.get_max_steps() == std::optional<unsigned>{4u});
    REQUIRE(leg.get_nseg() == 3u);
    REQUIRE(leg.get_nseg_fwd() == 2u);
    REQUIRE(leg.get_nseg_bck() == 1u);

    // Another successful set attempt
    leg.set(states, controls, tgrid, 0.25);
    REQUIRE(leg.get_states() == states);
    REQUIRE(leg.get_controls() == controls);
    REQUIRE(leg.get_tgrid() == tgrid);
    REQUIRE(leg.get_cut() == 0.25);
    REQUIRE_FALSE(leg.get_max_steps().has_value());
    REQUIRE(leg.get_nseg() == 2u);
    REQUIRE(leg.get_nseg_fwd() == 0u);
    REQUIRE(leg.get_nseg_bck() == 2u);
}

TEST_CASE("compute_defects") {
    // We test forward and backward defects on a nonuniform-grid 2-state, 1-control leg.
    auto leg = make_test_leg_21();
    const auto defects = leg.compute_defects();
    const std::vector<double> expected{4.487212707001282, 16.901100113049495, 9.653047597773553, 13.203406906114177, 15.10383276937845, 16.952903708643585};
    REQUIRE(kep3_tests::L_infinity_norm_rel(defects, expected) < 1e-12);
}

TEST_CASE("get_state_info")
{
    // We test per-segment sampling order on a nonuniform-grid 1-state, 1-control leg.
    auto leg = make_test_leg_11();
    const auto [state_fwd, state_bck, success] = leg.get_state_info(3u);

    REQUIRE(success);
    REQUIRE(state_fwd.size() == 1u);
    REQUIRE(state_fwd[0].size() == 3u);
    REQUIRE(state_bck.size() == 2u);
    REQUIRE(state_bck[0].size() == 3u);
    REQUIRE(state_bck[1].size() == 3u);

    const std::vector<double> actual_fwd{state_fwd[0][0][0], state_fwd[0][1][0], state_fwd[0][2][0]};
    const std::vector<double> actual_bck{state_bck[0][0][0], state_bck[0][1][0], state_bck[0][2][0],
                                         state_bck[1][0][0], state_bck[1][1][0], state_bck[1][2][0]};
    const std::vector<double> expected_fwd{1., 1.426038125031612, 1.9730819060501923};
    const std::vector<double> expected_bck{8., 0.6197365214100832, -1.0270228505052925, 4.,
                                           1.3618327637050736, 0.11565080074214906};
    REQUIRE(kep3_tests::L_infinity_norm_rel(actual_fwd, expected_fwd) < 1e-12);
    REQUIRE(kep3_tests::L_infinity_norm_rel(actual_bck, expected_bck) < 1e-12);
}

TEST_CASE("set_initial_guess propagates with stored controls")
{
    // We test control-driven propagation on a nonuniform-grid 1-state, 1-control leg.
    auto leg = make_test_leg_11();
    const std::vector<double> controls{0.5, 1., 1.5};
    const std::vector<double> tgrid{0., 0.5, 2., 5.};

    leg.set_initial_guess();

    const std::vector<double> expected_states{1., 1.9730819060501923, -1.0270228505052925, 8.};
    REQUIRE(kep3_tests::L_infinity_norm_rel(leg.get_states(), expected_states) < 1e-12);
    REQUIRE(leg.get_controls() == controls);
    REQUIRE(leg.get_tgrid() == tgrid);
}

TEST_CASE("set_initial_guess propagates ballistically")
{
    // We test ballistic propagation on a nonuniform-grid 1-state, 1-control leg.
    auto leg = make_test_leg_11();
    const std::vector<double> controls{0.5, 1., 1.5};
    const std::vector<double> tgrid{0., 0.5, 2., 5.};

    leg.set_initial_guess(true);

    const std::vector<double> expected_states{1., 1.6487212707001282, 0.39829654694291156, 8.};
    REQUIRE(kep3_tests::L_infinity_norm_rel(leg.get_states(), expected_states) < 1e-12);
    REQUIRE(leg.get_controls() == controls);
    REQUIRE(leg.get_tgrid() == tgrid);
}