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
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <heyoka/kw.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/leg/zoh_ms.hpp>
#include <kep3/detail/s11n.hpp>

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

integrator make_variational_integrator()
{
    auto ta = make_integrator(1u, 2u);
    const auto &sys = ta.get_sys();
    auto var_sys = heyoka::var_ode_sys(sys, {sys[0].first, heyoka::par[0]}, 1u);
    return integrator{std::move(var_sys)};
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

TEST_CASE("zoh_ms cached defect gradient sparsity")
{
    // A two-state leg makes the different forward and backward state blocks visible.
    auto leg = make_test_leg_21();
    const auto &[sp_states, sp_controls, sp_tgrid] = leg.defects_grad_sparsity();
    const kep3::leg::zoh_ms::sparsity_pattern expected_states{
        {0u, 0u}, {0u, 1u}, {0u, 2u}, {1u, 0u}, {1u, 1u}, {1u, 3u},
        {2u, 2u}, {2u, 4u}, {2u, 5u}, {3u, 3u}, {3u, 4u}, {3u, 5u},
        {4u, 4u}, {4u, 6u}, {4u, 7u}, {5u, 5u}, {5u, 6u}, {5u, 7u}};
    const kep3::leg::zoh_ms::sparsity_pattern expected_controls{
        {0u, 0u}, {1u, 0u}, {2u, 1u}, {3u, 1u}, {4u, 2u}, {5u, 2u}};
    const kep3::leg::zoh_ms::sparsity_pattern expected_tgrid{
        {0u, 0u}, {0u, 1u}, {1u, 0u}, {1u, 1u}, {2u, 1u}, {2u, 2u},
        {3u, 1u}, {3u, 2u}, {4u, 2u}, {4u, 3u}, {5u, 2u}, {5u, 3u}};
    REQUIRE(sp_states == expected_states);
    REQUIRE(sp_controls == expected_controls);
    REQUIRE(sp_tgrid == expected_tgrid);

    // Value-only setters must leave the structural coordinates unchanged.
    leg.set_states(std::vector<double>(leg.get_states().size(), 1.0));
    leg.set_controls(std::vector<double>(leg.get_controls().size(), 2.0));
    leg.set_tgrid({0.0, 1.0, 2.0, 3.0});
    leg.set_max_steps(10u);
    const auto &[same_states, same_controls, same_tgrid] = leg.defects_grad_sparsity();
    REQUIRE(same_states == expected_states);
    REQUIRE(same_controls == expected_controls);
    REQUIRE(same_tgrid == expected_tgrid);

    // Changing the cut updates which node contributes the full state block.
    leg.set_cut(0.0);
    REQUIRE(std::get<0>(leg.defects_grad_sparsity())[0] == std::pair<std::size_t, std::size_t>{0u, 0u});
    REQUIRE(std::get<0>(leg.defects_grad_sparsity())[1] == std::pair<std::size_t, std::size_t>{0u, 2u});
    leg.set_cut(1.0);
    REQUIRE(std::get<0>(leg.defects_grad_sparsity())[0] == std::pair<std::size_t, std::size_t>{0u, 0u});
    REQUIRE(std::get<0>(leg.defects_grad_sparsity())[3] == std::pair<std::size_t, std::size_t>{1u, 0u});

    // Resizing the mesh rebuilds all three blocks for the new segment count.
    leg.set({0.0, 1.0, 2.0, 3.0, 4.0, 5.0}, {0.1, 0.2}, {0.0, 1.0, 2.0}, 0.5);
    const auto &[resized_states, resized_controls, resized_tgrid] = leg.defects_grad_sparsity();
    REQUIRE(leg.get_nseg() == 2u);
    REQUIRE(resized_states.size() == 12u);
    REQUIRE(resized_controls.size() == 4u);
    REQUIRE(resized_tgrid.size() == 8u);

    // The zero-segment case has valid, empty coordinate patterns.
    const auto ta = make_integrator(2u, 1u);
    kep3::leg::zoh_ms empty_leg{{0.0, 0.0}, {}, {0.0}, 0.5, {ta, std::nullopt}, std::nullopt, 2u, 1u};
    const auto &[empty_states, empty_controls, empty_tgrid] = empty_leg.defects_grad_sparsity();
    REQUIRE(empty_states.empty());
    REQUIRE(empty_controls.empty());
    REQUIRE(empty_tgrid.empty());

    // Loading an archive reconstructs the derived cache from the restored structure.
    std::stringstream archive_data;
    {
        boost::archive::binary_oarchive archive(archive_data);
        archive << leg;
    }
    kep3::leg::zoh_ms restored_leg{};
    {
        boost::archive::binary_iarchive archive(archive_data);
        archive >> restored_leg;
    }
    REQUIRE(restored_leg.defects_grad_sparsity() == leg.defects_grad_sparsity());
}

TEST_CASE("compute_defects") {
    // We test forward and backward defects on a nonuniform-grid 2-state, 1-control leg.
    auto leg = make_test_leg_21();
    const auto defects = leg.compute_defects();
    const std::vector<double> expected{4.487212707001282, 16.901100113049495, 9.653047597773553, 13.203406906114177, 15.10383276937845, 16.952903708643585};
    REQUIRE(kep3_tests::L_infinity_norm_rel(defects, expected) < 1e-12);
}

TEST_CASE("compute_defects_grad")
{
    // We test that sparse gradients match finite differences for states, controls and node times.
    const auto ta = make_integrator(1u, 2u);
    const auto ta_var = make_variational_integrator();
    kep3::leg::zoh_ms leg{{0.2, 0.5, -0.1}, {0.3, -0.2}, {0.0, 0.2, 0.7}, 0.5,
                          {ta, ta_var}, std::nullopt, 1u, 1u};
    const auto [grad_states, grad_controls, grad_tgrid] = leg.compute_defects_grad();
    const auto [sp_states, sp_controls, sp_tgrid] = leg.defects_grad_sparsity();
    const auto states = leg.get_states();
    const auto controls = leg.get_controls();
    const auto tgrid = leg.get_tgrid();
    const double step = 1e-5;

    // Perturb one input at a time and compare the matching sparse entries.
    const auto check_gradient = [&](const std::vector<double> &values, const auto &sparsity,
                                    const std::vector<double> &gradient, const auto &set_values) {
        REQUIRE(gradient.size() == sparsity.size());
        for (std::size_t column = 0u; column < values.size(); ++column) {
            auto plus = values;
            auto minus = values;
            plus[column] += step;
            minus[column] -= step;
            set_values(plus);
            const auto defects_plus = leg.compute_defects();
            set_values(minus);
            const auto defects_minus = leg.compute_defects();

            for (std::size_t entry = 0u; entry < sparsity.size(); ++entry) {
                if (sparsity[entry].second == column) {
                    const auto row = sparsity[entry].first;
                    const auto estimate = (defects_plus[row] - defects_minus[row]) / (2.0 * step);
                    REQUIRE(std::abs(gradient[entry] - estimate) < 2e-6);
                }
            }
        }
        set_values(values);
    };

    check_gradient(states, sp_states, grad_states,
                   [&](const auto &values) { leg.set_states(values); });
    check_gradient(controls, sp_controls, grad_controls,
                   [&](const auto &values) { leg.set_controls(values); });
    check_gradient(tgrid, sp_tgrid, grad_tgrid,
                   [&](const auto &values) { leg.set_tgrid(values); });

    // A nominal-only leg must reject gradient requests.
    auto nominal_only = make_test_leg_11();
    REQUIRE_THROWS_AS(nominal_only.compute_defects_grad(), std::logic_error);
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