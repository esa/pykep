// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <heyoka/expression.hpp>
#include <heyoka/kw.hpp>
#include <heyoka/taylor.hpp>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/ta/kep.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_kep;
using kep3::ta::get_ta_kep_cache_dim;
using kep3::ta::get_ta_kep_var;
using kep3::ta::get_ta_kep_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_kep_cache_dim() == 0u);
    auto ta_cached = get_ta_kep(1e-16);
    REQUIRE(get_ta_kep_cache_dim() == 1u);
    ta_cached = get_ta_kep(1e-16);
    REQUIRE(get_ta_kep_cache_dim() == 1u);
    ta_cached = get_ta_kep(1e-8);
    REQUIRE(get_ta_kep_cache_dim() == 2u);
    //
    //// The variational integrator.
    REQUIRE(get_ta_kep_var_cache_dim() == 0u);
    ta_cached = get_ta_kep_var(1e-16);
    REQUIRE(get_ta_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_kep_var(1e-16);
    REQUIRE(get_ta_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_kep_var(1e-8);
    REQUIRE(get_ta_kep_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_kep(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 6);
    REQUIRE(ta_cached.get_pars().size() == 1);

    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        const std::array<std::array<double, 3>, 2> pos_vel_s
            = {{{0.43038358727124, -1.64650668902846, 0.10271923139472},
                {-0.9315629872575, -0.42680151362818, 0.22257221768767}}};
        std::ranges::copy(pos_vel_s | std::views::join, ta.get_state_data());
        // Assign value for cis-lunar case. (note: not the official pykep on in this case)
        double tof = 5.7856656782589234;
        double mu = 1.215058560962404;
        *ta.get_pars_data() = mu;
        auto out = ta.propagate_until(tof);
        // We do the same but with the propagate lagrangian function
        auto const propagated_posvel = kep3::propagate_lagrangian(pos_vel_s, tof, mu);
        std::vector<double> res_pl(6);
        std::ranges::copy(propagated_posvel.first | std::views::join, res_pl.begin());
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), res_pl) <= 1e-13);
    }
}

TEST_CASE("variational_propagation")
{
    auto ta_cached = get_ta_kep_var(1e-16);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 42);
    REQUIRE(ta_cached.get_pars().size() == 1);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 6);

    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        const std::array<std::array<double, 3>, 2> pos_vel_s
            = {{{0.43038358727124, -1.64650668902846, 0.10271923139472},
                {-0.9315629872575, -0.42680151362818, 0.22257221768767}}};
        std::ranges::copy(pos_vel_s | std::views::join, ta.get_state_data());
        // Assign value for cis-lunar case. (note: not the official pykep on in this case)
        double tof = 5.7856656782589234;
        double mu = 1.215058560962404;
        *ta.get_pars_data() = mu;
        auto out = ta.propagate_until(tof);
        // We do the same but with the propagate lagrangian function
        auto const propagated_posvel = kep3::propagate_lagrangian(pos_vel_s, tof, mu, true);
        std::vector<double> posvel_pl(6);
        std::ranges::copy(propagated_posvel.first | std::views::join, posvel_pl.begin());
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(std::vector<double>(ta.get_state().begin(), ta.get_state().begin() + 6),
                                                posvel_pl)
                <= 1e-13);
        // We now check the STM
        std::vector<double> stm_pl(propagated_posvel.second.value().begin(), propagated_posvel.second.value().end());
        REQUIRE(kep3_tests::L_infinity_norm_rel(std::vector<double>(ta.get_state().begin()+6, ta.get_state().end()),
                                                stm_pl)
                <= 1e-13);
    }
}