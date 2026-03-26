// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <heyoka/taylor.hpp>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/zero_hold_kep.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_zero_hold_kep;
using kep3::ta::get_ta_zero_hold_kep_cache_dim;
using kep3::ta::get_ta_zero_hold_kep_var;
using kep3::ta::get_ta_zero_hold_kep_var_cache_dim;

using kep3_tests::L_infinity_norm_rel;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_zero_hold_kep_cache_dim() == 0u);
    auto ta_cached = get_ta_zero_hold_kep(1e-16);
    REQUIRE(get_ta_zero_hold_kep_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_kep(1e-16);
    REQUIRE(get_ta_zero_hold_kep_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_kep(1e-8);
    REQUIRE(get_ta_zero_hold_kep_cache_dim() == 2u);

    // The variational integrator.
    REQUIRE(get_ta_zero_hold_kep_var_cache_dim() == 0u);
    ta_cached = get_ta_zero_hold_kep_var(1e-16);
    REQUIRE(get_ta_zero_hold_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_kep_var(1e-16);
    REQUIRE(get_ta_zero_hold_kep_var_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_kep_var(1e-8);
    REQUIRE(get_ta_zero_hold_kep_var_cache_dim() == 2u);
}

TEST_CASE("dynamics")
{
    auto ta_cached = get_ta_zero_hold_kep(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 7);
    REQUIRE(ta_cached.get_pars().size() == 5);

    {
        // We test a ballistic circular orbit for mu=1 and one period.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1., 0., 0., 0., 1., 0., 1.};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{1., 1., 0., 0., 0.};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(2. * kep3::pi);
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(L_infinity_norm_rel(ta.get_state(), ic) <= 1e-13);
    }
    {
        // We test a generic case.
        // The ground truth were computed with these values of the astro constants, hence we cannot use here the pykep
        // ones.
        double MU_OLD = 1.32712440018e20;
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{164557657760.1, 179517444829.19998, 47871318621.12, 32763.159550000004,
                               23827.7524,     9531.10096,         1200.};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{MU_OLD, 3000. * kep3::G0, 0.05, 0., 0.032};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(3888000.0);
        std::vector<double> const ground_truth
            = {284296823432.0578,  263961690798.0665, 82814214381.94377, 29341.50292902231,
               20219.294034700008, 8592.028822618351, 1192.1548315009777};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics")
{
    auto ta_cached = get_ta_zero_hold_kep_var(1e-16);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 77);
    REQUIRE(ta_cached.get_pars().size() == 5);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 10);

    {
        // We test a ballistic circular orbit for mu=1 and one period.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1., 0., 0., 0., 1., 0., 1.};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{1., kep3::G0, 0., 0., 1e-14};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(2. * kep3::pi);
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(L_infinity_norm_rel(std::vector<double>(ta.get_state().begin(), ta.get_state().begin() + 7), ic)
                <= 1e-13);
    }
    {
        // We test a generic case.
        // The ground truth were computed with these values of the astro constants, hence we cannot use here the pykep
        // ones.
        double MU_OLD = 1.32712440018e20;
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{164557657760.1, 179517444829.19998, 47871318621.12, 32763.159550000004,
                               23827.7524,     9531.10096,         1200.};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{MU_OLD, 3000. * kep3::G0, 0.05, 0., 0.032};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(3888000.0);
        std::vector<double> const ground_truth
            = {284296823432.0578,  263961690798.0665, 82814214381.94377, 29341.50292902231,
               20219.294034700008, 8592.028822618351, 1192.1548315009777};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(std::vector<double>(ta.get_state().begin(), ta.get_state().begin() + 7),
                                                ground_truth)
                <= 1e-13);
    }
}
