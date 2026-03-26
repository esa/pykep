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
#include <kep3/ta/cr3bp.hpp>
#include <kep3/ta/zero_hold_cr3bp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_cr3bp;
using kep3::ta::get_ta_cr3bp_var;
using kep3::ta::get_ta_zero_hold_cr3bp;
using kep3::ta::get_ta_zero_hold_cr3bp_cache_dim;
using kep3::ta::get_ta_zero_hold_cr3bp_var;
using kep3::ta::get_ta_zero_hold_cr3bp_var_cache_dim;

using kep3_tests::L_infinity_norm_rel;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_zero_hold_cr3bp_cache_dim() == 0u);
    auto ta_cached = get_ta_zero_hold_cr3bp(1e-16);
    REQUIRE(get_ta_zero_hold_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_cr3bp(1e-16);
    REQUIRE(get_ta_zero_hold_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_cr3bp(1e-8);
    REQUIRE(get_ta_zero_hold_cr3bp_cache_dim() == 2u);

    // The variational integrator.
    REQUIRE(get_ta_zero_hold_cr3bp_var_cache_dim() == 0u);
    ta_cached = get_ta_zero_hold_cr3bp_var(1e-16);
    REQUIRE(get_ta_zero_hold_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_cr3bp_var(1e-16);
    REQUIRE(get_ta_zero_hold_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_zero_hold_cr3bp_var(1e-8);
    REQUIRE(get_ta_zero_hold_cr3bp_var_cache_dim() == 2u);
}

TEST_CASE("dynamics")
{
    auto ta_cached = get_ta_zero_hold_cr3bp(1e-16);
    auto ta_cached_cr3bp = get_ta_cr3bp(1e-16);

    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 7);
    REQUIRE(ta_cached.get_pars().size() == 5);

    {
        // We test a ballistic circular orbit for mu=1 and one period.
        taylor_adaptive<double> ta(ta_cached);             // making a copy as to be able to modify the object.
        taylor_adaptive<double> ta_cr3bp(ta_cached_cr3bp); // making a copy as to be able to modify the object.

        std::vector<double> ic{1., 0.1, 0.2, 0.3, 1.1, 0.1, 1.23};
        std::vector<double> pars{1.234, 1.234, 0., 0., 0.};

        // We set the zero_hold version
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        // We set the cr3bp version
        ta_cr3bp.set_time(0.);
        std::copy(ic.begin(), ic.end() - 1l, ta_cr3bp.get_state_data());
        *ta_cr3bp.get_pars_data() = pars[0]; // mu

        // We call propagation
        auto out_zero_hold = ta.propagate_until(2. * kep3::pi);
        auto out_cr3bp = ta_cr3bp.propagate_until(2. * kep3::pi);

        REQUIRE(std::get<0>(out_zero_hold) == taylor_outcome::time_limit);
        REQUIRE(std::get<0>(out_cr3bp) == taylor_outcome::time_limit);

        // Remove the mass
        std::vector<double> state_no_mass(ta.get_state().begin(), ta.get_state().end() - 1l);

        // Check states are the same (up to the mass which is not there in the cr3bp case)
        REQUIRE(L_infinity_norm_rel(state_no_mass, ta_cr3bp.get_state()) <= 1e-13);

        // Check that the mass is unchanged
        REQUIRE(std::abs(ta.get_state()[6] - ic[6]) <= 1e-10);
    }
    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1., 0.1, 0.2, 0.3, 1.1, 0.1, 1.23};
        ;
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{1.234, 1.234, 0.0123, -0.134, 0.0122};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(2.34);
        std::vector<double> const ground_truth
            = {2.4229189094753427, -3.747329291698182,  0.3346665501262248, -2.5629530910977696,
               -3.655884167388002, 0.05495744837535643, 0.9737847016016398};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics")
{
    auto ta_cached = get_ta_zero_hold_cr3bp_var(1e-16);
    auto ta_cached_cr3bp = get_ta_cr3bp_var(1e-16);

    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 77);
    REQUIRE(ta_cached.get_pars().size() == 5);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 10);

    {
        // We test a ballistic circular orbit for mu=1 and one period.
        taylor_adaptive<double> ta(ta_cached);             // making a copy as to be able to modify the object.
        taylor_adaptive<double> ta_cr3bp(ta_cached_cr3bp); // making a copy as to be able to modify the object.

        std::vector<double> ic{1., 0.1, 0.2, 0.3, 1.1, 0.1, 1.23};
        std::vector<double> pars{1.234, 1.234, 0., 0., 1e-22};

        // We set the zero_hold version
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        // We set the cr3bp version
        ta_cr3bp.set_time(0.);
        std::copy(ic.begin(), ic.end() - 1l, ta_cr3bp.get_state_data());
        *ta_cr3bp.get_pars_data() = pars[0]; // mu

        // We call propagation
        auto out_zero_hold = ta.propagate_until(2. * kep3::pi);
        auto out_cr3bp = ta_cr3bp.propagate_until(2. * kep3::pi);
        REQUIRE(std::get<0>(out_zero_hold) == taylor_outcome::time_limit);
        REQUIRE(std::get<0>(out_cr3bp) == taylor_outcome::time_limit);

        // Remove the mass
        std::vector<double> state_no_mass(ta.get_state().begin(), ta.get_state().begin() + 6l);

        // Check states are the same (up to the mass which is not there in the cr3bp case)
        REQUIRE(L_infinity_norm_rel(state_no_mass, ta_cr3bp.get_state()) <= 1e-13);

        // Check that the mass is unchanged (theres a small thrust effect that avoids singularities in variational eqs)
        REQUIRE(std::abs(ta.get_state()[6] - ic[6]) <= 1e-13);
    }
    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1., 0.1, 0.2, 0.3, 1.1, 0.1, 1.23};
        ;
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::vector<double> pars{1.234, 1.234, 0.0123, -0.134, 0.0122};
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());
        auto out = ta.propagate_until(2.34);
        // clang-format off
        std::vector<double> const ground_truth
            = {2.4229189094753414, -3.747329291698182, 0.33466655012622487, -2.562953091097771, -3.6558841673880016, 0.05495744837535647, 0.9737847016016395, 0.9614045765001359, 1.7825151782583069, -0.05380751373522139, -1.4916942841350285, 1.6060590962977148, 0.0024713194199781135, 0.23498459426780896, 0.0914800320274582, 2.1645517884972714, -0.008714400385924139, -3.528505676262339, 0.8303241271221435, -0.23770254041640587, -1.9169791848431046, -1.8967323269032756, -0.07633496547661953, 0.028075797602765806, -2.1249597040958346, 0.060224039780949784, -0.02673814429310827, 0.2643581134550415, 0.05807472956155947, 0.49713225005173545, 0.05380237317326432, 0.052367695746998164, 2.153245008056082, -0.023992310361915153, 0.01661325075820763, -0.0051764833633344155, 2.3452909827985953, -2.729442109347333, 1.225646733237797, -0.25367743960762396, -2.4916597576377355, -1.1888115995264377, -0.06881197651631461, 0.1979604842362697, -1.3440572093079388, 1.6907243082970393, -0.0329864446164595, -2.3392350042527363, -1.1874243102825033, -0.07083004120954858, 0.6411399930116608, -2.5168292443865474, -0.052113020181042066, -0.12697424228277854, -1.6848782320089704, -1.322401263483442, -0.024545003709288412, 0.14140398386331368, 0.03111548201964529, -0.2633220832330632, 0.0311564231196948, 0.03732800082256391, 0.8671386699376833, -0.022244893804737154, 0.015004028067529112, -0.013356525038244298, 2.0808930720878154, 0, 0, 0, 0, 0, 0, 1, -0.17262410874045223, 1.8806203716439507, -0.17122066070191191};
        // clang-format on
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}
