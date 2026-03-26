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
#include <kep3/ta/cr3bp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_cr3bp;
using kep3::ta::get_ta_cr3bp_cache_dim;
using kep3::ta::get_ta_cr3bp_var;
using kep3::ta::get_ta_cr3bp_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    // The non variational one.
    REQUIRE(get_ta_cr3bp_cache_dim() == 0u);
    auto ta_cached = get_ta_cr3bp(1e-16);
    REQUIRE(get_ta_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_cr3bp(1e-16);
    REQUIRE(get_ta_cr3bp_cache_dim() == 1u);
    ta_cached = get_ta_cr3bp(1e-8);
    REQUIRE(get_ta_cr3bp_cache_dim() == 2u);

    // The variational integrator.
    REQUIRE(get_ta_cr3bp_var_cache_dim() == 0u);
    ta_cached = get_ta_cr3bp_var(1e-16);
    REQUIRE(get_ta_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_cr3bp_var(1e-16);
    REQUIRE(get_ta_cr3bp_var_cache_dim() == 1u);
    ta_cached = get_ta_cr3bp_var(1e-8);
    REQUIRE(get_ta_cr3bp_var_cache_dim() == 2u);
}

TEST_CASE("propagation")
{
    auto ta_cached = get_ta_cr3bp(1e-16);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 6);
    REQUIRE(ta_cached.get_pars().size() == 1);

    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1.01238082345234, -0.0423523523454,  0.22634376321,
                               -0.1232623614,    0.123462698209365, 0.123667064622};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        // Assign value for cis-lunar case. (note: not the official pykep on in this case)
        *ta.get_pars_data() = 0.01215058560962404;
        auto out = ta.propagate_until(5.7856656782589234);
        std::vector<double> const ground_truth = {0.43038358727124, -1.64650668902846, 0.10271923139472,
                                                  -0.9315629872575, -0.42680151362818, 0.22257221768767};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_propagation")
{
    auto ta_cached = get_ta_cr3bp_var(1e-16);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 42);
    REQUIRE(ta_cached.get_pars().size() == 1);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 6);

    {
        // We test a generic case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        ta.set_time(0.);
        std::vector<double> ic{1.01238082345234, -0.0423523523454,  0.22634376321,
                               -0.1232623614,    0.123462698209365, 0.123667064622};
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        // Assign value for cis-lunar case. (note: not the official pykep on in this case)
        *ta.get_pars_data() = 0.01215058560962404;
        auto out = ta.propagate_until(5.7856656782589234);
        std::vector<double> const ground_truth = {0.43038358727124, -1.64650668902846, 0.10271923139472,
                                                  -0.9315629872575, -0.42680151362818, 0.22257221768767};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(std::vector<double>(ta.get_state().begin(), ta.get_state().begin() + 6),
                                                ground_truth)
                <= 1e-13);
    }
}

TEST_CASE("expression_factory")
{
    // Get the symbolic heyoka expressions.
    auto [x, y, z, vx, vy, vz] = heyoka::make_vars("x", "y", "z", "vx", "vy", "vz");
    auto dyn = kep3::ta::cr3bp_dyn();
    auto C = kep3::ta::cr3bp_jacobi_C();
    auto U = kep3::ta::cr3bp_effective_potential_U();
    // Compile them into functions.
    auto C_cfunc = heyoka::cfunc<double>({C}, {x, y, z, vx, vy, vz});
    auto U_cfunc = heyoka::cfunc<double>({U}, {x, y, z, vx, vy, vz});
    // We build the Taylor integrator from the dynamics.
    auto ta = heyoka::taylor_adaptive<double>{dyn};
    // We initialize it with some initial conditions.
    ta.set_time(0.);
    std::vector<double> ic{1.01238082345234, -0.0423523523454,  0.22634376321,
                           -0.1232623614,    0.123462698209365, 0.123667064622};
    std::copy(ic.begin(), ic.end(), ta.get_state_data());
    // Assign value for cis-lunar case.
    *ta.get_pars_data() = kep3::CR3BP_MU_EARTH_MOON;
    // Propagate to some future time.
    auto out = ta.propagate_until(2.45);
    // Test Jacobi constant being constant along a trajectory
    std::array<double, 1> C0 = {0.}, Cf = {0.}, U0 = {0.};
    C_cfunc(C0, ic, heyoka::kw::pars = {kep3::CR3BP_MU_EARTH_MOON});
    U_cfunc(U0, ic, heyoka::kw::pars = {kep3::CR3BP_MU_EARTH_MOON});
    C_cfunc(Cf, ic, heyoka::kw::pars = {kep3::CR3BP_MU_EARTH_MOON});
    REQUIRE(C0[0] == Approx(Cf[0]).epsilon(1e-14));
    // Test Jacobi constant being defined as C = 2*U - v^2
    REQUIRE(2. * U0[0] - ic[3] * ic[3] - ic[4] * ic[4] - ic[5] * ic[5] == Approx(C0[0]).epsilon(1e-14));
}
