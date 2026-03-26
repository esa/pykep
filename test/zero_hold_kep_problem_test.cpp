// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <iostream>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/zero_hold_kep_problem.hpp>

#include "catch.hpp"
#include "kep3/core_astro/constants.hpp"
#include "test_helpers.hpp"

using kep3_tests::L_infinity_norm_rel;

TEST_CASE("construct")
{
    // Here we test construction for a simple geometry
    REQUIRE_NOTHROW(kep3::zero_hold_kep_problem(1., 10., 1e-12));
    auto s = kep3::zero_hold_kep_problem(1.23, 10.12, 1e-12);
    REQUIRE(s.get_ta().get_pars()[0] == 1.23);
    REQUIRE(s.get_ta().get_pars()[1] == 10.12);
    REQUIRE(s.get_ta_var().get_pars()[0] == 1.23);
    REQUIRE(s.get_ta_var().get_pars()[1] == 10.12);
}

TEST_CASE("propagate")
{
    // The ground truth were computed with these values of the astro constants, hence we cannot use here the pykep ones.
    double AU_OLD = 149597870700.0;
    double EV_OLD = 29784.691831696804;
    double MU_OLD = 1.32712440018e20;

    auto s = kep3::zero_hold_kep_problem(MU_OLD, 3000 * kep3::G0, 1e-16);
    // Mass is zero
    REQUIRE_THROWS_AS(s.propagate({AU_OLD, 1000, 01000, 0.123, EV_OLD, 10.1234, 0.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC), std::domain_error);
    // Radius is zero
    REQUIRE_THROWS_AS(s.propagate({0., 0., 0., 0.123, EV_OLD, 10.1234, 0.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC), std::domain_error);
    // Thrust is zero
    REQUIRE_NOTHROW(s.propagate({AU_OLD, 1000, 01000, 0.123, EV_OLD, 10.1234, 1000.}, {0., 0., 0.}, 300. * kep3::DAY2SEC));

    auto res = s.propagate({AU_OLD, 1000, 1000, 0.123, EV_OLD, 10.1234, 1000.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    std::array<double, 7> const ground_truth{86447318382.73862, -109196478678.49327, 9774837.693084948, 25535.89383875655, 18933.302734531168, -47.00078590023982, 954.2200884785473};
    REQUIRE(L_infinity_norm_rel(res, ground_truth) <=1e-13);
}
TEST_CASE("propagate_var")
{
    // The ground truth were computed with these values of the astro constants, hence we cannot use here the pykep ones.
    double AU_OLD = 149597870700.0;
    double EV_OLD = 29784.691831696804;
    double MU_OLD = 1.32712440018e20;

    auto s = kep3::zero_hold_kep_problem(MU_OLD, 3000 * kep3::G0, 1e-16);
    // Mass is zero
    REQUIRE_THROWS_AS(s.propagate_var({AU_OLD, 1000, 01000, 0.123, EV_OLD, 10.1234, 0.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC), std::domain_error);
    // Radius is zero
    REQUIRE_THROWS_AS(s.propagate_var({0., 0., 0., 0.123, EV_OLD, 10.1234, 0.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC), std::domain_error);
    // Thrust is zero
    REQUIRE_THROWS_AS(s.propagate_var({AU_OLD, 1000, 01000, 0.123, EV_OLD, 10.1234, 1000.}, {0., 0., 0.}, 300. * kep3::DAY2SEC), std::domain_error);

    auto res = s.propagate_var({AU_OLD, 1000, 1000, 0.123, EV_OLD, 10.1234, 1000.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    std::array<double, 7> const ground_truth{86447318382.73862, -109196478678.49327, 9774837.693084948, 25535.89383875655, 18933.302734531168, -47.00078590023982, 954.2200884785473};
    REQUIRE(L_infinity_norm_rel(std::get<0>(res), ground_truth) <=1e-13);
}
TEST_CASE("getters_setters")
{
    auto s = kep3::zero_hold_kep_problem(kep3::MU_SUN, 3000 * kep3::G0, 1e-16);
    REQUIRE(s.get_mu() == kep3::MU_SUN);
    REQUIRE(s.get_veff() == 3000 * kep3::G0);
    REQUIRE(s.get_tol() == 1e-16);
    s.set_mu(12.34);
    REQUIRE(s.get_mu() == 12.34);
    REQUIRE(s.get_ta().get_pars()[0] == 12.34);
    REQUIRE(s.get_ta_var().get_pars()[0] == 12.34);
    s.set_veff(1234);
    REQUIRE(s.get_veff() == 1234);
    REQUIRE(s.get_ta().get_pars()[1] == 1234);
    REQUIRE(s.get_ta_var().get_pars()[1] == 1234);
}

TEST_CASE("serialization_test")
{
    // Instantiate a planet
    kep3::zero_hold_kep_problem s1{kep3::MU_SUN, 3000 * kep3::G0, 1e-16};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(s1);
    // Now serialize, deserialize and compare the result.
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << s1;
    }
    // Create a new planet object
    auto s2 = kep3::zero_hold_kep_problem{};
    boost::lexical_cast<std::string>(s2); // triggers the streaming operator
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> s2;
    }
    auto after = boost::lexical_cast<std::string>(s2);
    REQUIRE(before == after);
    auto res1 = s1.propagate({kep3::AU, 21000, 1221000, 0.321, kep3::EARTH_VELOCITY, 10.1234, 1300.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    auto res2 = s2.propagate({kep3::AU, 21000, 1221000, 0.321, kep3::EARTH_VELOCITY, 10.1234, 1300.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    REQUIRE(L_infinity_norm_rel(res1, res2) < 1e-14);
    auto res1_var = s1.propagate_var({kep3::AU, 21000, 1221000, 0.321, kep3::EARTH_VELOCITY, 10.1234, 1300.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    auto res2_var = s2.propagate_var({kep3::AU, 21000, 1221000, 0.321, kep3::EARTH_VELOCITY, 10.1234, 1300.}, {0.05, 0.01, 0.01}, 300. * kep3::DAY2SEC);
    REQUIRE(L_infinity_norm_rel(std::get<0>(res1_var), std::get<0>(res2_var)) < 1e-14);
    REQUIRE(L_infinity_norm_rel(std::get<1>(res1_var), std::get<1>(res2_var)) < 1e-14);
    REQUIRE(L_infinity_norm_rel(std::get<2>(res1_var), std::get<2>(res2_var)) < 1e-14);
}