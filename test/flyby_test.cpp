// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cmath>
#include <stdexcept>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <boost/math/constants/constants.hpp>

#include <kep3/core_astro/flyby.hpp>
#include <kep3/udpla/keplerian.hpp>

#include "catch.hpp"
#include "kep3/core_astro/constants.hpp"
#include "test_helpers.hpp"

TEST_CASE("fb_con")
{
    const double mu = 398600441800000.0;
    const double safe_radius = 7015800.000000001;
    // Call from doubles
    {
        const std::array<double, 3> v_rel_in = {7000., 0., 0.};
        const std::array<double, 3> v_rel_out = {7100, 0., 0.};

        const std::pair<double, double> ground_truth{-1410000.0, -1.1335787944199696};
        auto res = kep3::fb_con(v_rel_in, v_rel_out, mu, safe_radius);
        REQUIRE(kep3_tests::floating_point_error(res.first, ground_truth.first) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error(res.second, ground_truth.second) < 1e-13);
    }
    {
        const std::array<double, 3> v_rel_in = {7200., -4567.7655, 1234.4233};
        const std::array<double, 3> v_rel_out = {7100., 220.123, -144.432};
        const std::pair<double, double> ground_truth{23748967.80882015, -0.2801587391727032};
        auto res = kep3::fb_con(v_rel_in, v_rel_out, mu, safe_radius);
        REQUIRE(kep3_tests::floating_point_error(res.first, ground_truth.first) < 1e-14);
        REQUIRE(kep3_tests::floating_point_error(res.second, ground_truth.second) < 1e-14);
    }
    {
        // Call from planet
        // We use 2000-01-01 as a reference epoch for all these tests
        const double ref_epoch = 0.;
        // This is a circular orbit at 1 AU.
        const std::array<std::array<double, 3>, 2> pos_vel_0{{{kep3::AU, 0., 0.}, {0., kep3::EARTH_VELOCITY, 0.}}};
        // A keplerian planet orbiting the Sun on such a perfectly circular orbit.
        const kep3::udpla::keplerian udpla{
            kep3::epoch(ref_epoch), pos_vel_0, kep3::MU_SUN, "test", {mu, safe_radius, safe_radius}};
        const kep3::planet pl_test{udpla};
        const std::array<double, 3> v_rel_in = {7200., -4567.7655, 1234.4233};
        const std::array<double, 3> v_rel_out = {7100., 220.123, -144.432};
        const std::pair<double, double> ground_truth{23748967.80882015, -0.2801587391727032};
        auto res = kep3::fb_con(v_rel_in, v_rel_out, pl_test);
        REQUIRE(kep3_tests::floating_point_error(res.first, ground_truth.first) < 1e-14);
        REQUIRE(kep3_tests::floating_point_error(res.second, ground_truth.second) < 1e-14);
    }
}

TEST_CASE("fb_dv")
{
    const double mu = 398600441800000.0;
    const double safe_radius = 7015800.000000001;
    // Call from doubles
    {
        const std::array<double, 3> v_rel_in = {7000., 0., 0.};
        const std::array<double, 3> v_rel_out = {7100, 0., 0.};
        const double ground_truth = 100.;
        auto res = kep3::fb_dv(v_rel_in, v_rel_out, mu, safe_radius);
        REQUIRE(kep3_tests::floating_point_error(res, ground_truth) < 1e-13);
    }
    {
        const std::array<double, 3> v_rel_in = {7200., -4567.7655, 1234.4233};
        const std::array<double, 3> v_rel_out = {7100., 220.123, -144.432};
        const double ground_truth = 1510.704060449003;
        auto res = kep3::fb_dv(v_rel_in, v_rel_out, mu, safe_radius);
        REQUIRE(kep3_tests::floating_point_error(res, ground_truth) < 1e-14);
    }
    {
        // Call from planet
        // We use 2000-01-01 as a reference epoch for all these tests
        const double ref_epoch = 0.;
        // This is a circular orbit at 1 AU.
        const std::array<std::array<double, 3>, 2> pos_vel_0{{{kep3::AU, 0., 0.}, {0., kep3::EARTH_VELOCITY, 0.}}};
        // A keplerian planet orbiting the Sun on such a perfectly circular orbit.
        const kep3::udpla::keplerian udpla{
            kep3::epoch(ref_epoch), pos_vel_0, kep3::MU_SUN, "test", {mu, safe_radius, safe_radius}};
        const kep3::planet pl_test{udpla};
        const std::array<double, 3> v_rel_in = {7200., -4567.7655, 1234.4233};
        const std::array<double, 3> v_rel_out = {7100., 220.123, -144.432};
        const double ground_truth = 1510.704060449003;
        auto res = kep3::fb_dv(v_rel_in, v_rel_out, mu, safe_radius);
        REQUIRE(kep3_tests::floating_point_error(res, ground_truth) < 1e-14);
    }
}

TEST_CASE("fb_vout")
{
    {
        const double rp = 2.;
        const double beta = 3.1415 / 3.;
        const double mu = 1.;
        const std::array<double, 3> v_in = {1., 0., 0.};
        const std::array<double, 3> v_pla = {0., 1., 0.};
        const std::array<double, 3> ground_truth = {0.5805947972994727, -0.2594052027005272, 0.27714295365321495};
        auto res = kep3::fb_vout(v_in, v_pla, rp, beta, mu);
        REQUIRE(kep3_tests::floating_point_error_vector(ground_truth, res) < 1e-14);
    }
    {
        const double rp = 7600000.;
        const double beta = 3.1415 / 3. * 2.;
        const double mu = kep3::MU_EARTH;
        const std::array<double, 3> v_in = {20000., 1234.5678, -12098.123565234};
        const std::array<double, 3> v_pla = {25000.243234322, 1., 0.};
        const std::array<double, 3> ground_truth = {15773.610960458058, 3862.490616725501, -8535.037364832286};
        auto res = kep3::fb_vout(v_in, v_pla, rp, beta, mu);
        REQUIRE(kep3_tests::floating_point_error_vector(ground_truth, res) < 1e-14);
    }
}