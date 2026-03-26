// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <stdexcept>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/keplerian.hpp>

#include "catch.hpp"

#include "test_helpers.hpp"

using kep3::epoch;
using kep3::udpla::keplerian;

TEST_CASE("constructor")
{
    REQUIRE_NOTHROW(keplerian{});
    kep3::epoch ref_epoch{12.22, kep3::epoch::julian_type::MJD2000};
    // From posvel
    REQUIRE_NOTHROW(keplerian{ref_epoch, {{{0.3, 1., 0.2}, {0.0, 1.12, 0.}}}, 1., "unknown"});
    REQUIRE_NOTHROW(keplerian{ref_epoch, {{{0.3, 1., 0.2}, {0.0, 1.12, 0.}}}, 1., "unknown", {-1, -1, -1}});
    // From parameters kep3::elements_type::KEP_F
    std::array<double, 6> par0{{1., 0., 0., 0., 0., 0.}};
    REQUIRE_NOTHROW(keplerian{ref_epoch, par0, 1., "unknown"});
    REQUIRE_NOTHROW(keplerian{ref_epoch, par0, 1., "unknown", {-1, -1, -1}});
    // Checking the data members initializations:
    keplerian udpla{ref_epoch, par0, 1.1, "unknown", {1.2, 2.2, 1.9}};
    REQUIRE(udpla.get_ref_epoch() == ref_epoch);
    REQUIRE(udpla.get_name() == "unknown");
    REQUIRE(udpla.get_mu_central_body() == 1.1);
    REQUIRE(udpla.get_mu_self() == 1.2);
    REQUIRE(udpla.get_radius() == 2.2);
    REQUIRE(udpla.get_safe_radius() == 1.9);
    REQUIRE(udpla.period() == 2 * kep3::pi * std::sqrt(1. / 1.1));
    // Calling constructor with different elements type
    {
        std::array<double, 6> par{{1., 0., 0., 0., 0., 0.}};
        REQUIRE_NOTHROW(keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::KEP_F});
    }
    {
        std::array<double, 6> par{{1., 0., 0., 0., 0., 0.}};
        REQUIRE_NOTHROW(keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::KEP_M});
    }
    {
        std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
        REQUIRE_NOTHROW(keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::MEE});
    }
    {
        std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
        REQUIRE_NOTHROW(keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::MEE_R});
    }
    { // hyperbola and mean anomaly????
        std::array<double, 6> par{{-10., 10., 0., 1., 0., 0.}};
        REQUIRE_THROWS_AS((keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::KEP_M}),
                          std::logic_error);
    }
    { // posvel as 1x6 orbital parameters????
        std::array<double, 6> par{{1., 0., 0., 1., 0., 0.}};
        REQUIRE_THROWS_AS((keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::POSVEL}),
                          std::logic_error);
    }
    { // negative a but ecc < 1????
        std::array<double, 6> par{{-10., 0., 0., 1., 0., 0.}};
        REQUIRE_THROWS_AS((keplerian{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::KEP_F}),
                          std::domain_error);
    }
    { // We construct an hyperbolic planet
        std::array<double, 6> par{{-10., 10., 0., 0., 0., 0.}};
        keplerian udpla2{ref_epoch, par, 1., "unknown", {-1, -1, -1}, kep3::elements_type::KEP_F};
        REQUIRE(!std::isfinite(udpla2.period()));
    }
    { // We construct an hyperbolic planet
        std::array<std::array<double, 3>, 2> posvel{{{1, 0, 0}, {0, 10, 0}}};
        keplerian udpla2{ref_epoch, posvel, 1., "unknown", {-1, -1, -1}};
        REQUIRE(!std::isfinite(udpla2.period()));
    }
}

TEST_CASE("eph")
{
    // We use 2000-01-01 as a reference epoch for all these tests
    double ref_epoch = 0.;
    // This is a circular orbit at 1 AU.
    std::array<std::array<double, 3>, 2> pos_vel_0{{{kep3::AU, 0., 0.}, {0., kep3::EARTH_VELOCITY, 0.}}};
    // A keplerian planet orbiting the Sun on such a perfectly circular orbit.
    keplerian udpla{kep3::epoch(ref_epoch), pos_vel_0, kep3::MU_SUN};
    double period_in_days = (2. * kep3::pi * kep3::AU / kep3::EARTH_VELOCITY) * kep3::SEC2DAY;
    auto [r, v] = udpla.eph(period_in_days);

    // Testing after a period - same eph
    REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 1e-13);

    // Testing ephs are different after 3 months.
    auto pos_vel = udpla.eph(period_in_days + 90.);
    fmt::print("{}\n{}", pos_vel, pos_vel_0);
    REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[0], pos_vel_0[0]) > 1e-4);
    REQUIRE(kep3_tests::floating_point_error_vector(pos_vel[1], pos_vel_0[1]) > 1e-4);
}

TEST_CASE("elements")
{
    kep3::epoch ref_epoch{12.22, kep3::epoch::julian_type::MJD2000};
    // Non singular elements
    std::array<std::array<double, 3>, 2> pos_vel{{{1., 0.1, 0.1}, {0.1, 1., 0.1}}};
    keplerian udpla{ref_epoch, pos_vel};
    kep3::planet pla{udpla};
    // Test on various element types
    {
        auto par = pla.elements(ref_epoch.mjd2000(), kep3::elements_type::KEP_F);
        auto [r, v] = kep3::par2ic(par, 1.);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = pla.elements(ref_epoch.mjd2000(), kep3::elements_type::KEP_M);
        par[5] = kep3::m2f(par[5], par[1]);
        auto [r, v] = kep3::par2ic(par, 1.);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = pla.elements(ref_epoch.mjd2000(), kep3::elements_type::MEE);
        auto [r, v] = kep3::mee2ic(par, 1.);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = pla.elements(ref_epoch.mjd2000(), kep3::elements_type::MEE_R);
        auto [r, v] = kep3::mee2ic(par, 1., true);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        REQUIRE_THROWS_AS(pla.elements(ref_epoch.mjd2000(), kep3::elements_type::POSVEL), std::logic_error);
    }
    // We test on hyperbolas
    std::array<std::array<double, 3>, 2> pos_vel_h{{{1., 0.1, 0.1}, {0.1, 3., 0.1}}};
    keplerian udpla_h{kep3::epoch(ref_epoch), pos_vel_h, kep3::MU_SUN};

    {
        REQUIRE_THROWS_AS(pla.elements(ref_epoch.mjd2000(), kep3::elements_type::POSVEL), std::logic_error);
    }
}

TEST_CASE("stream_operator")
{
    REQUIRE_NOTHROW((std::cout << keplerian{} << '\n'));
}

TEST_CASE("serialization_test")
{
    // Instantiate a generic udpla
    kep3::epoch ref_epoch{2423.4343, kep3::epoch::julian_type::MJD2000};
    keplerian udpla{ref_epoch, {{{0.33, 1.3, 0.12}, {0.01, 1.123, 0.2}}}, 1.12, "enterprise", {12.32, 44.6, 98.23}};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(udpla);
    // Now serialize
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << udpla;
    }

    // Deserialize
    // Create a new udpla object
    keplerian udpla2{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> udpla2;
    }
    auto after = boost::lexical_cast<std::string>(udpla2);
    // Compare the string represetation
    REQUIRE(before == after);
}

TEST_CASE("serialization_test_2")
{
    // Instantiate a planet with jpl_lp udpla
    kep3::epoch ref_epoch{2423.4343, kep3::epoch::julian_type::MJD2000};
    kep3::udpla::keplerian udpla{
        ref_epoch, {{{0.33, 1.3, 0.12}, {0.01, 1.123, 0.2}}}, 1.12, "enterprise", {12.32, 44.6, 98.23}};
    kep3::planet pla{udpla};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(pla);
    // Now serialize, deserialize and compare the result.
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << pla;
    }
    // Create a new planet object
    auto pla2 = kep3::planet{kep3::detail::null_udpla{}};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> pla2;
    }
    auto after = boost::lexical_cast<std::string>(pla2);
    REQUIRE(before == after);
    // Check explicitly that the properties of base_p where restored as well.
    REQUIRE(pla.get_mu_central_body() == pla2.get_mu_central_body());
    REQUIRE(pla.get_mu_self() == pla2.get_mu_self());
    REQUIRE(pla.get_radius() == pla2.get_radius());
    REQUIRE(pla.get_safe_radius() == pla2.get_safe_radius());
}