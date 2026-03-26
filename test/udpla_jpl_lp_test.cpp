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
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/jpl_lp.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using kep3::udpla::jpl_lp;

TEST_CASE("constructor")
{
    REQUIRE_NOTHROW(jpl_lp{});
    kep3::epoch ref_epoch{12.22, kep3::epoch::julian_type::MJD2000};
    // From name
    REQUIRE_NOTHROW(jpl_lp{"Mars"});
    REQUIRE_NOTHROW(jpl_lp{"mars"});
    REQUIRE_NOTHROW(kep3::planet{jpl_lp{"neptune"}});

    REQUIRE_THROWS_AS(jpl_lp{"gigi"}, std::logic_error);
    jpl_lp udpla{"earth"};
    kep3::planet earth{udpla};
    REQUIRE(kep3_tests::floating_point_error(earth.period() * kep3::SEC2DAY, 365.25) < 0.01);
}

TEST_CASE("eph")
{
    // We use 2020-01-01 as a reference epoch for all these tests. To compute the ground truth 
    // we queried JPL Horizon. Since the queries are made at different times the eph used are also
    // different. As a consequence if you repeat the query today you may get a different result.
    // Note also the low error tolerance requested as low precision eph are indeed low precision.
    double ref_epoch = kep3::epoch(2458849.5, kep3::epoch::julian_type::JD).mjd2000();
    {
        // This is Mercury w.r.t. the Sun queried from JPL Horizon (DE441) at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{-1.004313454665499E+10, -6.782852744259514E+10, -4.760875668975301E+09},
             {3.847265152100766E+04, -4.158689617302953E+03, -3.869763814829368E+03}}};
        // Mercury in jpl_lp mode
        jpl_lp udpla{"mercury"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.05);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.05);
    }
    {
        // This is Venus w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{1.081892249749067E+11, 7.861125522626230E+09, -6.135421905733132E+09},
             {-2.679023504807336E+03, 3.476995213020635E+04, 6.316923820826013E+02}}};
        // Venus in jpl_lp mode
        jpl_lp udpla{"venus"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.02);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.02);
    }
    {
        // This is the Earth-Moon w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{-2.488023054631234E+10, 1.449771522542222E+11, -6.590293144971132E+02},
             {-2.984589828430694E+04, -5.151004951052294E+03, 3.108878527788850E-01}}};
        // The Earth-Moon in jpl_lp mode
        jpl_lp udpla{"earth"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
    }

    {
        // This is Mars w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{-1.974852868472516E+11, -1.325074305699912E+11, 2.068800235454373E+09},
             {1.440720082303430E+04, -1.804659323991406E+04, -7.316474757792575E+02}}};
        // Mars in jpl_lp mode
        jpl_lp udpla{"mars"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
    }
    {
        // This is Jupiter w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{7.871048884706007E+10, -7.780620023532844E+11, 1.470618693758428E+09},
             {1.285491315086331E+04, 1.933721291664580E+03, -2.956500740610059E+02}}};
        // Jupiter in jpl_lp mode
        jpl_lp udpla{"jupiter"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
    }
    {
        // This is Saturn w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{5.680597453102431E+11, -1.389479460523918E+12, 1.545819892540634E+09},
             {8.420955066542843E+03, 3.631222339233865E+03, -3.987639953503348E+02}}};
        // Uranus in jpl_lp mode
        jpl_lp udpla{"sAtURN"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
    }
    {
        // This is Uranus w.r.t. the Sun queried from JPL Horizon at
        // 2020-01-01
        std::array<std::array<double, 3>, 2> pos_vel_0{
            {{2.427299831783689E+12, 1.702279661193254E+12, -2.511573638659978E+09},
             {-3.947889815690411E+03, 5.259970919483185E+03, 7.055304073010626E+01}}};
        // Uranus in jpl_lp mode
        jpl_lp udpla{"uranus"};
        auto [r, v] = udpla.eph(ref_epoch);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel_0[0]) < 0.01);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel_0[1]) < 0.01);
    }
    jpl_lp udpla{"uranus"};
    REQUIRE_THROWS_AS(udpla.eph(5347534.), std::domain_error);
}

TEST_CASE("elements")
{
    double ref_epoch = 12.22;
    // We use Neptune
    jpl_lp udpla{"nePTUne"}; // casing is not important
    auto pos_vel = udpla.eph(ref_epoch);
    // Test on various element types
    {
        auto par = udpla.elements(ref_epoch, kep3::elements_type::KEP_F);
        auto [r, v] = kep3::par2ic(par, kep3::MU_SUN);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = udpla.elements(ref_epoch, kep3::elements_type::KEP_M);
        par[5] = kep3::m2f(par[5], par[1]);
        auto [r, v] = kep3::par2ic(par, kep3::MU_SUN);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = udpla.elements(ref_epoch, kep3::elements_type::MEE);
        auto [r, v] = kep3::mee2ic(par, kep3::MU_SUN);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        auto par = udpla.elements(ref_epoch, kep3::elements_type::MEE_R);
        auto [r, v] = kep3::mee2ic(par, kep3::MU_SUN, true);
        REQUIRE(kep3_tests::floating_point_error_vector(r, pos_vel[0]) < 1e-13);
        REQUIRE(kep3_tests::floating_point_error_vector(v, pos_vel[1]) < 1e-13);
    }
    {
        REQUIRE_THROWS_AS(udpla.elements(ref_epoch, kep3::elements_type::POSVEL), std::logic_error);
    }
}

TEST_CASE("getters-setters")
{
    jpl_lp udpla{"nePTUne"}; // casing is not important
    REQUIRE(udpla.get_name() == "neptune(jpl_lp)");
    REQUIRE(udpla.get_mu_central_body() == kep3::MU_SUN);
    REQUIRE(udpla.get_mu_self() == 6836529e9);
    REQUIRE(udpla.get_radius() == 24622000.);
    REQUIRE(udpla.get_safe_radius() == 1.1 * 24622000.);
    udpla.set_safe_radius(1.1 * 24622000 + 1000);
    REQUIRE(udpla.get_safe_radius() == 1.1 * 24622000 + 1000);
    REQUIRE_THROWS_AS(udpla.set_safe_radius(12340.1234), std::domain_error);
}

TEST_CASE("stream_operator")
{
    REQUIRE_NOTHROW((std::cout << jpl_lp{} << '\n'));
}

TEST_CASE("serialization_test")
{
    // Instantiate a generic jpl_lp
    kep3::epoch ref_epoch{2423.4343, kep3::epoch::julian_type::MJD2000};
    jpl_lp udpla{"neptune"};

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
    jpl_lp udpla2{};
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
    kep3::planet pla{kep3::udpla::jpl_lp{"earth"}};

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