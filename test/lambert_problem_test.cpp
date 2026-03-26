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

#include <kep3/lambert_problem.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

TEST_CASE("construct")
{
    // Here we test construction for a simple geometry
    REQUIRE_NOTHROW(kep3::lambert_problem{{1., 0., 0.}, {0., 1., 0.}, 3 * kep3::pi / 2, 1., true, 100});
    // And we test the throws
    REQUIRE_THROWS_AS((kep3::lambert_problem{{1., 0., 0.}, {0., 1., 0.}, 3 * kep3::pi / 2, -1.2, true, 100}),
                      std::domain_error);
    REQUIRE_THROWS_AS((kep3::lambert_problem{{1., 0., 0.}, {0., 1., 0.}, -3 * kep3::pi / 2, 1.2, true, 100}),
                      std::domain_error);
    REQUIRE_THROWS_AS((kep3::lambert_problem{{0, 0., 1.}, {0., 1., 0.}, 3 * kep3::pi / 2, 1.2, true, 100}),
                      std::domain_error);
}

TEST_CASE("delta_guidance")
{
    // Here we test that in a number of randomly generated Lambert Problems
    // The boundary data must satisfy the Delta Guidance law

    // Preamble
    std::array<double, 3> r0{{0, 0, 0}}, r1{{0, 0, 0}};
    double tof = 0.;

    // NOLINTNEXTLINE(cert-msc32-c, cert-msc51-cpp)
    std::mt19937 rng_engine(12201203u);
    std::uniform_int_distribution<unsigned> cw_d(0, 1);
    std::uniform_real_distribution<double> r_d(-2, 2);
    std::uniform_real_distribution<double> tof_d(2., 40.);
    std::uniform_real_distribution<double> mu_d(0.9, 1.1);
    unsigned revs_max = 20u;

    unsigned trials = 10000u;

    for (auto i = 0u; i < trials; ++i) {
        // 1 - generate a random problem geometry
        r0[0] = r_d(rng_engine);
        r0[1] = r_d(rng_engine);
        r0[2] = r_d(rng_engine);
        r1[0] = r_d(rng_engine);
        r1[1] = r_d(rng_engine);
        r1[2] = r_d(rng_engine);
        tof = tof_d(rng_engine);
        bool cw = static_cast<bool>(cw_d(rng_engine));
        double mu = mu_d(rng_engine);
        // 2 - Solve the lambert problem
        kep3::lambert_problem lp(r0, r1, tof, mu, cw, revs_max);

        // 3 - Check the Delta guidance error
        for (const auto &v1 : lp.get_v0()) {
            double dg_err = kep3_tests::delta_guidance_error(r0, r1, v1, mu);
            if (!(dg_err < 1e-12)) {
                std::cout << lp << std::endl;
                fmt::print("\nr1= {}\nr2= {}\ntof= {}\nmu= {}\ncw= {}\nrevs_max= {}", r0, r1, tof, mu, cw, revs_max);
            }
            REQUIRE(dg_err < 1e-12);
        }
    }
}

TEST_CASE("methods")
{
    // Here we test construction for a simple geometry
    kep3::lambert_problem lp{{1., 0., 0.}, {0., 1., 0.}, 3. * kep3::pi / 2., 1., true, 5};
    auto v1 = lp.get_v0()[0];
    auto v2 = lp.get_v1()[0];
    REQUIRE(kep3_tests::floating_point_error_vector(v1, {0, -1, 0}) < 1e-13);
    REQUIRE(kep3_tests::floating_point_error_vector(v2, {1, 0, 0}) < 1e-13);
    auto r1 = lp.get_r0();
    auto r2 = lp.get_r1();
    REQUIRE(r1 == std::array<double, 3>{1, 0, 0});
    REQUIRE(r2 == std::array<double, 3>{0, 1, 0});
    REQUIRE(lp.get_tof() == 3 * kep3::pi / 2);
    REQUIRE(lp.get_mu() == 1.);
    REQUIRE(kep3_tests::floating_point_error(lp.get_x()[0], -0.3826834323650896) < 1e-13);
    REQUIRE(lp.get_iters()[0] == 3u);
    REQUIRE(lp.get_Nmax() == 0u);
    REQUIRE(lp.get_cw() == true);
}

TEST_CASE("serialization_test")
{
    // Instantiate a generic lambert problem
    kep3::lambert_problem lp{
        {1.23, 0.1253232342323, 0.57235553354}, {0.234233423, 1.8645645645, 0.234234234}, 25.254856435, 1., true, 10};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(lp);
    // Now serialize
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << lp;
    }
    // Deserialize
    // Create a new lambert problem object
    kep3::lambert_problem lp2{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> lp2;
    }
    auto after = boost::lexical_cast<std::string>(lp2);
    // Compare the string represetation
    REQUIRE(before == after);
}