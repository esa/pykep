// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <optional>
#include <stdexcept>
#include <string>

#include <fmt/core.h>

#include "catch.hpp"
#include <boost/core/demangle.hpp>
#include <boost/lexical_cast.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/jpl_lp.hpp>
#include <kep3/udpla/keplerian.hpp>

#include "test_helpers.hpp"

using kep3::epoch;
using kep3::planet;
using kep3::detail::null_udpla;

struct simple_udpla {
    simple_udpla() = default;
    static std::array<std::array<double, 3>, 2> eph(double)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever!";
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

KEP3_S11N_EXPORT_WRAP(simple_udpla, kep3::detail::planet_iface)

struct simple_udpla_mu {
    simple_udpla_mu() = default;
    static std::array<std::array<double, 3>, 2> eph(double)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet with mu";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever but with mu!";
    }
    static double get_mu_central_body()
    {
        return 1.;
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

KEP3_S11N_EXPORT_WRAP(simple_udpla_mu, kep3::detail::planet_iface)

struct simple_udpla_mu_h {
    simple_udpla_mu_h() = default;
    static std::array<std::array<double, 3>, 2> eph(double)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 10., 0.};
        return {pos, vel};
    };
    static std::string get_name()
    {
        return "A simple planet with mu, but hyperbolic conditions";
    }
    static std::string get_extra_info()
    {
        return "The simplest planet ever but with mu!";
    }
    static double get_mu_central_body()
    {
        return 1.;
    }

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

KEP3_S11N_EXPORT_WRAP(simple_udpla_mu_h, kep3::detail::planet_iface)

struct complete_udpla {
    explicit complete_udpla(std::array<double, 4> physical_properties = {-1., -1., -1., -1.})
        : m_name("A complete, albeit simple Planet"), m_mu_central_body(physical_properties[0]),
          m_mu_self(physical_properties[1]), m_radius(physical_properties[2]), m_safe_radius(physical_properties[3]) {};
    static std::array<std::array<double, 3>, 2> eph(double)
    {
        std::array<double, 3> pos = {1., 0., 0.};
        std::array<double, 3> vel = {0., 1., 0.};
        return {pos, vel};
    };

    static std::vector<double> eph_v(const std::vector<double> &)
    {
        return {1., 0., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0.};
    };

    static std::array<double, 3> acc(double &)
    {
        return std::array<double, 3>{0., 0., 0.};
    };

    [[nodiscard]] std::string get_name() const
    {
        return m_name;
    }

    [[nodiscard]] double get_mu_central_body() const
    {
        return m_mu_central_body;
    }
    [[nodiscard]] double get_mu_self() const
    {
        return m_mu_self;
    }
    [[nodiscard]] double get_radius() const
    {
        return m_radius;
    }
    [[nodiscard]] double get_safe_radius() const
    {
        return m_safe_radius;
    }
    [[nodiscard]] double period(double) const
    {
        return m_radius - m_radius;
    }

    [[nodiscard]] std::array<double, 6> elements(double, kep3::elements_type) const
    {
        return {m_radius, m_radius, m_radius, m_radius, m_radius, m_radius};
    }

    [[nodiscard]] std::string get_extra_info() const
    {
        return fmt::format("Gravitational parameter of main attracting body: {}\nGravitational "
                           "parameter of body: {}\nBody radius: {}\nBody safe radius: {}\n",
                           get_mu_central_body(), get_mu_self(), get_radius(), get_safe_radius());
    }

    std::string m_name;
    double m_mu_central_body;
    double m_mu_self;
    double m_radius;
    double m_safe_radius;

private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar & m_name;
        ar & m_mu_central_body;
        ar & m_mu_self;
        ar & m_radius;
        ar & m_safe_radius;
    }
};

KEP3_S11N_EXPORT_WRAP(complete_udpla, kep3::detail::planet_iface)

TEST_CASE("construction")
{
    {
        kep3::detail::null_udpla udpla{};
        udpla.eph(0.);
        // Default constructor (a null planet)
        REQUIRE_NOTHROW(planet());
        auto pla = planet();
        REQUIRE_NOTHROW(pla.eph(0.));
        REQUIRE(!kep3::detail::udpla_has_acc<kep3::detail::null_udpla>);
        auto pos_vel = pla.eph(0.);
        REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
        REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
        REQUIRE(pla.get_name() == boost::core::demangle(typeid(null_udpla).name()));
        REQUIRE(pla.get_extra_info() == std::string(""));
        REQUIRE(pla.get_mu_central_body() == -1);
        REQUIRE(pla.get_mu_self() == -1);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == -1);
        REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.elements()), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.acc(0.)), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.acc(kep3::epoch(0.))), kep3::not_implemented_error);

        auto eph1 = pla.eph(0.);
        auto eph2 = pla.eph(2.);
        REQUIRE(pla.eph_v(std::vector<double>{1., 2.})
                == std::vector<double>{eph1[0][0], eph1[0][1], eph1[0][2], eph1[1][0], eph1[1][1], eph1[1][2],
                                       eph2[0][0], eph2[0][1], eph2[0][2], eph2[1][0], eph2[1][1], eph2[1][2]});
        REQUIRE(!pla.extract<complete_udpla>());
        REQUIRE(pla.extract<null_udpla>());
    }
    {
        // Constructor from a simple udpla
        simple_udpla udpla{};
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        auto pos_vel = pla.eph(epoch(0.));
        REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
        REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
        REQUIRE(pla.get_name() == "A simple planet");
        REQUIRE(pla.get_extra_info() == "The simplest planet ever!");
        REQUIRE(pla.get_mu_central_body() == -1);
        REQUIRE(pla.get_mu_self() == -1);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == -1);
        REQUIRE_THROWS_AS((pla.period()), kep3::not_implemented_error);
        REQUIRE_THROWS_AS((pla.elements()), kep3::not_implemented_error);
        auto eph1 = pla.eph(0.);
        auto eph2 = pla.eph(2.);
        REQUIRE(pla.eph_v(std::vector<double>{1., 2.})
                == std::vector<double>{eph1[0][0], eph1[0][1], eph1[0][2], eph1[1][0], eph1[1][1], eph1[1][2],
                                       eph2[0][0], eph2[0][1], eph2[0][2], eph2[1][0], eph2[1][1], eph2[1][2]});
    }
    {
        // Constructor from a more complete udpla
        complete_udpla udpla({1., 2., -1., 4.});
        REQUIRE_NOTHROW(planet(udpla));
        REQUIRE(kep3::detail::udpla_has_acc<complete_udpla>);
        planet pla(udpla);
        auto pos_vel = pla.eph(epoch(0.));
        REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
        REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
        REQUIRE(pla.acc(0.) == std::array<double, 3>{0., 0., 0.});
        REQUIRE(pla.acc(kep3::epoch(0.)) == std::array<double, 3>{0., 0., 0.});
        REQUIRE(pla.get_name() == "A complete, albeit simple Planet");
        REQUIRE(pla.get_mu_central_body() == 1.);
        REQUIRE(pla.get_mu_self() == 2.);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == 4.);
        REQUIRE(std::isfinite(pla.period()));
        REQUIRE(pla.elements() == std::array<double, 6>{-1, -1, -1, -1, -1, -1});
        auto eph1 = pla.eph(0.);
        auto eph2 = pla.eph(2.);
        REQUIRE(pla.eph_v(std::vector<double>{1., 2.})
                == std::vector<double>{eph1[0][0], eph1[0][1], eph1[0][2], eph1[1][0], eph1[1][1], eph1[1][2],
                                       eph2[0][0], eph2[0][1], eph2[0][2], eph2[1][0], eph2[1][1], eph2[1][2]});
    }
    {
        // Constructor from a simple udpla with mu and hyperbolic orbit
        simple_udpla_mu_h udpla{};
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        REQUIRE(!std::isfinite(pla.period()));
        REQUIRE(pla.elements()[0] < 0);
        REQUIRE_THROWS_AS(pla.elements(kep3::epoch(), kep3::elements_type::KEP_M), std::logic_error);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEE)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEE_R)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_F)[0] < 0);
    }
    {
        // Constructor from a simple udpla with mu and elliptical orbit
        simple_udpla_mu udpla{};
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        REQUIRE(kep3_tests::floating_point_error(pla.period(), kep3::pi * 2.) < 1e-14);
        REQUIRE(pla.elements()[0] > 0);
        REQUIRE(pla.elements()[1] < 1);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_M)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEE)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::MEE_R)[0] > 0);
        REQUIRE(pla.elements(kep3::epoch(), kep3::elements_type::KEP_F)[0] > 0);
    }
    {
        // Constructor from a more complete udpla
        complete_udpla udpla({1., 2., -1., 4.});
        REQUIRE_NOTHROW(planet(udpla));
        planet pla(udpla);
        auto pos_vel = pla.eph(epoch(0.));
        REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
        REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
        REQUIRE(pla.get_name() == "A complete, albeit simple Planet");
        REQUIRE(pla.get_mu_central_body() == 1.);
        REQUIRE(pla.get_mu_self() == 2.);
        REQUIRE(pla.get_radius() == -1);
        REQUIRE(pla.get_safe_radius() == 4.);
    }
    // Check copy semantics.
    auto p0 = planet();
    planet p1{p0};
    REQUIRE(p0.extract<null_udpla>());
    REQUIRE(p1.extract<null_udpla>());
    planet p2{simple_udpla{}};
    p2 = p1;
    REQUIRE(p2.extract<null_udpla>());
    //  Move semantics.
    planet p3{std::move(p0)};
    REQUIRE(p3.extract<null_udpla>());
    planet p4{simple_udpla{}};
    p4 = std::move(p2);
    REQUIRE(p4.extract<null_udpla>());
    //  Check we can revive moved-from objects.
    p0 = p4;
    REQUIRE(p0.extract<null_udpla>());
    p2 = std::move(p4);
    REQUIRE(p2.extract<null_udpla>());
}

TEST_CASE("planet_extract_test")
{
    // We instantiate a planet
    planet pla{complete_udpla({1., 2., -1., 4.})};

    auto *p0 = pla.extract<complete_udpla>();

    // We check thet we can access to public data members
    REQUIRE(p0->m_mu_central_body == 1.);
    REQUIRE(p0->m_mu_self == 2.);
    REQUIRE(p0->m_radius == -1.);
    REQUIRE(p0->m_safe_radius == 4.);

    // We check the is method
    REQUIRE(pla.extract<complete_udpla>());
    REQUIRE(!pla.extract<simple_udpla>());
}

TEST_CASE("generic_assignment")
{
    REQUIRE((!std::is_assignable<planet, void>::value));
    REQUIRE((!std::is_assignable<planet, int &>::value));
    REQUIRE((!std::is_assignable<planet, const int &>::value));
    REQUIRE((!std::is_assignable<planet, int &&>::value));
}

TEST_CASE("stream_operator")
{
    REQUIRE_NOTHROW((std::cout << planet{} << '\n'));
}

TEST_CASE("planet_astro_methods_test")
{
    // We instantiate a planet
    planet pla{complete_udpla({1.1, 2.3, 4.02, 4.5})};
    // Test eph
    auto pos_vel = pla.eph(epoch(0.));
    REQUIRE(pos_vel[0] == std::array<double, 3>{1., 0., 0.});
    REQUIRE(pos_vel[1] == std::array<double, 3>{0., 1., 0.});
    REQUIRE(pla.get_name() == "A complete, albeit simple Planet");
    REQUIRE(pla.get_mu_central_body() == 1.1);
    REQUIRE(pla.get_mu_self() == 2.3);
    REQUIRE(pla.get_radius() == 4.02);
    REQUIRE(pla.get_safe_radius() == 4.5);
    REQUIRE(pla.period(kep3::epoch(0.)) == 0.);
}

TEST_CASE("acc_test")
{
    // A keplerian udpla
    kep3::udpla::keplerian udpla_kep{
        kep3::epoch(0.), {{{0.33, 1.3, 0.12}, {0.01, 1.123, 0.2}}}, 1.12, "enterprise", {12.32, 44.6, 98.23}};
    planet pla_kep{udpla_kep};
    // A custom udpla with acc implemented
    planet pla_complete{complete_udpla({1.1, 2.3, 4.02, 4.5})};
    // A custom udpla with mu implemented but not acc
    planet pla_mu{simple_udpla_mu{}};
    // A custom udpla with nothing
    planet pla_simple{simple_udpla{}};
    // 1 - Throws
    REQUIRE_THROWS_AS((pla_simple.acc(12.3)), kep3::not_implemented_error);
    REQUIRE_NOTHROW(pla_kep.acc(12.3));
    REQUIRE_NOTHROW(pla_complete.acc(12.3));
    REQUIRE_NOTHROW(pla_mu.acc(12.3));
    // 2 - The user implements acc ... we check that the returned value is the one implemented by the user.
    auto acc_complete = pla_complete.acc(0.);
    REQUIRE(acc_complete == std::array<double, 3>{0., 0., 0.});
    // 2 - The user implements mu ... we check that the returned value is the default keplerian.
    {
        auto mjd2000 = 212323423.23234234234;
        const auto &[pos, vel] = pla_mu.eph(mjd2000);
        auto acc_mu = pla_mu.acc(mjd2000);
        auto mu = pla_mu.get_mu_central_body();
        auto mur3 = -mu / std::pow(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2], 1.5); // -mu/r^3
        std::array<double, 3> gt = {mur3 * pos[0], mur3 * pos[1], mur3 * pos[2]};
        REQUIRE(acc_mu == gt);
    }
    // 3 - The user an udpla shipped in pykep ... we check that the returned value is the default.
    {
        auto mjd2000 = 212323423.23234234234;
        const auto &[pos, vel] = pla_kep.eph(mjd2000);
        auto acc_mu = pla_kep.acc(mjd2000);
        auto mu = pla_kep.get_mu_central_body();
        auto mur3 = -mu / std::pow(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2], 1.5); // -mu/r^3
        std::array<double, 3> gt = {mur3 * pos[0], mur3 * pos[1], mur3 * pos[2]};
        REQUIRE(acc_mu == gt);
    }
    // 4 - The vectorized (defaulted) implementation
    std::vector<double> mjd2000s = {12.,34.,0.03,-323.231};
    auto res = pla_kep.acc_v(mjd2000s);
    for (auto i = 0u; i < 4u; ++i) {
        auto gt = pla_kep.acc(mjd2000s[i]);
        REQUIRE(res[3*i] == gt[0]);
        REQUIRE(res[3*i+1] == gt[1]);
        REQUIRE(res[3*i+2] == gt[2]);
    }
}

TEST_CASE("serialization_test")
{
    // Instantiate a planet
    planet pla{complete_udpla({1.1, 2.3, 4.02, 4.5})};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(pla);
    // Now serialize, deserialize and compare the result.
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << pla;
    }
    // Create a new planet object
    auto pla2 = planet{simple_udpla{}};
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

TEST_CASE("serialization_test_2")
{
    // Instantiate a planet
    planet pla{null_udpla{}};

    // Store the string representation.
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(pla);
    // Now serialize, deserialize and compare the result.
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << pla;
    }
    // Create a new planet object
    auto pla2 = planet{simple_udpla{}};
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
