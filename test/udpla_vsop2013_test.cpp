// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/algorithm/string/predicate.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/udpla/jpl_lp.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "catch.hpp"

using kep3::planet;
using kep3::udpla::vsop2013;

TEST_CASE("basic")
{
    using Catch::Matchers::Message;

    {
        std::ostringstream oss;

        vsop2013 pl1;
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "mercury"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
    }

    {
        std::ostringstream oss;

        vsop2013 pl1("vEnUs");
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
    }

    {
        std::ostringstream oss;

        vsop2013 pl1("venus", 0.5);
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }

    // Run it again to hit the cache.
    {
        std::ostringstream oss;

        vsop2013 pl1("venus", 0.5);
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }

    REQUIRE_THROWS_MATCHES(vsop2013("venus", -1), std::invalid_argument,
                           Message("The threshold value must be finite and non-negative, but it is -1 instead"));

    REQUIRE_THROWS_MATCHES(vsop2013("goofy"), std::invalid_argument, Message("The planet name 'goofy' is invalid"));

    // Move ctor.
    {
        std::ostringstream oss;

        vsop2013 pl0("venus", 0.5);
        auto pl1 = std::move(pl0);
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }

    // Copy ctor.
    {
        std::ostringstream oss;

        vsop2013 pl0("venus", 0.5);
        auto pl1 = pl0;
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));

        oss.str("");

        oss << planet(pl0);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }

    // Move assignment.
    {
        std::ostringstream oss;

        vsop2013 pl0("venus", 0.5), pl1;
        pl1 = std::move(pl0);
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }

    // Copy assignment.
    {
        std::ostringstream oss;

        vsop2013 pl0("venus", 0.5), pl1;
        pl1 = pl0;
        oss << planet(pl1);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));

        oss.str("");

        oss << planet(pl0);

        REQUIRE(boost::contains(oss.str(), "venus"));
        REQUIRE(boost::contains(oss.str(), "thresh"));
        REQUIRE(boost::contains(oss.str(), "0.5"));
    }
}

TEST_CASE("s11n")
{
    std::stringstream ss;

    vsop2013 pl0("jupiter", 0.5);

    {
        boost::archive::binary_oarchive oa(ss);
        oa << pl0;
    }

    pl0 = vsop2013{};

    {
        boost::archive::binary_iarchive ia(ss);
        ia >> pl0;
    }

    std::ostringstream oss;
    oss << planet(pl0);

    REQUIRE(boost::contains(oss.str(), "jupiter"));
    REQUIRE(boost::contains(oss.str(), "thresh"));
    REQUIRE(boost::contains(oss.str(), "0.5"));
}

TEST_CASE("s11n_2")
{
    // Instantiate a planet with jpl_lp udpla
    kep3::epoch ref_epoch{2423.4343, kep3::epoch::julian_type::MJD2000};
    vsop2013 udpla("jupiter", 0.5);
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

TEST_CASE("eph")
{
    planet p{vsop2013{"venus", 1e-9}};

    auto eph = p.eph(123.);

    // NOTE: these values have been manually checked against DE440.
    REQUIRE(eph[0][0] == Approx(103304986899.7981));
    REQUIRE(eph[0][1] == Approx(32220404104.1199));
    REQUIRE(eph[0][2] == Approx(7957719449.51538));

    REQUIRE(eph[1][0] == Approx(-10696.505905435035));
    REQUIRE(eph[1][1] == Approx(30061.035989651813));
    REQUIRE(eph[1][2] == Approx(14201.00090492195));
}
