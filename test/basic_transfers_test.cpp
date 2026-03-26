// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <kep3/core_astro/basic_transfers.hpp>

#include "catch.hpp"

TEST_CASE("hohmann")
{
    // Ground truths from (https://www.omnicalculator.com/physics/hohmann-transfer)
    {
        auto [dv_total, t_transfer, dvs] = kep3::hohmann(1., 2., 1.);
        REQUIRE(dv_total == Approx(0.28445705).epsilon(1e-6));
        REQUIRE(dvs[0] == Approx(0.15470054).epsilon(1e-6));
        REQUIRE(dvs[1] == Approx(0.1297565).epsilon(1e-6));
    }
    {
        auto [dv_total, t_transfer, dvs] = kep3::hohmann(1.1, 2.2, 1.3);
        REQUIRE(dv_total == Approx(0.3092374).epsilon(1e-6));
        REQUIRE(dvs[0] == Approx(0.1681772).epsilon(1e-6));
        REQUIRE(dvs[1] == Approx(0.1410602).epsilon(1e-6));
    }
}

TEST_CASE("bielliptic")
{
    // Ground truths from (https://www.omnicalculator.com/physics/hohmann-transfer)
    {
        auto [dv_total, t_transfer, dvs] = kep3::bielliptic(1., 2., 2., 1.);
        REQUIRE(dv_total == Approx(0.28445705).epsilon(1e-6));
        REQUIRE(dvs[0] == Approx(0.15470054).epsilon(1e-6));
        REQUIRE(dvs[1] == Approx(0.1297565).epsilon(1e-6));
    }
    {
        auto [dv_total, t_transfer, dvs] = kep3::bielliptic(1.1, 2.2, 2.2, 1.3);
        REQUIRE(dv_total == Approx(0.3092374).epsilon(1e-6));
        REQUIRE(dvs[0] == Approx(0.1681772).epsilon(1e-6));
        REQUIRE(dvs[1] == Approx(0.1410602).epsilon(1e-6));
    }
}