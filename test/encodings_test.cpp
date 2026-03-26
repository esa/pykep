// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <kep3/core_astro/encodings.hpp>
#include <stdexcept>

#include "catch.hpp"

TEST_CASE("encodings")
{
    // We only call these here: their correctness is tested in python
    REQUIRE_NOTHROW(kep3::alpha2direct({0.1,0.2,0.3},  21.2));
    REQUIRE_NOTHROW(kep3::direct2alpha({0.1,0.2,0.3}));
    REQUIRE_NOTHROW(kep3::eta2direct({0.1,0.2,0.3},  21.2));
    REQUIRE_NOTHROW(kep3::direct2eta({0.1,0.2,0.3},  21.2));
    // We check throws for invalid dimensions
    REQUIRE_THROWS_AS(kep3::direct2eta({},  21.2), std::invalid_argument);
    REQUIRE_THROWS_AS(kep3::eta2direct({},  21.2), std::invalid_argument);
}