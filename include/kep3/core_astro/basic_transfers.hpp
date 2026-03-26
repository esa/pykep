// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_IMPULSIVE_TRANSFERS_H
#define kep3_IMPULSIVE_TRANSFERS_H

#include <array>   // For std::array
#include <tuple> // For std::tuple

#include <kep3/detail/visibility.hpp>

namespace kep3
{
kep3_DLL_PUBLIC std::tuple<double, double, std::array<double, 2>> hohmann(double r1, double r2, double mu);

kep3_DLL_PUBLIC std::tuple<double, double, std::array<double, 3>> bielliptic(double r1, double r2, double rb, double mu);
} // namespace kep3
#endif // kep3_IMPULSIVE_TRANSFERS_H