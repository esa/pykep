// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_SPECIAL_FUNCTIONS_H
#define kep3_SPECIAL_FUNCTIONS_H

#include <cmath>

namespace kep3
{

inline double stumpff_s(const double x)
{
    if (x > 0) {
        return (std::sqrt(x) - std::sin(std::sqrt(x))) / std::pow(std::sqrt(x), 3);
    } else if (x < 0) {
        return (std::sinh(std::sqrt(-x)) - std::sqrt(-x)) / std::pow(-x, 3. / 2);
    } else {
        return (1. / 6.);
    }
}

inline double stumpff_c(const double x)
{
    if (x > 0) {
        return (1 - std::cos(std::sqrt(x))) / x;
    } else if (x < 0) {
        return (std::cosh(std::sqrt(-x)) - 1) / (-x);
    } else {
        return 0.5;
    }
}
} // namespace kep3

#endif // kep3_SPECIAL_FUNCTIONS_H