// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>

namespace kep3
{

std::array<double, 6> mee2par(const std::array<double, 6> &mee, bool retrogade)
{
    std::array<double, 6> retval{};
    int I = 1;
    if (retrogade) {
        I = -1;
    }
    double ecc = std::sqrt(mee[1] * mee[1] + mee[2] * mee[2]);
    double tmp = std::sqrt(mee[3] * mee[3] + mee[4] * mee[4]);
    double zita = std::atan2(mee[2] / ecc, mee[1] / ecc); // [-pi, pi]
    if (zita < 0) {
        zita += 2 * pi; // [0, 2*pi]
    }

    retval[1] = ecc;
    retval[0] = mee[0] / (1. - ecc * ecc);
    retval[2] = half_pi * (1. - I) + 2. * I * std::atan(tmp);
    retval[3] = std::atan2(mee[4] / tmp, mee[3] / tmp); // [-pi, pi]
    if (retval[3] < 0) {
        retval[3] += 2 * pi; // [0, 2*pi]
    }
    retval[4] = zita - I * retval[3]; //
    if (retval[4] < 0) {
        retval[4] += 2 * pi;
    } else if (retval[4] > 2 * pi) {
        retval[4] -= 2 * pi;
    }
    retval[5] = mee[5] - I * retval[3] - retval[4];
    return retval;
}

std::array<double, 6> par2mee(const std::array<double, 6> &par, bool retrogade)
{
    std::array<double, 6> mee{};
    int I = 0;
    if (retrogade) {
        I = -1;
        mee[3] = 1. / std::tan(par[2] / 2) * std::cos(par[3]);
        mee[4] = 1. / std::tan(par[2] / 2) * std::sin(par[3]);
    } else {
        I = 1;
        mee[3] = std::tan(par[2] / 2) * std::cos(par[3]);
        mee[4] = std::tan(par[2] / 2) * std::sin(par[3]);
    }
    mee[0] = par[0] * (1 - par[1] * par[1]);
    mee[1] = par[1] * std::cos(par[4] + I * par[3]);
    mee[2] = par[1] * std::sin(par[4] + I * par[3]);
    mee[5] = par[5] + par[4] + I * par[3];
    return mee;
}
} // namespace kep3