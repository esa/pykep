// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>

namespace kep3
{

std::array<double, 6> mee2par(const std::array<double, 6> &eq, bool retrogade)
{
    std::array<double, 6> retval{};
    int I = 1;
    if (retrogade) {
        I = -1;
    }
    double ecc = std::sqrt(eq[1] * eq[1] + eq[2] * eq[2]);
    double tmp = std::sqrt(eq[3] * eq[3] + eq[4] * eq[4]);
    double zita = std::atan2(eq[2] / ecc, eq[1] / ecc); // [-pi, pi]
    if (zita < 0) {
        zita += 2 * pi; // [0, 2*pi]
    }

    retval[1] = ecc;
    retval[0] = eq[0] / (1. - ecc * ecc);
    retval[2] = half_pi * (1. - I) + 2. * I * std::atan(tmp);
    retval[3] = std::atan2(eq[4] / tmp, eq[3] / tmp); // [-pi, pi]
    if (retval[3] < 0) {
        retval[3] += 2 * pi; // [0, 2*pi]
    }
    retval[4] = zita - I * retval[3]; //
    if (retval[4] < 0) {
        retval[4] += 2 * pi;
    } else if (retval[4] > 2 * pi) {
        retval[4] -= 2 * pi;
    }
    retval[5] = eq[5] - I * retval[3] - retval[4];
    return retval;
}

std::array<double, 6> par2mee(const std::array<double, 6> &par, bool retrogade)
{
    std::array<double, 6> eq{};
    int I = 0;
    if (retrogade) {
        I = -1;
        eq[3] = 1. / std::tan(par[2] / 2) * std::cos(par[3]);
        eq[4] = 1. / std::tan(par[2] / 2) * std::sin(par[3]);
    } else {
        I = 1;
        eq[3] = std::tan(par[2] / 2) * std::cos(par[3]);
        eq[4] = std::tan(par[2] / 2) * std::sin(par[3]);
    }
    eq[0] = par[0] * (1 - par[1] * par[1]);
    eq[1] = par[1] * std::cos(par[4] + I * par[3]);
    eq[2] = par[1] * std::sin(par[4] + I * par[3]);
    eq[5] = par[5] + par[4] + I * par[3];
    return eq;
}
} // namespace kep3