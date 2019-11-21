/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef KEP_TOOLBOX_EQ2IC_H
#define KEP_TOOLBOX_EQ2IC_H

#include <cmath>
#include <keplerian_toolbox/io.hpp>

namespace kep_toolbox
{

/// From equinoctial to cartesian
/**
* Transforms equinoctial elements (p,f,g,h,k,L) to cartesian.
* Note that we use the true longitude LF (i.e. L + om + Om)
*
* @param[in] EQ the the equinoctial elements (p,f,g,h,k,L) in SI.
* @param[in] mu gravitational parameter, (m^3/s^2)
*
* @param[out] r the cartesian position corresponding to the input elements
* @param[out] v the cartesian velocity corresponding to the input elements
*
*/
template <class vettore3D, class vettore6D>
void eq2ic(const vettore6D &EQ, const double &mu, vettore3D &r, vettore3D &v, const bool retrogade = false)
{
    int I;
    if (retrogade) {
        I = -1;
    } else {
        I = 1;
    }

    // p = a (1-e^2) will be negative for eccentricities > 1, we here need a positive number for the following
    // computations
    // to make sense
    auto par = std::abs(EQ[0]);
    auto f = EQ[1];
    auto g = EQ[2];
    auto h = EQ[3];
    auto k = EQ[4];
    auto L = EQ[5];

    // We compute the equinoctial reference frame
    auto den = k * k + h * h + 1;
    auto fx = (1 - k * k + h * h) / den;
    auto fy = (2 * k * h) / den;
    auto fz = (-2 * I * k) / den;

    auto gx = (2 * I * k * h) / den;
    auto gy = (1 + k * k - h * h) * I / den;
    auto gz = (2 * h) / den;

    auto wx = (2 * k) / den;
    auto wy = (-2 * h) / den;
    auto wz = (1 - k * k - h * h) * I / den;

    // Auxiliary
    auto b = std::sqrt(1 - g * g - f * f);
    auto radius = par / (1 + g * std::sin(L) + f * std::cos(L));
    // In the equinoctial reference frame
    auto X = radius * std::cos(L);
    auto Y = radius * std::sin(L);
    auto VX = -std::sqrt(mu / par) * (g + std::sin(L));
    auto VY = std::sqrt(mu / par) * (f + std::cos(L));

/*print("r is: ", radius, "\n");
print("f is: [", fx, ", ", fy, ", ", fz, "]\n");
print("g is: [", gx, ", ", gy, ", ", gz, "]\n");
print("X is: ", X , "\n");
print("Y is: ", Y , "\n");*/
    // Results
    r[0] = X * fx + Y * gx;
    r[1] = X * fy + Y * gy;
    r[2] = X * fz + Y * gz;

    v[0] = VX * fx + VY * gx;
    v[1] = VX * fy + VY * gy;
    v[2] = VX * fz + VY * gz;

/*print("rx: ", r[0], "\n");
print("ry: ", r[1], "\n\n");*/

    return;
}
}
#endif // KEP_TOOLBOX_PAR2IC_H
