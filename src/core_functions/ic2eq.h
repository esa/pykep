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

#ifndef KEP_TOOLBOX_IC2EQ_H
#define KEP_TOOLBOX_IC2EQ_H

#include <cmath>

#include "array3D_operations.h"
#include "../io.h"

namespace kep_toolbox
{
/// From cartesian to modified equinoctial
/**
* Transforms cartesian coordinates (r,v) to modified equinoctial elements (p,h,k,p,q,L).
* Note that we use the true longitude L (i.e. f + om + Om).
*
* @param[in] r0 cartesian position vector.
* @param[in] v0 cartesian velocity vector.
* @param[in] mu gravitational parameter.
* @param[in] retrogade forces retrograde elements to be used
*
* @param[out] EQ the modified equinoctial elements (p,h,k,p,q,L).
*
*/
template <class vettore3D, class vettore6D>
void ic2eq(const vettore3D &r0, const vettore3D &v0, const double &mu, vettore6D &EQ, const bool retrogade = false)
{
    int I;
    if (retrogade) {
        I = -1;
    } else {
        I = 1;
    }
    // The equinoctial reference frame
    vettore3D f = {0.0, 0.0, 0.0};
    vettore3D g = {0.0, 0.0, 0.0};
    vettore3D w = {0.0, 0.0, 0.0};
    // The eccentricity vector
    vettore3D e = {0.0, 0.0, 0.0};
    // angular momentum
    vettore3D ang = {0.0, 0.0, 0.0};
    cross(ang, r0, v0);

    // 0 - We compute the semi-major axis
    auto R0 = norm(r0);
    auto V0 = norm(v0);
    auto a = std::abs(1. / (2. / R0 - V0 * V0 / mu));

    // 1 - We compute the equinoctial frame
    cross(w, r0, v0);
    vers(w, w);

    auto p = w[0] / (1 + I * w[2]);
    auto q = -w[1] / (1 + I * w[2]);
    auto den = p * p + q * q + 1;
    f[0] = (1. - p * p + q * q) / den;
    f[1] = (2. * p * q) / den;
    f[2] = (-2. * I * p) / den;

    g[0] = (2. * I * p * q) / den;
    g[1] = (1. + p * p - q * q) * I / den;
    g[2] = (2. * q) / den;

    // 2 - We compute evett: the eccentricity vector
    vettore3D Dum_Vec = {0.0, 0.0, 0.0};
    cross(Dum_Vec, v0, ang);
    for (auto i = 0u; i < 3; ++i) {
        e[i] = Dum_Vec[i] / mu - r0[i] / R0;
    }
    auto h = dot(e, g);
    auto k = dot(e, f);
    auto ecc = norm(e);

    // 3 - We compute the true longitude L
    // This solution is certainly not the most elegant, but it works and will never be singular.

    auto det1 = (g[1]*f[0]-f[1]*g[0]); // xy
    auto det2 = (g[2]*f[0]-f[2]*g[0]); // xz
    auto det3 = (g[2]*f[1]-f[2]*g[1]); // yz
    auto max = std::max({std::abs(det1), std::abs(det2), std::abs(det3)});

    double X, Y;
    if (std::abs(det1) == max) {
        X = (g[1]*r0[0] - g[0] * r0[1]) / det1;
        Y = (-f[1]*r0[0]+f[0]*r0[1]) / det1;
    } else if (std::abs(det2) == max) {
        X = (g[2]*r0[0] - g[0] * r0[2]) / det2;
        Y = (-f[2]*r0[0]+f[0]*r0[2]) / det2;
    } else {
        X = (g[2]*r0[1] - g[1] * r0[2]) / det3;
        Y = (-f[2]*r0[1]+f[1]*r0[2]) / det3;
    }

    auto L = std::atan2(Y/R0, X/R0);

    // 5 - We assign the results
    EQ[0] = a * (1. - ecc * ecc);
    EQ[1] = h;
    EQ[2] = k;
    EQ[3] = p;
    EQ[4] = q;
    EQ[5] = L;
}
} // namespace kep_toolbox end

#endif // KEP_TOOLBOX_IC2EQ_H
