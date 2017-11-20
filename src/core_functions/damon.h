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
#ifndef KEP_TOOLBOX_DAMON_APPROX_H
#define KEP_TOOLBOX_DAMON_APPROX_H

#include <boost/math/special_functions/sign.hpp>
#include <cmath>

#include "../astro_constants.h"
#include "../exceptions.h"

namespace kep_toolbox
{

/// Computes Damon's linear approximation
/**
 * This function uses the model developed by Damon Landau (JPL) during GTOC7 to
 * compute an approximation to the low-thrust transfer
 *
 * \param[in] v1 starting velocity relative to the departure body. This is
 *               computed from the Lambert solution linking initial and
 *               final position in tof
 * \param[in] v2 ending velocity relative to the arrival body. This is
 *               computed from the Lambert solution linking initial and
 *               final position in tof
 * \param[in] tof time of flight.
 *
 * \param[out] a1 acceleration direction in the first part of the transfer
 * \param[out] a2 acceleration direction in the second part of the transfer
 * \param[out] tau duration of the first part of the transfer
 * \param[out] dv total estimated DV
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
inline void damon_approx(const array3D &v1, const array3D &v2, double tof, array3D &a1, array3D &a2, double &tau,
                         double &dv)
{
    array3D v1_nd, v2_nd, U, V;
    double alpha;

    // We compute preliminary quantities
    for (int i = 0; i < 3; ++i)
        v1_nd[i] = v1[i] * tof;
    for (int i = 0; i < 3; ++i)
        v2_nd[i] = v2[i] * tof;
    for (int i = 0; i < 3; ++i)
        U[i] = v2_nd[i] + v1_nd[i];
    for (int i = 0; i < 3; ++i)
        V[i] = v2_nd[i] - v1_nd[i];
    alpha = (U[0] * U[0] + U[1] * U[1] + U[2] * U[2]) / (U[0] * V[0] + U[1] * V[1] + U[2] * V[2]);

    // We then compute a1, a2, tau and dv
    tau = (alpha + 1. - boost::math::sign(alpha) * sqrt(1. + alpha * alpha)) / 2.;
    for (int i = 0; i < 3; ++i)
        a1[i] = V[i] - U[i] / tau;
    for (int i = 0; i < 3; ++i)
        a2[i] = V[i] - U[i] / (tau - 1.);
    dv = (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    dv += (U[0] * U[0] + U[1] * U[1] + U[2] * U[2]) / tau / tau;
    dv -= 2. * (U[0] * V[0] + U[1] * V[1] + U[2] * V[2]) / tau;
    dv = sqrt(dv);

    // We put units back
    for (int i = 0; i < 3; ++i)
        a1[i] = a1[i] / tof / tof;
    for (int i = 0; i < 3; ++i)
        a2[i] = a2[i] / tof / tof;
    tau = tau * tof;
    dv = dv / tof;
}

/// Computes maximum starting mass to convert a Lambert arc to low-thrust
/**
 * This function uses Damon's model to estimate the maximum
 * mass a NEP spacececraft can have in order for a given Lambert arc to be
 * convertable into a valid low-thrust transfer
 *
 * \param[in] a acceleration magnitude as computed from Damon's model
 * \param[in] dv dv magnitude as computed from Damon's model
 * \param[in] T_max maximum thrust of the spacecarft NEP propulsion system
 * \param[in] Isp specific impulse of the spacecarft NEP propulsion system
 *
 * \param[out] a1 acceleration direction in the first part of the transfer
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
inline double max_start_mass(double a, double dv, double T_max, double Isp)
{
    return 2 * T_max / a / (1 + exp(-dv / Isp / ASTRO_G0));
}

} // namespace

#endif // KEP_TOOLBOX_DAMON_APPROX_H
