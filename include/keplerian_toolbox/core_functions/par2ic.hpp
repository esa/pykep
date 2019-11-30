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

#ifndef KEP_TOOLBOX_PAR2IC_H
#define KEP_TOOLBOX_PAR2IC_H

#include <cmath>

#include <keplerian_toolbox/io.hpp>

namespace kep_toolbox
{

/// From osculating Keplerian to cartesian
/**
* Transforms Keplerian elements (a,e,i,W,w,E) to cartesian.
* Note that we use the eccentric anomaly (or Gudermannian)
*
* @param[in] E the osculating Keplerian elements (a,e,i,W,w,E) in SI
* @param[in] mu gravitational parameter, (m^3/s^2)
*
* @param[out] r the cartesian position corresponding to the input elements
* @param[out] v the cartesian velocity corresponding to the input elements
*
* Note: both a and e must be positive, otherwise undefined behaviour
*/
template <class vettore3D, class vettore6D>
void par2ic(const vettore6D &E, const double &mu, vettore3D &r0, vettore3D &v0)
{
    double a = E[0];
    double e = E[1];
    double i = E[2];
    double omg = E[3];
    double omp = E[4];
    double EA = E[5];
    double b, n, xper, yper, xdotper, ydotper;
    double R[3][3];
    double cosomg, cosomp, sinomg, sinomp, cosi, sini;
    double dNdZeta;

    // semi-major axis is assumed to be positive here we apply the convention of having it negative as for
    // computations to result in higher elegance
    if (e > 1) {
        a = -a;
    }

    // 1 - We start by evaluating position and velocity in the perifocal reference system
    if (e < 1.0) // EA is the eccentric anomaly
    {
        b = a * sqrt(1 - e * e);
        n = sqrt(mu / (a * a * a));

        xper = a * (cos(EA) - e);
        yper = b * sin(EA);
        xdotper = -(a * n * sin(EA)) / (1 - e * cos(EA));
        ydotper = (b * n * cos(EA)) / (1 - e * cos(EA));
    } else // EA is the Gudermannian
    {
        b = -a * sqrt(e * e - 1);
        n = sqrt(-mu / (a * a * a));

        dNdZeta = e * (1 + tan(EA) * tan(EA)) - (0.5 + 0.5 * pow(tan(0.5 * EA + M_PI_4), 2)) / tan(0.5 * EA + M_PI_4);

        xper = a / cos(EA) - a * e;
        yper = b * tan(EA);

        xdotper = a * tan(EA) / cos(EA) * n / dNdZeta;
        ydotper = b / pow(cos(EA), 2) * n / dNdZeta;
    }

    // 2 - We then built the rotation matrix from perifocal reference frame to inertial

    cosomg = cos(omg);
    cosomp = cos(omp);
    sinomg = sin(omg);
    sinomp = sin(omp);
    cosi = cos(i);
    sini = sin(i);

    R[0][0] = cosomg * cosomp - sinomg * sinomp * cosi;
    R[0][1] = -cosomg * sinomp - sinomg * cosomp * cosi;
    R[0][2] = sinomg * sini;
    R[1][0] = sinomg * cosomp + cosomg * sinomp * cosi;
    R[1][1] = -sinomg * sinomp + cosomg * cosomp * cosi;
    R[1][2] = -cosomg * sini;
    R[2][0] = sinomp * sini;
    R[2][1] = cosomp * sini;
    R[2][2] = cosi;

    // 3 - We end by transforming according to this rotation matrix

    double temp[3] = {xper, yper, 0.0};
    double temp2[3] = {xdotper, ydotper, 0};

    for (int j = 0; j < 3; j++) {
        r0[j] = 0.0;
        v0[j] = 0.0;
        for (int k = 0; k < 3; k++) {
            r0[j] += R[j][k] * temp[k];
            v0[j] += R[j][k] * temp2[k];
        }
    }
    return;
}
}
#endif // KEP_TOOLBOX_PAR2IC_H
