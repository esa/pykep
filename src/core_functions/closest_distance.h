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

#ifndef KEP_TOOLBOX_CLOSEST_DISTANCE_H
#define KEP_TOOLBOX_CLOSEST_DISTANCE_H

#include<cmath>

/// Lagrangian propagation
/**
 * This template function returns the closest distance along a keplerian orbit defined by r1,v1,r2,v2. This can be or not the
 * periplanet, depending on the geometry.
 *
 * \param[out] d_min minimum distance along the defined orbit
 * \param[out] ra apoapsis radius
 * \param[in] r0 initial position vector (of dimension 3)
 * \param[in] v0 initial velocity vector (of dimension 3)
 * \param[in] r1 final position vector (of dimension 3)
 * \param[in] v1 final velocity vector (of dimension 3)
 * \param[in] mu central body gravitational parameter
 *
 * NOTE: multiple revolutions are not accounted for
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

namespace kep_toolbox {
template<class vettore3D>
void closest_distance(double& d_min, double& ra, const vettore3D& r0, const vettore3D& v0, const vettore3D& r1, const vettore3D& v1, const double &mu)
{
    double h[3];
    double Dum_Vec[3];
    double evett[3];

    double p = 0.0;
    double temp =0.0;
    double R0, R1, ni0, ni1, rp, e;
    int i;
    
    //1 - We compute h: the orbital angular momentum vector
    h[0] = r0[1]*v0[2] - r0[2]*v0[1];
    h[1] = r0[2]*v0[0] - r0[0]*v0[2];
    h[2] = r0[0]*v0[1] - r0[1]*v0[0];

    //2 - We compute p: the orbital parameter
    p = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    p /= mu;

    //3 - We compute evett: the eccentricity vector
    R0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] +r0[2]*r0[2]);
    R1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] +r1[2]*r1[2]);
    Dum_Vec[0] = v0[1]*h[2] - v0[2]*h[1];
    Dum_Vec[1] = v0[2]*h[0] - v0[0]*h[2];
    Dum_Vec[2] = v0[0]*h[1] - v0[1]*h[0];
    for (i=0; i<3; i++)
        evett[i] = Dum_Vec[i]/mu - r0[i]/R0;

    //4 - The eccentricity is calculated and stored as the second orbital element
    e  = sqrt(evett[0]*evett[0] + evett[1]*evett[1] +evett[2]*evett[2]);
    
    //5 - In a circular orbit everything is easier
    if (e < 1e-12) {
	d_min = R0;
	ra = R0;
	return;
    }

    //6 - The periplanet is now calculated 
    rp = p/(1+e);
    ra = p/(1-e);

    //7 - We now compute the true anomaly relative to r0 (in 0-2pi)
    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=evett[i]*r0[i];
    ni0 = acos(temp/e/R0);
    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=r0[i]*v0[i];
    if (temp<0.0) ni0 = 2*M_PI - ni0;
    
    //8 - We now compute the true anomaly relative to r1 (in 0-2pi)
    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=evett[i]*r1[i];
    ni1 = acos(temp/e/R1);
    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=r1[i]*v1[i];
    if (temp<0.0) ni1 = 2*M_PI - ni1;

    if (ni0>ni1) //the periplanet is transversed in the zero-th revolution
    {
	d_min = rp;
    } else {//the periplanet is not transversed in the zero-th revolution
	d_min = std::min(R0,R1);
    }
    return;
}
}
#endif // KEP_TOOLBOX_CLOSEST_DISTANCE_H
