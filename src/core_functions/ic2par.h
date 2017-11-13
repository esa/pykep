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

#ifndef KEP_TOOLBOX_IC2PAR_H
#define KEP_TOOLBOX_IC2PAR_H

#include <cmath>

#include "array3D_operations.h"

namespace kep_toolbox {
template<class vettore3D, class vettore6D>
void ic2par(const vettore3D& r0, const vettore3D& v0, const double &mu, vettore6D& E)
{
    vettore3D k = { 0.0, 0.0, 1.0 };
    // build generic arrays to size - init values don't matter:
    vettore3D h = { 0.0, 0.0, 0.0 };
    vettore3D Dum_Vec = { 0.0, 0.0, 0.0 };
    vettore3D n = { 0.0, 0.0, 0.0 };
    vettore3D evett = { 0.0, 0.0, 0.0 };

    double p = 0.0;
    double temp = 0.0;
    double R0, ni;
    int i;

    ///1 - We compute h: the orbital angular momentum vector
	cross( h, r0, v0 );

    ///2 - We compute p: the orbital parameter
	p = dot( h, h ) / mu; // h^2 / mu

    ///3 - We compute n: the vector of the node line
    ///This operation is singular when inclination is zero, in which case the orbital parameters
    ///are not defined
	cross( n, k, h );
	vers( n, n ); // vers(x, y) = unit vector of y -> x

    ///4 - We compute evett: the eccentricity vector
    ///This operation is singular when eccentricity is zero, in which case the orbital parameters
    ///are not defined
	R0 = norm( r0 );
	cross( Dum_Vec, v0, h );
    for (i=0; i<3; i++)
        evett[i] = Dum_Vec[i]/mu - r0[i]/R0;

    ///The eccentricity is calculated and stored as the second orbital element
	E[1] = norm( evett );

    ///The semi-major axis is calculated and stored as the first orbital element
    E[0] = p/(1-E[1]*E[1]);

    ///Inclination is calculated and stored as the third orbital element
    E[2] = acos( h[2]/norm(h) );

    ///Argument of pericentrum is calculated and stored as the fifth orbital element
	temp = dot( n, evett );
    E[4] = acos(temp/E[1]);
    if (evett[2] < 0) E[4] = 2*M_PI - E[4];

    ///Argument of longitude is calculated and stored as the fourth orbital element
    E[3] = acos(n[0]);
    if (n[1] < 0) E[3] = 2*M_PI-E[3];

	temp = dot( evett, r0 );
	
    ///4 - We compute ni: the true anomaly (in 0, 2*PI)
    ni = acos(temp/E[1]/R0);

	temp = dot( r0, v0 );

    if (temp<0.0) ni = 2*M_PI - ni;

    ///Eccentric anomaly or the gudermannian is calculated and stored as the sixth orbital element
    if (E[1]<1.0)
        E[5] = 2.0*atan(sqrt((1-E[1])/(1+E[1]))*tan(ni/2.0));  // algebraic kepler's equation
    else
        E[5] =2.0*atan(sqrt((E[1]-1)/(E[1]+1))*tan(ni/2.0));   // algebraic equivalent of kepler's equation in terms of the Gudermannian
    return;
}
}

#endif // KEP_TOOLBOX_IC2PAR_H
