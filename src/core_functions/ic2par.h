/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://keptoolbox.sourceforge.net/index.html                            *
 *   http://keptoolbox.sourceforge.net/credits.html                          *
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


namespace kep_toolbox {
template<class vettore3D, class vettore6D>
void ic2par(const vettore3D& r0, const vettore3D& v0, const double &mu, vettore6D& E)
{
    double k[3];
    double h[3];
    double Dum_Vec[3];
    double n[3];
    double evett[3];

    double p = 0.0;
    double temp =0.0;
    double R0, ni;
    int i;

    ///1 - We compute h: the orbital angular momentum vector
    h[0] = r0[1]*v0[2] - r0[2]*v0[1];
    h[1] = r0[2]*v0[0] - r0[0]*v0[2];
    h[2] = r0[0]*v0[1] - r0[1]*v0[0];

    ///2 - We compute p: the orbital parameter
    p = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    p /= mu;

    //3 - We compute n: the vector of the node line
    //This operation is singular when inclination is zero, in which case the orbital parameters
    //are not defined
    k[0] = 0; k[1] = 0; k[2] = 1;
    n[0] = k[1]*h[2] - k[2]*h[1];
    n[1] = k[2]*h[0] - k[0]*h[2];
    n[2] = k[0]*h[1] - k[1]*h[0];
    temp = sqrt(n[0]*n[0] + n[1]*n[1] +n[2]*n[2]);
    for (i=0; i<3; i++)
        n[i] /= temp;

    //4 - We compute evett: the eccentricity vector
    //This operation is singular when eccentricity is zero, in which case the orbital parameters
    //are not defined
    R0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] +r0[2]*r0[2]);
    Dum_Vec[0] = v0[1]*h[2] - v0[2]*h[1];
    Dum_Vec[1] = v0[2]*h[0] - v0[0]*h[2];
    Dum_Vec[2] = v0[0]*h[1] - v0[1]*h[0];
    for (i=0; i<3; i++)
        evett[i] = Dum_Vec[i]/mu - r0[i]/R0;

    //The eccentricity is calculated and stored as the second orbital element
    E[1]  = sqrt(evett[0]*evett[0] + evett[1]*evett[1] +evett[2]*evett[2]);

    //The semi-major axis is calculated and stored as the first orbital element
    E[0] = p/(1-E[1]*E[1]);

    //Inclination is calculated and stored as the third orbital element
    E[2] = acos( h[2]/sqrt(h[0]*h[0] + h[1]*h[1] +h[2]*h[2]) );

    //Argument of pericentrum is calculated and stored as the fifth orbital element
    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=n[i]*evett[i];
    E[4] = acos(temp/E[1]);
    if (evett[2] < 0) E[4] = 2*M_PI - E[4];

    //Argument of longitude is calculated and stored as the fourth orbital element
    E[3] = acos(n[0]);
    if (n[1] < 0) E[3] = 2*M_PI-E[3];

    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=evett[i]*r0[i];

    //4 - We compute ni: the true anomaly (in 0, 2*PI)
    ni = acos(temp/E[1]/R0);

    temp = 0.0;
    for (i=0; i<3; i++)
        temp+=r0[i]*v0[i];

    if (temp<0.0) ni = 2*M_PI - ni;

    //Eccentric anomaly or the gudermannian is calculated and stored as the sixth orbital element
    if (E[1]<1.0)
        E[5] = 2.0*atan(sqrt((1-E[1])/(1+E[1]))*tan(ni/2.0));  // algebraic kepler's equation
    else
        E[5] =2.0*atan(sqrt((E[1]-1)/(E[1]+1))*tan(ni/2.0));   // algebraic equivalent of kepler's equation in terms of the Gudermannian
    return;
}
}

#endif // KEP_TOOLBOX_IC2PAR_H
