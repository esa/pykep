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

#ifndef KEPLERIAN_TOOLBOX_LAMBERT_3D_H
#define KEPLERIAN_TOOLBOX_LAMBERT_3D_H

#include<cmath>

#include"lambert_2d.h"
#include"../astro_constants.h"
#include"array3D_operations.h"

namespace kep_toolbox {


/// Lambert solver (3 dimensional)
/**
 * This function solves a Lambert problem embedded into a 3D space. It thus
 * will be singular for 180 degrees problems as the transfer plane is, in that case, undetermined.
 *
 * \param[out] v1 velocity vector at r1 (preallocated)
 * \param[out] v2 velocity vector at r2 (preallocated)
 * \param[out] a semi major axis of the solution (negative is hyperbolae)
 * \param[out] p parameter of the solution (p = a * (1-e^2))
 *
 *\param[in] r1 first position vector
 *\param[in] r2 second position vector
 *\param[in] tof time of flight
 *\param[in] mu gravity parameter
 *\param[in] lw if == 1 long way is selected
 *\param[in] N number of multiple revolutions (default is 0)
 *\param[in] branch selects the tof branch in case N>0 (default is left branch)
 *
 * \return number of iterations to solve the tof equation. If 50, regula falsi algorithm has not converged
 *
 *NOTE 1: when N>0 there may just not be any solution to the lambert problem in which
 *case this function does not detect it. Reaching a maximum number of iteration in the regula falsi
 *method is an indicator that no solution exists.
 *
 *NOTE 2: the function works in all units as long as fed with consistent ones
 *
 * @see kep_toolbox::lambert_2d
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
template <class vettore3D>
inline int lambert_3d (vettore3D &v1, vettore3D &v2, double &a, double &p,
            const vettore3D &r1, const vettore3D &r2, const double &tof, const double &mu,
            const int &lw, const int &N=0, const char branch='l') {

    // 0 - Sanity checks
    if (tof <= 0)
    {
        throw_value_error("time of flight in lambert_3d is negative");
    }
    if (mu <= 0)
    {
        throw_value_error("gravity parameter in lambert_3d is negative");
    }

    if ( branch!='l' && branch!='r' ) {
        throw_value_error("Select either 'r' or 'l' branch for multiple revolutions");
    }

    // 1 - Computing non dimensional units
    double R = norm(r1);
    double V = sqrt(mu / R);
    double T = R/V;

    // 2 - Computing geometry of transfer
    double R2 = norm(r2);
    double costheta = dot(r1,r2);
    costheta /= R*R2;

    double r2_mod = R2 / R;
    double c = sqrt(1 + r2_mod*(r2_mod - 2.0 * costheta));
    double s = (1 + r2_mod + c)/2.0;

    // 3 - Solving the problem in 2D
    double vr1,vr2,vt1,vt2;
    int retval;
    retval = lambert_2d(vr1,vt1,vr2,vt2,a,p,s,c,tof / T,lw,N,branch);

    // 4 - Reconstructing the velocity vectors in three dimensions
    array3D ir1,it1,ir2,it2,ih;
    vers(ir1,r1);
    vers(ir2,r2);
    (lw ? cross(ih,ir2,ir1) : cross(ih,ir1,ir2));     //here is the singularity: as when ir1||ir2 this plane is not defined!!
	if (norm(ih) == 0) {
        throw_value_error("lambert problem is singular in 3D as the transfer angle is 180*n degrees");
    }
    vers(ih,ih);
    cross(it1,ir1,ih);
    cross(it2,ir2,ih);
    for (int i=0;i<3;++i) v1[i] = vr1*ir1[i] - vt1*it1[i];
    for (int i=0;i<3;++i) v2[i] = vr2*ir2[i] - vt2*it2[i];

    // 5 - Putting back dimensions
    for (int i=0;i<3;++i) v1[i] *= V;
    for (int i=0;i<3;++i) v2[i] *= V;
    a *= R;
    p *= R;
    return retval;
}
} //namespaces

#endif // KEPLERIAN_TOOLBOX_LAMBERT_3D_H
