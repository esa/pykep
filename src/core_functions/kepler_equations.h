/*****************************************************************************
 *   Copyright (C) 2004-2012 The PyKEP development team,                     *
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

#ifndef KEPLER_EQUATIONS_H
#define KEPLER_EQUATIONS_H

#include<cmath>
#include"stumpff.h"
#include"array3D_operations.h"

namespace kep_toolbox {
    //With the eccentric anomaly (E)
    inline double kepE(const double& E, const double& M, const double& eccentricity ){
        return ( E - eccentricity*sin(E) - M );
    }
    inline double d_kepE(const double& E, const double& eccentricity){
        return ( 1 - eccentricity*cos(E) );
    }

    //With the eccentric anomaly difference (DE)
    inline double kepDE(const double& DE, const double& DM, const double& sigma0, const double& sqrta, const double& a, const double& R){
        return ( (-DM + DE + sigma0 / sqrta * (1 - cos(DE)) - (1 - R / a) * sin(DE)) );
    }

    inline double d_kepDE(const double& DE, const double& sigma0, const double& sqrta, const double& a, const double& R){
        return ( (1 + sigma0 / sqrta * sin(DE) - (1 - R / a) * cos(DE)) );
    }

    //With the hyperbolic anomaly difference (DH)
    inline double kepDH(const double& DH, const double& DN, const double& sigma0, const double& sqrta, const double& a, const double& R){
        return ( -DN -DH + sigma0/sqrta * (cosh(DH) - 1) + (1 - R / a) * sinh(DH) );
    }

	//Barker time of flight equation
	inline double barker(const array3D& r1, const array3D& r2, const double mu){
		double R1 = norm(r1);
		double R2 = norm(r2);
		array3D r21;
		diff(r21,r2,r1);
		double R21 = norm(r21);
		double x = r1[0]*r2[1]-r1[1]*r2[0];
		double sigma;
		sigma = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
		return ( pow(R1+R2+R21,1.5) - sigma*pow(R1+R2-R21,1.5) )/(6.0*sqrt(mu));
	}

    inline double d_kepDH(const double& DH, const double& sigma0, const double& sqrta, const double& a, const double& R){
        return ( -1 + sigma0 / sqrta * sinh(DH) + (1 - R / a) * cosh(DH) );
    }
    //With the universal anomaly difference (DS)
    inline double kepDS(const double& DS, const double& DT, const double& r0, const double& vr0, const double& alpha, const double& mu){
        double S = stumpff_s(alpha*DS*DS);
        double C = stumpff_c(alpha*DS*DS);
        double retval = -sqrt(mu)*DT + r0*vr0*DS*DS*C/sqrt(mu) + (1-alpha*r0)*DS*DS*DS*S + r0*DS;
        return ( retval );
    }
    inline double d_kepDS(const double& DS, const double& r0, const double& vr0, const double& alpha, const double& mu){
        double S = stumpff_s(alpha*DS*DS);
        double C = stumpff_c(alpha*DS*DS);
        double retval = r0*vr0/sqrt(mu)*DS * (1-alpha*DS*DS*S) + (1-alpha*r0)*DS*DS*C + r0;
        return ( retval );
    }
}
#endif // KEPLER_EQUATIONS_H
