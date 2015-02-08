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

#ifndef KEPLERIAN_TOOLBOX_LAMBERT_2D_H
#define KEPLERIAN_TOOLBOX_LAMBERT_2D_H

#include<boost/bind.hpp>
#include<cmath>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/tools/roots.hpp>


#include"../numerics/regula_falsi.h"
#include"../astro_constants.h"
#include"x2tof.h"
#include"../exceptions.h"
#define D_ALP 1.5

namespace kep_toolbox {

inline double tof_curve(const double& ix,const double& s,const double& c,const double& tof,const int &lw){
    return ( log( kep_toolbox::x2tof(exp(ix)-1,s,c,lw,0) )-log(tof) );
}

inline double tof_curve_multi_rev(const double& ix,const double& s,const double& c,const double& tof,const int &lw, const int &N){
    return ( log( kep_toolbox::x2tof( (exp(ix*D_ALP)-1.0) / (exp(ix*D_ALP)+1.0),s,c,lw,N) ) - log(tof) );
}



/// Lambert solver (2 dimensional)
/**
 * This function solves a Lambert problem in its 'minimal' two-dimensional formulation.
 * It makes use of the Battin's variable \f$x$\f and the 'Izzo' plane to rectify
 * the time of flight curves (i.e. [log(x+1),log(tof)] and [log[(x+1)/(1-x)],log(tof)] for the mutirev case)
 *
 * \param[out] vr1 radial component of the velocity at r1
 * \param[out] vt1 tangential component of the velocity at r1
 * \param[out] vr2 radial component of the velocity at r2
 * \param[out] vt2 tangential component of the velocity at r2
 * \param[out] a semi major axis of the solution (negative is hyperbolae)
 * \param[out] p parameter of the solution (p = a * (1-e^2))
 *
 * \param[in] s semi perimeter of the triangle formed by r1=1,r2
 * \param[in] c chord joining r1=1 and r2
 * \param[in] tof time of flight in units R=r1, MU=1.
 * \param[in] lw when 1 the transfer with theta > pi is selected
 * \param[in] N number of revolutions (no multirev is default)
 * \param[in] branch selects the right or left branch of the tof curve when N>0 (default is the left branch)
 *
 * \return number of iterations to solve the tof equation. If 50, regula falsi algorithm has not converged
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

inline int lambert_2d(double &vr1, double &vt1, double &vr2, double &vt2, double &a, double &p,
            const double &s,const double &c, const double &tof, const int &lw, const int &N=0, const char &branch = 'l') {

    //Sanity checks
    if (tof <= 0) {
        throw_value_error("Time of flight needs to be positive");
    }

    if (c > s) {
        throw_value_error("The chord needs to be smaller than the semiperimeter");
    }
    
    if (c < s - 1) {
        throw_value_error("The chord needs to be larger than s - 1");
    }
    
    if (s < 1) {
        throw_value_error("The semiperimeter st be larger than 1 as r1 is used as unit");
    }

    if ( branch!='l' && branch!='r' ) {
        throw_value_error("Select either 'r' or 'l' branch for multiple revolutions");
    }

    //0 - Some geometry
    const double am = s/2.0;					//semi-major axis of the minimum energy ellipse
    const double r2 = 2.0 * s - c - 1.0;				//r2 in r1 units
    double tmp = acos ( (1.0 - c * c) / r2 / 2.0 + r2 / 2.0 );
    const double theta = (lw ? 2*M_PI-tmp : tmp);			//transfer angle
    const double lambda = sqrt (r2) * cos (theta/2.0) / s;

    int retval;

    //1 - We solve the tof equation
    double x;
    if (N==0) { //no multi-rev
        double ia = -0.69314718056; //log(1 - .5);
        double ib = 0.4054651081;   //log(1 + .5);
        retval = regula_falsi(ia,ib, boost::bind(kep_toolbox::tof_curve,_1,s,c,tof,lw), ASTRO_MAX_ITER,1e-12);

        x = exp(ia)-1.0;
    }
    else {	//multiple revolutions solution
        //left branch by default
	double ia = -1.0;
	double ib = -0.5;
	if (branch=='r') {
		ia = 1.0;
		ib = 0.5;
	}
        retval = regula_falsi(ia,ib, boost::bind(kep_toolbox::tof_curve_multi_rev,_1,s,c,tof,lw,N), ASTRO_MAX_ITER,1e-12);
        x = (exp(ia*D_ALP)-1.0) / (exp(ia*D_ALP)+1.0);
    }

    //3 - Using the Battin variable we recover all our outputs
    a = am/(1.0 - x*x);

    double beta,alfa,eta2,eta,psi,sigma1;

    if (x < 1.0)	// ellipse
    {
        beta = 2.0 * asin (sqrt( (s-c)/(2.0*a) ));
        if (lw) beta = -beta;
        alfa=2.0*acos(x);
        psi=(alfa-beta)/2.0;
        eta2=2.0*a*pow(sin(psi),2.0)/s;
        eta=sqrt(eta2);
    }
    else		// hyperbola
    {
        beta = 2.0*boost::math::asinh(sqrt((c-s)/(2.0*a)));
        if (lw) beta = -beta;
        alfa = 2.0*boost::math::acosh(x);
        psi = (alfa-beta)/2.0;
        eta2 = -2.0 * a * pow(sinh(psi),2.0)/s;
        eta = sqrt(eta2);
    }

    p = ( r2 / (am * eta2) ) * pow (sin (theta/2.0),2.0);
    sigma1 = (1.0/(eta * sqrt(am)) )* (2.0 * lambda * am - (lambda + x * eta));
    vr1 = sigma1;
    vt1 = sqrt(p);
    vt2 = vt1 / r2;
    vr2 = -vr1 + (vt1 - vt2)/tan(theta/2.0);
    return retval;
}
} //namespaces

#endif // KEPLERIAN_TOOLBOX_LAMBERT_2D_H
