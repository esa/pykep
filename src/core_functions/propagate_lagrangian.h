/*****************************************************************************
 *   Copyright (C) 2004-2014 The PyKEP development team,                     *
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


#ifndef PROPAGATE_LAGRANGIAN_H
#define PROPAGATE_LAGRANGIAN_H

#include<boost/bind.hpp>
#include<boost/math/tools/roots.hpp>

#include"../astro_constants.h"
#include"../numerics/newton_raphson.h"
#include"kepler_equations.h"



namespace kep_toolbox {

/// Lagrangian propagation
/**
 * This template function propagates an initial state for a time t assuming a central body and a keplerian
 * motion. Lagrange coefficients are used as basic numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 *
 * \param[in,out] r0 initial position vector. On output contains the propagated position. (r0[1],r0[2],r0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in,out] v0 initial velocity vector. On output contains the propagated velocity. (v0[1],v0[2],v0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in] t propagation time (can be negative)
 * \param[in] mu central body gravitational parameter
 *
 * NOTE: The solver used for the kepler equation is a derivative free solver from the boost libraries.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
template<class T>
void propagate_lagrangian(T& r0, T& v0, const double &t, const double &mu)
{
    double R = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
    double V = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
    double energy = (V*V/2 - mu/R);
    double a = - mu / 2.0 / energy;
    double sqrta;
    double F,G,Ft,Gt;

    double sigma0 = (r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2]) / sqrt(mu);

    if (a > 0){	//Solve Kepler's equation, elliptical case
        sqrta = sqrt(a);
        double DM = sqrt(mu / pow(a,3)) * t;
        double DE = DM;

        //Solve Kepler Equation for ellipses in DE (eccentric anomaly difference)
        //newton_raphson(DE,boost::bind(kepDE,_1,DM,sigma0,sqrta,a,R),boost::bind(d_kepDE,_1,sigma0,sqrta,a,R),100,ASTRO_TOLERANCE);
        std::pair<double, double> result;
        boost::uintmax_t iter = ASTRO_MAX_ITER;
        boost::math::tools::eps_tolerance<double> tol(64);
        result = boost::math::tools::bracket_and_solve_root(boost::bind(kepDE,_1,DM,sigma0,sqrta,a,R),DE,2.0,true,tol,iter);
        DE = (result.first + result.second) / 2;
        double r = a + (R - a) * cos(DE) + sigma0 * sqrta * sin(DE);

        //Lagrange coefficients
        F  = 1 - a / R * (1 - cos(DE));
        G  = a * sigma0 / sqrt(mu) * (1 - cos(DE)) + R * sqrt(a / mu) * sin(DE);
        Ft = -sqrt(mu * a) / (r * R) * sin(DE);
        Gt = 1 - a / r * (1 - cos(DE));
    }
    else{	//Solve Kepler's equation, hyperbolic case
        sqrta = sqrt(-a);
        double DN = sqrt(-mu / pow(a,3)) * t;
        double DH;
        t > 0 ? DH = 1 : DH = -1; // TODO: find a better initial guess. I tried with 0 and D (both have numercial problems and result in exceptions)

        //Solve Kepler Equation for hyperbolae in DH (hyperbolic anomaly difference)
        //newton_raphson(DH,boost::bind(kepDH,_1,DN,sigma0,sqrta,a,R),boost::bind(d_kepDH,_1,sigma0,sqrta,a,R),100,ASTRO_TOLERANCE);
        std::pair<double, double> result;
        boost::uintmax_t iter = ASTRO_MAX_ITER;
        boost::math::tools::eps_tolerance<double> tol(64);
        result = boost::math::tools::bracket_and_solve_root(boost::bind(kepDH,_1,DN,sigma0,sqrta,a,R),DH,2.0,true,tol,iter);
        DH = (result.first + result.second) / 2;
        double r = a + (R - a) * cosh(DH) + sigma0 * sqrta * sinh(DH);

        //Lagrange coefficients
        F  = 1 - a / R * (1 - cosh(DH));
        G  = a * sigma0 / sqrt(mu) * (1 - cosh(DH)) + R * sqrt(-a / mu) * sinh(DH);
        Ft = -sqrt(-mu * a) / (r * R) * sinh(DH);
        Gt = 1 - a / r * (1 - cosh(DH));
    }

    double temp[3] = {r0[0],r0[1],r0[2]};
    for (int i=0;i<3;i++){
        r0[i] = F * r0[i] + G * v0[i];
        v0[i] = Ft * temp[i] + Gt * v0[i];
    }
}
}

#endif // PROPAGATE_LAGRANGIAN_H
