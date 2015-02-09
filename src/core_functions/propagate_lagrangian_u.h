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


#ifndef KEP_TOOLBOX_PROPAGATE_LAGRANGIAN_U_H
#define KEP_TOOLBOX_PROPAGATE_LAGRANGIAN_U_H

#include <boost/bind.hpp>
#include <boost/math/tools/roots.hpp>

#include "../astro_constants.h"
#include "../numerics/newton_raphson.h"
#include "kepler_equations.h"
#include "stumpff.h"



namespace kep_toolbox {

/// Lagrangian propagation using the universal anomaly
/**
 * This template function propagates an initial state for a time t assuming a central body and a keplerian
 * motion. Lagrange coefficients are used as basic numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 *
 * \param[in,out] r0 initial position vector. On output contains the propagated position. (r0[1],r0[2],r0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in,out] v0 initial velocity vector. On output contains the propagated velocity. (v0[1],v0[2],v0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in] t propagation time
 * \param[in] mu central body gravitational parameter
 *
 * NOTE: Negative times are dealt by inverting the time sign and the initial conditions
 *
 * @see http://www.google.it/url?sa=t&source=web&cd=1&ved=0CBYQFjAA&url=http%3A%2F%2Fwww3.uta.edu%2Ffaculty%2Fsubbarao%2FMAE3304Astronautics%2FSampleStuff%2Fappend-d.pdf&ei=8eL0TKDUKMrrOcj2ybMI&usg=AFQjCNFLBgLMvPWSDsCvZMVOW3kJV9uh-Q
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
template<class T>
void propagate_lagrangian_u(T& r0, T& v0, const double &t, const double &mu = 1)
{
    //If time is negative we need to invert time and velocities. Unlike the other formulation
    //of the propagate lagrangian we cannot rely on negative times to automatically mean back-propagation
    double t_copy = t;
    if (t < 0)
    {
        t_copy = -t;
        v0[0] = -v0[0]; v0[1] = -v0[1]; v0[2] = -v0[2];
    }

    double F,G,Ft,Gt;
    double R0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
    double V0 = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
    //the reciprocal of the semi-major axis
    double alpha = 2/R0 - V0*V0/mu;
    //initial radial velocity
    double VR0 = (r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2]) / R0;

    //solve kepler's equation in universal variables
    double DS = 1;
    alpha > 0 ? DS = sqrt(mu)*t_copy*std::abs(alpha) : DS = 1; //initial guess for the universal anomaly. For hyperbolas it is 1.... can be better?
    //newton_raphson(DS,boost::bind(kepDS,_1,t_copy,R0,VR0,alpha,mu),boost::bind(d_kepDS,_1,R0,V0,alpha,mu),100,ASTRO_TOLERANCE);
    std::pair<double, double> result;
    boost::uintmax_t iter = ASTRO_MAX_ITER;
    boost::math::tools::eps_tolerance<double> tol(64);
    result = boost::math::tools::bracket_and_solve_root(boost::bind(kepDS,_1,t_copy,R0,VR0,alpha,mu),DS,2.0,true,tol,iter);
    DS = (result.first + result.second) / 2;

    //evaluate the lagrangian coefficients F and G
    double S = stumpff_s(alpha*DS*DS);
    double C = stumpff_c(alpha*DS*DS);
    //
    double z = alpha*DS*DS;
    F = 1 - DS*DS/R0*C;
    G = t_copy - 1/sqrt(mu)*DS*DS*DS*S;

    //compute the final position
    T rf;
    rf[0] = F*r0[0] + G*v0[0];
    rf[1] = F*r0[1] + G*v0[1];
    rf[2] = F*r0[2] + G*v0[2];
    double RF = sqrt(rf[0]*rf[0] + rf[1]*rf[1] + rf[2]*rf[2]);

    //compute the lagrangian coefficients Ft, Gt
    Ft = sqrt(mu)/RF/R0*(z*S - 1)*DS;
    Gt = 1 - DS*DS/RF*C;

    //compute the final velocity
    T vf;
    vf[0] = Ft*r0[0] + Gt*v0[0];
    vf[1] = Ft*r0[1] + Gt*v0[1];
    vf[2] = Ft*r0[2] + Gt*v0[2];

    r0=rf;
    v0=vf;
    if (t < 0)
    {
        v0[0] = -v0[0];v0[1] = -v0[1]; v0[2] = -v0[2];
    }
}
}

#endif // KEP_TOOLBOX_PROPAGATE_LAGRANGIAN_U_H
