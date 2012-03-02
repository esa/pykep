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


#ifndef PROPAGATE_TAYLOR_S_H
#define PROPAGATE_TAYLOR_S_H

#include<algorithm>
#include<cmath>
#include<boost/array.hpp>

#include"../exceptions.h"

namespace kep_toolbox {

template<class T>
double propagate_taylor_s_step(T& r0, T& v0, double &m0, double &t0 , const double &sf, const int &order, const T &thrust, const double &mu, const double &sundmann_alpha, const double &sundmann_c, const double &veff, const double &xm, const double &eps_a, const double &eps_r, std::vector< boost::array<double,8> > &x, std::vector< boost::array<double,25> > &u){

    double sqrtT = sqrt(thrust[0]*thrust[0] + thrust[1]*thrust[1] + thrust[2]*thrust[2]);

    //We initialize the initial conditions
    x[0][0] = r0[0];
    x[0][1] = r0[1];
    x[0][2] = r0[2];
    x[0][3] = v0[0];
    x[0][4] = v0[1];
    x[0][5] = v0[2];
    x[0][6] = m0;
    x[0][7] = t0;

    //We compute all needed Taylor coefficients
    int n=0;
    double gamma = sundmann_alpha/2;        //Exponent for u[12]
    double sigma = (sundmann_alpha-3.0)/2; 	//Exponent for u[13]
    while (n<order) {
        u[n][0] = x[n][0];  	//x
        u[n][1] = x[n][1];  	//y
        u[n][2] = x[n][2];  	//z
        u[n][3] = x[n][3];  	//vx
        u[n][4] = x[n][4];  	//vy
        u[n][5] = x[n][5];  	//vz
        u[n][6] = x[n][6];  	//m
        u[n][7] = x[n][7];  	//time

        for (int j=0;j<=n;j++) u[n][8]  += u[j][0]*u[n-j][0];  	//x^2
        for (int j=0;j<=n;j++) u[n][9]  += u[j][1]*u[n-j][1];  	//y^2
        for (int j=0;j<=n;j++) u[n][10] += u[j][2]*u[n-j][2]; 	//z^2
        u[n][11] = u[n][8] + u[n][9] + u[n][10];  		//r^2


        if (n==0){
            u[n][12] = pow(u[n][11],gamma);
        } else {
            for (int j=0;j<n;++j) u[n][12] += (gamma*n - j*(gamma+1))*u[n-j][11]*u[j][12];
            u[n][12] = u[n][12] / n / u[0][11];
        } //r^alpha

        if (n==0){
            u[n][13]= pow(u[n][11],sigma);
        } else {
            for (int j=0;j<n;++j) u[n][13] += (sigma*n - j*(sigma+1))*u[n-j][11]*u[j][13];
            u[n][13] = u[n][13] / n / u[0][11];
        } //1 / r^(3-alpha)

        for (int j=0;j<=n;j++) u[n][14] += u[j][13]*u[n-j][0];	// 1 / r^(3-alpha) x
        for (int j=0;j<=n;j++) u[n][15] += u[j][13]*u[n-j][1];	// 1 / r^(3-alpha) y
        for (int j=0;j<=n;j++) u[n][16] += u[j][13]*u[n-j][2];	// 1 / r^(3-alpha) z

        if (n==0){
            u[n][17] = u[n][12] / u[n][6];
        } else {
            for (int j=1;j<=n;j++) u[n][17] += u[j][6] * u[n-j][17];
            u[n][17] = 1.0 / u[0][6] * (u[n][12] - u[n][17]);
        } // r^alpha/m

        for (int j=0;j<=n;j++) u[n][18] += u[j][3]*u[n-j][12];				// eq1
        for (int j=0;j<=n;j++) u[n][19] += u[j][4]*u[n-j][12];				// eq2
        for (int j=0;j<=n;j++) u[n][20] += u[j][5]*u[n-j][12];				// eq3
        u[n][21] = - mu * u[n][14] + u[n][17] * thrust[0];  				// eq4
        u[n][22] = - mu * u[n][15] + u[n][17] * thrust[1];  				// eq5
        u[n][23] = - mu * u[n][16] + u[n][17] * thrust[2];  				// eq6
        u[n][24] = - sqrtT / veff * u[n][12];                               // eq7

        x[n+1][0] = 1./(n+1) * sundmann_c * u[n][18];
        x[n+1][1] = 1./(n+1) * sundmann_c * u[n][19];
        x[n+1][2] = 1./(n+1) * sundmann_c * u[n][20];
        x[n+1][3] = 1./(n+1) * sundmann_c * u[n][21];
        x[n+1][4] = 1./(n+1) * sundmann_c * u[n][22];
        x[n+1][5] = 1./(n+1) * sundmann_c * u[n][23];
        x[n+1][6] = 1./(n+1) * sundmann_c * u[n][24];
        x[n+1][7] = 1./(n+1) * sundmann_c * u[n][12];
        n++;
    }


    //We now compute the Taylor expansion, first the optimal step size and then the sum....

    //Determining the optimal step size (see Jorba's method)
    double step,rho_m;

    //The infinity norm of the highest order derivative
    double xm_n = std::max(std::abs(x[n][0]),std::abs(x[n][1]));
    xm_n = std::max(xm_n,std::abs(x[n][2]));
    xm_n = std::max(xm_n,std::abs(x[n][3]));
    xm_n = std::max(xm_n,std::abs(x[n][4]));
    xm_n = std::max(xm_n,std::abs(x[n][5]));
    xm_n = std::max(xm_n,std::abs(x[n][6]));
    xm_n = std::max(xm_n,std::abs(x[n][7]));

    //The infinity norm of the n-1 order derivative
    double xm_n1 = std::max(std::abs(x[n-1][0]),std::abs(x[n-1][1]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][2]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][3]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][4]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][5]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][6]));
    xm_n1 = std::max(xm_n1,std::abs(x[n-1][7]));

    if (eps_r*xm < eps_a) {
        rho_m = std::min(pow((1/xm_n),1./n),pow((1/xm_n1),1./(n-1)));
    } else {
        rho_m = std::min(pow((xm/xm_n),1./n),pow((xm/xm_n1),1./(n-1)));
    }
    step = rho_m/(M_E*M_E);

    if (sf<0) step = -step;
    if (std::abs(step) > std::abs(sf)) step = sf;

    //Horner method could be better here ...
    double steppow = step;
    for(int j=1; j<=order;++j){
        r0[0] += x[j][0]*steppow;
        r0[1] += x[j][1]*steppow;
        r0[2] += x[j][2]*steppow;
        v0[0] += x[j][3]*steppow;
        v0[1] += x[j][4]*steppow;
        v0[2] += x[j][5]*steppow;
        m0    += x[j][6]*steppow;
        t0    += x[j][7]*steppow;
        steppow*=step;
    }
    return step;
}

/// Taylor series propagation of a constant thrust arc using the Generalized Sundmann Transformation
/**
 * This template function propagates an initial state using the generalized Sundmann-Transformation. The independent variable thus is not
 * time, but \f$ ds = c r^\alpha dt \f$
 *
 * \param[in,out]	r0 initial position vector. On output contains the propagated position. (r0[1],r0[2],r0[3] need to be preallocated, suggested template type is boost::array<double>,3))
 * \param[in,out]	v0 initial velocity vector. On output contains the propagated velocity. (v0[1],v0[2],v0[3] need to be preallocated, suggested template type is boost::array<double>,3))
 * \param[in,out]	m0 initial mass. On output contains the propagated value
 * \param[in,out]	t0 initial time. on output contains the propagated time
 * \param[in]		u  thrust vector
 * \param[in,out]	sf propagation pesudo-time (can be negative). If the maximum number of iterations is reached, the pseudo-time is returned where the state is calculated for the last time
 * \param[in]		mu central body gravitational parameter
 * \param[in]		veff the product g0*Isp characteristic of the engine
 * \param[in]		c constant in the generalized Sundmann transform defaults to 1.0
 * \param[in]		alpha exponent in the generalized Sundmann transform (defaults to 1.5)
 * \param[in]		log10tolerance logarithm of the desired absolute tolerance
 * \param[in]		log10rtolerance logarithm of the desired relative tolerance
 * \param[in]		max_iter maximum number of iteration allowed
 * \param[in]		max_order maximum order for the polynomial expansion
 *
 * \throw value_error if max_iter is hit.....
 * \throw value_error if max_order is exceeded.....
 *
 * NOTE: Equations of motions are written and propagated in cartesian coordinates
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
template<class T>
void propagate_taylor_s(T& r0, T& v0, double &m0, double &t0, const T& thrust, const double &sf, const double &mu = 1, const double &veff = 1, const double &c = 1.0, const double &alpha = 1.5, const int &log10tolerance=-10, const int &log10rtolerance=-10, const int &max_iter = 10000, const int &max_order = 3000){

    boost::array<double,8> dumb;
    boost::array<double,25> dumb2;
    for (int i=0;i<8;++i) dumb[i]=0;
    for (int i=0;i<25;++i) dumb2[i]=0;

    std::vector< boost::array<double,8> > _x;   // x[order][var]
    std::vector< boost::array<double,25> > _u;   // u[order][var]

    double step = sf;
    double eps_a = pow(10.,log10tolerance);
    double eps_r = pow(10.,log10rtolerance);
    double eps_m,xm;
    int j;
    for (j=0; j< max_iter; ++j) {
        //We follow the method described by Jorba in "A software package ...."
        //1 - We determine eps_m from Eq. (7)
        xm = std::max(std::abs(r0[0]),std::abs(r0[1]));
        xm = std::max(xm,std::abs(r0[2]));
        xm = std::max(xm,std::abs(v0[0]));
        xm = std::max(xm,std::abs(v0[1]));
        xm = std::max(xm,std::abs(v0[2]));
        xm = std::max(xm,std::abs(m0));
        xm = std::max(xm,std::abs(t0));

        //2 - We evaluate the polynomial order
        (eps_r*xm < eps_a) ? eps_m = eps_a : eps_m = eps_r;
        int order = (int) ( ceil(-0.5*log(eps_m) + 1) );
        if (order > max_order) throw_value_error("Polynomial order is too high.....");

        //3 - We allocate or deallocate memory if necessary
        _x.assign(order+1,dumb);
        _u.assign(order,dumb2);
        double h = propagate_taylor_s_step(r0,v0,m0,t0,step,order,thrust,mu,alpha,c,veff,xm, eps_a, eps_r,_x,_u);
        if (std::abs(h)>=std::abs(step)) break; else {
            step = step - h;
        }
    }
    if (j>max_iter-1) throw_value_error("Maximum number of iteration reached in Taylor integration (sundmann)");
}

} //Namespace kep_toolbox

#endif // PROPAGATE_TAYLOR_S_H
