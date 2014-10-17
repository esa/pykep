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


#ifndef X2TOF_H
#define X2TOF_H

#include<cmath>

#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>

namespace kep_toolbox {
/// Convert the Battin's variable x to time of flight in a non dimensional Lambert's problem
/**
 * A Lambert problem is defined by \f$ r_1,r_2,t\f$ and \f$ \mu\f$. Assuming as length units \f$ |r_1|\f$ and as velocity
 * units \f$ \sqrt(\mu/|r_1|)\f$ (i.e. \f$ \mu = 1\f$) the Lambert's problem is only defined by \f$ r_2,t\f$. Employing these units
 * this function returns the time-of flight from the Battin's variable \f$ x\f$.
 *
 * \param[in] x Battin's variable
 * \param[in] s half perimeter of the triangle defined by r1,r2
 * \param[in] c chord defined by the r1,r2 triangle
 * \param[in] lw switch between long and short way solutions (when lw=1 the long way solution sselected)
 * \param[in] N multiple revolutions (default is 0)
 *
 * \see Battin "Introduction to the Mahematics and Methods of Astrodynamic"
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
inline double x2tof(const double &x,const double &s,const double &c,const int &lw, const int &N = 0)
{
    double am,a,alfa,beta;

    am = s/2;
    a = am/(1-x*x);
    if (x < 1)	//ellipse
	{
        beta = 2 * asin (sqrt((s - c)/(2*a)));
        if (lw) beta = -beta;
		alfa = 2 * acos(x);
    }
    else
    {
        alfa = 2 * boost::math::acosh(x);
        beta = 2 * boost::math::asinh(sqrt ((s - c)/(-2 * a)));
        if (lw) beta = -beta;
    }

    if (a > 0)
    {
        return (a * sqrt (a)* ( (alfa - sin(alfa)) - (beta - sin(beta)) + 2*M_PI*N));
    }
    else
    {
        return (-a * sqrt(-a)*( (sinh(alfa) - alfa) - ( sinh(beta) - beta)) );
    }

}

} //namespace

#endif // X2TOF_H
