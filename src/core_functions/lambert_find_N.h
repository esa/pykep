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

#ifndef KEPLERIAN_TOOLBOX_LAMBERT_FIND_N_H
#define KEPLERIAN_TOOLBOX_LAMBERT_FIND_N_H

#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>

#include "x2tof.h"


namespace kep_toolbox {

/// Finds multi rev. number in Lambert's problem
/**
 * The function returns the maximum number of revolutions for a particular Lambert's problem to
 * admit a solution. It tries to avoid the calculation of the minima of the relevant tof curve by
 * pruning out obvious cases
 *
 * \param[in] s semi perimeter of the triangle formed by r1,r2
 * \param[in] c chord
 * \param[in] tof time of flight. Units such as MU=1.
 * \param[in] lw when 1 the transfer with theta > pi is selected
 *
 * \return Maximum number of revolutions allowed
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

inline int lambert_find_N(const double &s,const double &c, const double &tof, const int &lw) {
    double Tm = M_PI/2 * sqrt(2*s*s*s);		//Minimum energy ellipse period
    int Ntemp = tof / Tm;				//It is either Nmax=Ntemp or Nmax = Ntemp-1
    if (Ntemp == 0) return 0;
    double Tmax = x2tof(0,s,c,lw,Ntemp);
    if (tof > Tmax) return Ntemp;
    std::pair<double,double> res = boost::math::tools::brent_find_minima(boost::bind(x2tof,_1,s,c,lw,Ntemp),0.0,0.5,8);
    if (res.second < tof) return Ntemp;
    return Ntemp - 1;
}

} //namespaces

#endif // KEPLERIAN_TOOLBOX_LAMBERT_FIND_N_H
