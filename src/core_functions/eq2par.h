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


#ifndef KEP_TOOLBOX_EQ2PAR_H
#define KEP_TOOLBOX_EQ2PAR_H

#include <cmath>

#include "../core_functions/convert_anomalies.h"
#include <boost/math/constants/constants.hpp>


namespace kep_toolbox {

/// From modified equinoctial elements to osculating parameters
/**
* Transforms modified equinoctial elements (p,h,k,p,q,L) to osculating Keplerian parameters (a,e,i,W,w,E). 
* Note that we use the eccentric anomaly E (or gudermannian) and the true longitude L (i.e. f + om + Om).
*
* Note that nan will be generated if the osculating parameters are not defined (e=0, i=0, etc..)
*
* @param[in] EQ the equinoctial elements (p,h,k,p,q,L), anomaly in rad.
* @param[in] retrogade forces retrograde elements to be used
* 
* @param[out] E the osculating Keplerian parameters (a,e,i,W,w,E).
*
*/
template<class vettore6D>
void eq2par(vettore6D& E, const vettore6D& EQ, const bool retrogade = false)
{
    int I = 1;
    if (retrogade) {
        I = -1;
    }
    auto ecc = std::sqrt(EQ[1]*EQ[1]+EQ[2]*EQ[2]);
    auto tmp = std::sqrt(EQ[3]*EQ[3]+EQ[4]*EQ[4]);
    auto zita = std::atan2(EQ[1] / ecc, EQ[2] / ecc);
    
    E[1] = ecc;
    E[0] = EQ[0] / (1-ecc*ecc);
    
    E[2] = boost::math::constants::pi<double>() / 2 * (1-I) + 2 * I * std::atan(tmp);
    E[3] = std::atan2(EQ[3] / tmp, EQ[4] / tmp);
    E[4] = zita - I * E[3];
    E[5] = EQ[5] - zita;
    if (E[1] < 1) {
        E[5] = f2e(E[5], E[1]);
    } else { // E[5] is the Gudermannian
        E[5] = f2zeta(E[5], E[1]);
    }
    return;
}
}
#endif // KEP_TOOLBOX_PAR2IC_H
