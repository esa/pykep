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

#ifndef KEP_TOOLBOX_PAR2EQ_H
#define KEP_TOOLBOX_PAR2EQ_H

#include <cmath>

#include "../core_functions/convert_anomalies.h"

namespace kep_toolbox
{

/// From osculating Keplerian to equinoctial
/**
* Transforms osculating Keplerian parameters (a,e,i,W,w,E) to
* modified equinoctial elements (p,f,g,h,k,L).
* Note that we use the eccentric anomaly E (or Gudermannian) and the true longitude L (i.e. f + om + Om)
*
* @param[in] E the osculating Keplerian parameters (a,e,i,W,w,E), anomalies in rad.
* @param[in] retrogade forces retrograde elements to be used
*
* @param[out] EQ the equinoctial elements (a,h,k,p,q,L).
*
*/
template <class vettore6D>
void par2eq(vettore6D &EQ, const vettore6D &E, const bool retrogade = false)
{
    int I;
    if (retrogade) {
        I = -1;
        EQ[3] = 1. / std::tan(E[2] / 2) * std::cos(E[3]);
        EQ[4] = 1. / std::tan(E[2] / 2) * std::sin(E[3]);
    } else {
        I = 1;
        EQ[3] = std::tan(E[2] / 2) * std::cos(E[3]);
        EQ[4] = std::tan(E[2] / 2) * std::sin(E[3]);
    }
    EQ[0] = E[0] * (1 - E[1] * E[1]);
    EQ[1] = E[1] * std::cos(E[4] + I * E[3]);
    EQ[2] = E[1] * std::sin(E[4] + I * E[3]);
    double f;
    if (E[1] < 1) {
        f = e2f(E[5], E[1]);
    } else { // E[5] is the Gudermannian
        f = zeta2f(E[5], E[1]);
    }
    EQ[5] = f + E[4] + I * E[3];
    return;
}
}
#endif // KEP_TOOLBOX_PAR2IC_H
