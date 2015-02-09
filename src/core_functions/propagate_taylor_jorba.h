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


#ifndef KEP_TOOLBOX_PROPAGATE_TAYLOR_JORBA_H
#define KEP_TOOLBOX_PROPAGATE_TAYLOR_JORBA_H

extern "C"{
    #include"jorba.h"
}

namespace kep_toolbox {

/// Taylor series propagation of a constant thrust segment
/**
 * This template function propagates an initial state for a time t assuming a central body and a keplerian
 * motion perturbed by an inertially constant thrust u
 *
 * \param[in,out] r0 initial position vector. On output contains the propagated position. (r0[1],r0[2],r0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in,out] v0 initial velocity vector. On output contains the propagated velocity. (v0[1],v0[2],v0[3] need to be preallocated, suggested template type is boost::array<double,3))
 * \param[in] T thrust vector (cartesian components)
 * \param[in] t propagation time (can be negative)
 * \param[in] mu central body gravitational parameter
 *
 * NOTE: The Taylor propagation was genrated using Jorba's tool available on-line
 *
 * @see http://www.maia.ub.es/~angel/taylor
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
template<class T>
void propagate_taylor_jorba(T& r0, T& v0, double &m0, const T& u, const double &t, const double &mu, const double &veff, const int &log10tolerance=-9, const int &log10rtolerance=-9){
    int i, order=20, itmp=0, direction;
    MY_FLOAT  startT, stopT, nextT;
    MY_FLOAT  xx[7];
    InitMyFloat(startT);
    InitMyFloat(stopT);
    InitMyFloat(nextT);
    for(i=0; i<7; i++) {InitMyFloat(xx[i]) };

    /* assign initials */
    xx[0] = r0[0];
    xx[1] = r0[1];
    xx[2] = r0[2];
    xx[3] = v0[0];
    xx[4] = v0[1];
    xx[5] = v0[2];
    xx[6] = m0;
    startT = 0;
    stopT = t;

    /* the main loop */
    T u_copy(u);
    if (t>0) {
        direction =1;
    } else {
        direction = -1;
        u_copy[0] = -u_copy[0];
        u_copy[1] = -u_copy[1];
        u_copy[2] = -u_copy[2];
    }

    do  {
        itmp = taylor_step_fixed_thrust( &startT, xx, direction, 1, log10tolerance, log10rtolerance, &stopT, &nextT, &order, mu, veff, u_copy[0], u_copy[1], u_copy[2]);
        } while(itmp == 0); /* while */
    r0[0] = xx[0];
    r0[1] = xx[1];
    r0[2] = xx[2];
    v0[0] = xx[3];
    v0[1] = xx[4];
    v0[2] = xx[5];
    m0 = xx[6];
}

}



#endif // KEP_TOOLBOX_PROPAGATE_TAYLOR_JORBA_H
