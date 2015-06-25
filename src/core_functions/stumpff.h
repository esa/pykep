/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
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


#ifndef KEP_TOOLBOX_STUMPFF_H
#define KEP_TOOLBOX_STUMPFF_H

#include<cmath>

namespace kep_toolbox {

inline double stumpff_s(const double x) {
    if (x > 0)
    {
        return (sqrt(x) - sin(std::sqrt(x)))/pow(sqrt(x),3);
    }
    else if (x < 0)
    {
        return (std::sinh(std::sqrt(-x)) - sqrt(-x))/pow(-x,3./2);
    }
    else
    {
        return (1./6.);
    }
}


inline double stumpff_c(const double x) {
    if (x > 0)
    {
        return (1 - cos(sqrt(x)))/x;
    }
    else if (x < 0)
    {
        return (std::cosh(sqrt(-x)) - 1)/(-x);
    }
    else
    {
        return 0.5;
    }

}

}

#endif // KEP_TOOLBOX_STUMPFF_H
