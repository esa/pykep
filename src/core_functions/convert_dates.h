/*****************************************************************************
 *   Copyright (C) 2004-2015 The pykep development team,                     *
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

#ifndef KEP_TOOLBOX_CONVERT_DATES_H
#define KEP_TOOLBOX_CONVERT_DATES_H

namespace kep_toolbox {
    inline double jd2mjd(const double & in){
        return ( in - 2400000.5 );
    }
    inline double jd2mjd2000(const double & in){
        return ( in - 2451544.5 );
    }
    inline double mjd2jd(const double & in){
        return ( in + 2400000.5 );
    }
    inline double mjd2mjd2000(const double & in){
        return ( in - 51544 );
    }
    inline double mjd20002jd(const double & in){
        return ( in + 2451544.5 );
    }
    inline double mjd20002mjd(const double & in){
        return ( in + 51544 );
    }
}

#endif // KEP_TOOLBOX_CONVERT_DATES_H
