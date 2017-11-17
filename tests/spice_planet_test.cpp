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

#include <iomanip>
#include <iostream>
#include <string>

#include "../src/planet/spice.h"

// In this test we test the functionality of the SPICE planet
using namespace kep_toolbox;

int main()
{
    // We start loading the kernel containing the 67P comet (by default its in planet::spice)
    util::load_spice_kernel("C_G_1000012_2012_2017.bsp");

    // We instantiate the object
    planet::spice pl1("CHURYUMOV-GERASIMENKO", "SUN", "ECLIPJ2000", "NONE");

    // We stream to screen the newly created planet
    std::cout << pl1 << std::endl;

    // We define the epoch to compute ephemeridess
    kep_toolbox::epoch when(kep_toolbox::epoch_from_string("2012-01-20 00:00:00.000"));
    array3D r, v;
    pl1.eph(when, r, v);

    std::cout << "67P eph at: " << when << std::endl;
    std::cout << r << v << std::endl;
    return failed_c();
}
