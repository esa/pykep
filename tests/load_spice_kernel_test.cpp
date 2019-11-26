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

#include <keplerian_toolbox/epoch.hpp>
#include <keplerian_toolbox/util/spice_utils.hpp>

// In this test we test the functionality of loading spice kernels
// Summary for: C_G_1000012_2012_2017.bsp

// Body: CHURYUMOV-GERASIMENKO (1000012) w.r.t. SUN (10)
//      Start of Interval (ET)              End of Interval (ET)
//      -----------------------------       -----------------------------
//      2012 JAN 01 00:00:00.000            2017 JAN 01 00:00:00.000

using namespace kep_toolbox::util;

std::string stream(double input[6])
{
    std::ostringstream retval;
    retval << "[";
    for (int i = 0; i < 6; ++i) {
        retval << std::setprecision(16) << input[i] << " ";
    }
    retval << "]";
    return retval.str();
}

int main()
{
    // We start loading the four kernels shipped with the kep_toolbox
    load_spice_kernel("pck00010.tpc");
    load_spice_kernel("gm_de431.tpc");
    int dim = 0;

    // Some definitions
    SpiceDouble radii[3];
    SpiceDouble mu_mars[1];

    // We define the epoch to compute ephemeridess
    kep_toolbox::epoch when(kep_toolbox::epoch_from_string("2012-01-20 00:00:00.000"));

    char stringa[] = "RETURN";

    // We set SPICE to allow error handling
    erract_c("SET", 0, stringa);

    // We check if the kernels have been loaded correctly by extracting a few
    // properties of mars
    bodvrd_c("MARS", "RADII", 3, &dim, radii);
    std::cout << "Mars Radius in km: " << std::setprecision(16) << radii[0] << std::endl;

    bodvrd_c("MARS", "GM", 1, &dim, mu_mars);
    std::cout << "Mars gravity parameter in km: " << std::setprecision(16) << mu_mars[0] << std::endl << std::endl;

    return failed_c();
}
