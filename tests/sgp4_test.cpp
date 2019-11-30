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

#include <keplerian_toolbox/core_functions/array3D_operations.hpp>
#include <keplerian_toolbox/planet/tle.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace kep_toolbox;

// In this test we read the historical two line elements of a GPS satellite
// from the file sgp4_test.txt
// These represent actual observations made on the satellite r and v.
// We then compute r and v propagating the TLE of the previous pobservation
// up to the current observation and measure the error. The test passes if
// the max error is less than 10 km on the position and 15 m/s/ on the velocity

int main()
{
    std::ifstream txtfile("sgp4_test.txt");
    std::string line1, line2, line21, line22;
    array3D r1, v1, r2, v2, err_r, err_v;
    std::vector<double> errors_r, errors_v;
    if (!txtfile.good()) {
        std::cout << "File sgp4_test.txt not found" << std::endl;
        return 1; // exit if file not found
    }
    getline(txtfile, line1);
    getline(txtfile, line2);
    while (!txtfile.eof()) {
        getline(txtfile, line21);
        getline(txtfile, line22);
        planet::tle sat1(line1.substr(0, 69), line2.substr(0, 69));
        planet::tle sat2(line21.substr(0, 69), line22.substr(0, 69));
        sat1.eph(sat2.get_ref_mjd2000(), r1, v1);
        sat2.eph(sat2.get_ref_mjd2000(), r2, v2);
        diff(err_r, r1, r2);
        diff(err_v, v1, v2);
        errors_r.push_back(norm(err_r));
        errors_v.push_back(norm(err_v));
        if (norm(err_r) > 10000) {
            std::cout << line1 << std::endl;
            std::cout << line2 << std::endl;
        }
        line1 = line21;
        line2 = line22;
    }
    double sum = std::accumulate(errors_r.begin(), errors_r.end(), 0.0);
    double mean_r = sum / static_cast<double>(errors_r.size());

    sum = std::accumulate(errors_v.begin(), errors_v.end(), 0.0);
    double mean_v = sum / static_cast<double>(errors_v.size());
    double max_r = *std::max_element(errors_r.begin(), errors_r.end());
    double max_v = *std::max_element(errors_v.begin(), errors_v.end());

    std::cout << "Mean error r (m): " << mean_r << std::endl;
    std::cout << "Mean error v (m/s): " << mean_v << std::endl;
    std::cout << "Max error r (m/s): " << max_r << std::endl;
    std::cout << "Max error v (m/s): " << max_v << std::endl;

    if (max_r < 10000 && max_v < 15) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}
