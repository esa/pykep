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

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <iomanip>
#include <iostream>

#include <keplerian_toolbox/core_functions/array3D_operations.hpp>
#include <keplerian_toolbox/core_functions/propagate_taylor_s.hpp>

using namespace std;
using namespace kep_toolbox;
int main()
{
    // Preamble
    array3D r0, v0, r0_cp, v0_cp, u;
    double t0, m0, s, t0_cp, m0_cp;
    boost::mt19937 rng;
    boost::uniform_int<> dist(0, 1);
    boost::variate_generator<boost::mt19937 &, boost::uniform_int<>> rand_bit(rng, dist);
    boost::uniform_real<> dist1(-1, 1);
    boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> drng(rng, dist1);
    double acc = 0, err_max = 0, err = 0, max_t0 = 0;
    int count = 0;

    // Experiment Settings
    unsigned int Ntrials = 30000;

    // Start Experiment
    for (unsigned int i = 0; i < Ntrials; ++i) {
        // 1 - generate a random propagation set-up
        r0[0] = drng() * 2;
        r0[1] = drng() * 2;
        r0[2] = drng() * 2;
        v0[0] = drng() * 1;
        v0[1] = drng() * 1;
        v0[2] = drng() * 1;
        m0 = (drng() + 1) * 500 + 1000;
        t0 = 0;
        u[0] = drng() * 1;
        u[1] = drng() * 1;
        u[2] = drng() * 1;
        s = drng() * 1;
        r0_cp = r0;
        v0_cp = v0;
        m0_cp = m0;
        t0_cp = t0;
        // 2 - propagate back and forth
        try {
            propagate_taylor_s(r0, v0, m0, t0, u, s, 1.1, 1.0, 1.0, 1.0, -14, -14);
            max_t0 = std::max(t0, max_t0);
            propagate_taylor_s(r0, v0, m0, t0, u, -s, 1.1, 1.0, 1.0, 1.0, -14, -14);
        } catch (...) {
            std::cout << "Exception raised: " << std::endl;
            std::cout << "r0: " << r0_cp << std::endl;
            std::cout << "v0: " << v0_cp << std::endl;
            std::cout << "m0: " << m0_cp << std::endl;
            std::cout << "t0: " << t0_cp << "tf: " << t0 << std::endl;
            ;
            std::cout << "u: " << u << std::endl;
            std::cout << "s: " << s << std::endl;
        }
        diff(r0_cp, r0, r0_cp);
        err = norm(r0_cp);
        err_max = std::max(err_max, err);
        acc += err;
        count++;
    }
    std::cout << "Maximum time integrated: " << max_t0 << std::endl;
    std::cout << "Max error: " << err_max << std::endl;
    std::cout << "Average Error: " << acc / count << std::endl;
    std::cout << "Number of Propagations Made: " << count << std::endl;
    if (err_max < 1e-7) {
        return 0;
    } else {
        return 1;
    }
}
