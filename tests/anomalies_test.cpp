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

#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <random>

#include "../src/core_functions/eq2par.h"
#include "../src/core_functions/par2eq.h"
#include "../src/core_functions/par2ic.h"
#include "../src/core_functions/eq2ic.h"
#include "../src/astro_constants.h"
#include "../src/io.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 2.0);
    std::uniform_real_distribution<> dis2(-1.0, 1.0);
    auto tol = 1e-10;
    bool fail = false;

    // We test that par2eq and eq2par are perfectly inverse of one another (wthin tol)
    for (auto i=0u; i<100000;++i){
        // We test on random arrays that the conversion is invertible
        std::vector<double> E = {dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis2(gen)};
        std::vector<double> EQ(6);
        std::vector<double> EINV(6);

        auto mu = dis(gen);
        // 1 - direct
        par2eq(EQ, E, false);
        eq2par(EINV, EQ, false);
        for (auto j = 0u; j<3;++j){
            if (std::abs(E[j] - EINV[j]) > tol) fail = true;
        }
        // for anomalies we need to check sin and cos
        for (auto j = 3u; j<6;++j){
            if ( std::abs(  std::sin(E[j]) - std::sin(EINV[j])) > tol) fail = true;
            if ( std::abs(std::cos(E[j]) - std::cos(EINV[j])) > tol) fail = true;
        }
        // 2 - retrograde
        par2eq(EQ, E, true);
        eq2par(EINV, EQ, true);
        for (auto j = 0u; j<3;++j){
            if (std::abs(E[j] - EINV[j]) > tol) fail = true;
        }
        // for anomalies we need to check sin and cos
        for (auto j = 3u; j<6;++j){
            if ( std::abs(  std::sin(E[j]) - std::sin(EINV[j])) > tol) fail = true;
            if ( std::abs(std::cos(E[j]) - std::cos(EINV[j])) > tol) fail = true;
        }
    }

    // We test that par2ic and eq2ic return the same values on the same orbit
    for (auto i=0u; i<100000;++i){
        auto mu = 1.;
        std::vector<double> E = {dis(gen), dis(gen), dis(gen), dis(gen), dis(gen), dis2(gen)};
        std::vector<double> EQ(6);
        std::vector<double> r(3),v(3),r1(3),v1(3);
        par2eq(EQ, E, true);
        par2ic(E, mu, r,v);
        eq2ic(EQ,mu, r1, v1, true);
        for (auto j = 0u; j<3;++j){
            if ( std::abs((r[j] - r1[j]) / std::max(std::abs(r[j]),1.)) > tol) {
                fail = true;
                print("Keplerian radius: ", r, "\n");
                print("Equinoctial radius: ", r1, "\n");
                print("Diff: ", std::abs((r[j] - r1[j]) / std::max(std::abs(r[j]),1.)), "\n");
                print("Keplerian Elements: ", E, "\n");
            }
            if ( std::abs((v[j] - v1[j]) / std::max(std::abs(v[j]),1.)) > tol) {
                fail = true;
                print("Keplerian velocity: ", v, "\n");
                print("Equinoctial velocity: ", v1, "\n");
                print("Diff: ", ((v[j] - v1[j]) / std::max(std::abs(v[j]),1.)), "\n");
                print("Keplerian Elements: ", E, "\n");
            }
        }
    }
   
    return fail;
}
