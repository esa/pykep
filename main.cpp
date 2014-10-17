/*****************************************************************************
 *   Copyright (C) 2004-2014 The PyKEP development team,                     *
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

#include <iostream>
#include <iomanip>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    //This is just an example on how to use the Keplerian toolbox in c++
    //Our recommendation is to use the Keplerian toolbox Python interface PyGMO
    //
    //You can compile this main linking it to the keplerian_toolbox library
    std::string line1("1 25544U 98067A   14266.21457330  .00014240  00000-0  25195-3 0  1363");
    std::string line2("2 25544  51.6476 333.7159 0002117 131.1360 342.8524 15.50569589906532");
	planet_tle sat1(line1,line2);
	array3D r,v;
	sat1.get_eph(sat1.get_ref_epoch() ,r,v);

	std::cout << sat1 << std::endl;
	std::cout << r << std::endl;
	std::cout << v << std::endl;
    return 0;
}
