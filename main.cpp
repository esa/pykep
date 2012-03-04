/*****************************************************************************
 *   Copyright (C) 2004-2012 The PyKEP development team,                     *
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
    //This is just an example on how to use the keplerian toolbox in c++
    //Our reccomendation is to use the keplerian toolbox Python interface PyGMO
    //
    //You can comile this main linking it to the keplerian_toolbox library
    

    array3D r0 = {{0.22, 1, 2}};
    array3D v0 = {{0.012, -0.2, 0.1}};
    array3D rcp(r0); array3D vcp(v0);
    //array3D v0 = {{0, 1, 0}};
    const array3D u = {{0.001,0.0023,-0.0023}};

    double m0=1000.0; double t0=0;
    std::cout << "r0: " << r0 << std::endl << "v0: "<< v0 << std::endl << "m0: " << m0 << std::endl << "t0: " << t0 << std::endl;
    double sf = 1.2;
    propagate_taylor_s(r0,v0,m0,t0,u,sf, 1.1, 1.1, 1.1, 1.5,-14,-14);
    std::cout << "r0: " << r0 << std::endl << "v0: "<< v0 << std::endl << "m0: " << m0 << std::endl << "t0: " << t0 << std::endl;
    propagate_taylor(r0,v0,m0,u,-t0, 1.1, 1.1,-14,-14);
    std::cout << "r0: " << r0 << std::endl << "v0: "<< v0 << std::endl << "m0: " << m0 << std::endl << "t0: " << t0 << std::endl;
    return 0;
}
