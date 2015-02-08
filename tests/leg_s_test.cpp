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

#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

#include "../src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;
int main() {
    int n_seg=15;
	double mu = ASTRO_MU_SUN;
	sims_flanagan::spacecraft sc = sims_flanagan::spacecraft(1000,0.1,2000);
    sims_flanagan::leg_s phase1(n_seg,pow(ASTRO_AU,-1.5), 1.5);
	phase1.set_mu(mu);
	phase1.set_sc(sc);
	planet_ss earth("earth");
	array3D r,v;
	earth.get_eph(epoch(0),r,v);
	sims_flanagan::sc_state x0(r,v,sc.get_mass());
	earth.get_eph(epoch(100),r,v);
	sims_flanagan::sc_state xf(r,v,sc.get_mass()/2);
	std::vector<double> throttles(n_seg*3,0.1423);
	phase1.set_leg(epoch(0),x0,throttles,epoch(100),xf,1.5*365.25*ASTRO_DAY2SEC,sc,mu);
    for (int i=0; i< 8;++i){
        std::cout << phase1.compute_mismatch_con()[i] << ", ";
    }
    std::cout<< std::endl;
    for (int i=0; i< n_seg;++i){
        std::cout << phase1.compute_throttles_con()[i] << ", ";
    }
    std::cout<< std::endl;
	return 0;
}

