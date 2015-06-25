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

#include <sstream>

#include "spacecraft.h"

namespace kep_toolbox{ namespace sims_flanagan{


std::string spacecraft::human_readable() const {
	std::ostringstream s;
	s << "NEP spacecraft:" << std::endl << std::endl;
	s << "mass: " << get_mass() << std::endl;
	s << "thrust: " << get_thrust() << std::endl;
	s << "isp: " << get_isp() << std::endl;
	return s.str();
};

std::ostream &operator<<(std::ostream &s, const spacecraft &in ) {
	s << "Spacecraft mass: " << in.get_mass() << std::endl;
	s << "Spacecraft thrust: " << in.get_thrust() << std::endl;
	s << "Spacecraft isp: " << in.get_isp();
	return s;
};

}} //end of namespaces
