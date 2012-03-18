/*****************************************************************************
*   Copyright (C) 2004-2009 The PaGMO development team,                     *
*   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
*   http://apps.sourceforge.net/mediawiki/pagmo                             *
*   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
*   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
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

#include <vector>
#include <numeric>

#include "leg_s.h"
#include "sc_state.h"
#include "../astro_constants.h"
#include "../core_functions/array3D_operations.h"
#include "../core_functions/propagate_lagrangian.h"
#include"../exceptions.h"

namespace kep_toolbox{ namespace sims_flanagan{


std::string leg::human_readable() const {
	std::ostringstream s;
	s << *this;
	return s.str();
}


/// Overload the stream operator for kep_toolbox::sims_flanagan::leg
/**
 * Streams out the leg object in a human readable format
 *
 * \param[in] s stream to which the planet will be sent
 * \param[in] in leg to be sent to stream
 *
 * \return reference to s
 *
 */
std::ostream &operator<<(std::ostream &s, const leg &in ){
	s << std::setprecision(15);
	s << "Number of segments: " << in.m_throttles.size() << std::endl << std::endl;
	s << in.get_spacecraft() << std::endl;
	s << "Central body gravitational parameter: " << in.get_mu() << std::endl << std::endl;
	s << "Departure date: " << in.get_ti() << ", mjd2000: " << in.get_ti().mjd2000() << std::endl;
	s << "Arrival date: " << in.get_tf() << ", mjd2000: " << in.get_tf().mjd2000() << std::endl;
	s << "Initial mass: " << in.get_xi().get_mass() << " kg" << std::endl;
	s << "Final mass: " << in.get_xf().get_mass() << " kg" << std::endl;
	s << "State at departure: " << in.get_xi() << std::endl;
	s << "State at arrival: " << in.get_xf() << std::endl;

	s << std::endl << "Throttles values: " << std::endl;
	for (size_t i=0; i<in.get_throttles_size(); i++) {
		s << "\t\t\t" << in.m_throttles[i].get_value()[0] << " " << in.m_throttles[i].get_value()[1] << " " << in.m_throttles[i].get_value()[2] << std::endl;
	}

	try
	{
		s << std::endl << "Mismatch at the midpoint: ";
		s << in.compute_mismatch_con() << std::endl;
	}
	catch (...)
	{
		s << std::endl << "Mismatch at the midpoint: NUMERICAL ERROR!! COULD NOT CALCULATE THE STATE MISMATCH, CHECK YOUR DATA" << std::endl;
	}

	s << "Throttle magnitude constraints (if negative satisfied): " << in.compute_throttles_con();
	return s;
}

}} //namespaces

