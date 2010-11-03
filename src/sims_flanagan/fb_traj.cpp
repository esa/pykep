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

#include "fb_traj.h"
#include "../astro_constants.h"
#include <algorithm>
#include <vector>
#include <numeric>
#include <boost/array.hpp>
#include "sc_state.h"
#ifdef USE_PAGMO
#include"../../exceptions.h"
#endif
#include"../core_functions/array3D_operations.h"

namespace kep_toolbox { namespace sims_flanagan{

/// Constructor.
/**
 * Constructs a fb_traj object from a sequence of planets (planet), the number of segments, and a spacecraft.
 * It also preallocate all necessary memory to store constraints violations and the "flight-plan"
 *
 * \param[in] sequence An STL vector containing the sequence of planets.
 * to the trajectory final goal. (e.g. Earth-Earth-Jupiter)
 * \param[in] n_seg An STL vector containing the number of segments to be used in each trajectory leg.
 * (e.g. 10-10-10)
 * \param[in] sc_ A spacecraft object containing details on the system engineering of the spacecraft.
 */

fb_traj::fb_traj(const std::vector<planet_ptr>& sequence, const std::vector<int>& n_seg, const spacecraft &sc_) : legs(sequence.size()-1),planets(sequence){
	//check consistency between sequence size and n_seg size
	if (sequence.size() -1 != n_seg.size()){
		throw_value_error("Inconsistency between planet sequence length and number of n_seg provided ");
	}
	//check consistency of n_seg
	for (unsigned int i=0;i<n_seg.size();i++){
		if (n_seg[i] <=1) throw_value_error("number of segments must be >=1");
	}

	//init spacecraft in all legs and allocate memory for throttles
	for (unsigned int i=0;i<n_seg.size();i++){
		legs[i].set_spacecraft(sc_);
		legs[i].set_throttles_size(n_seg[i]);
		legs[i].set_mu(sequence[i]->get_mu_central_body());
	}

	//init total_n_seg
	total_n_seg = std::accumulate(n_seg.begin(),n_seg.end(), 0);

}

/// Constructor.
/**
 * Constructs a fb_traj object from a sequence of planets (planet), the number of segments, and a spacecraft.
 * It also preallocate all necessary memory to store constraints violations and the "flight-plan"
 *
 * \param[in] sequence An STL vector containing the sequence of planet from the departure
 * to the trajectory final goal. (e.g. Earth-Earth-Jupiter)
 * \param[in] n_seg An STL vector containing the number of segments to be used in each trajectory leg.
 * (e.g. 10-10-10)
 * \param[in] mass_ The starting spacecraft wet_mass (kg)
 * \param[in] thrust_ The spacecraft maximum thrust (N)
 * \param[in] is_ The propulsion system specific impulse (sec)
 */

fb_traj::fb_traj(const std::vector<planet_ptr>& sequence, const std::vector<int>& n_seg, const double &mass_, const double &thrust_, const double &isp_) : legs(sequence.size()-1),planets(sequence){
	//check consistency between sequence size and n_seg size
	if (sequence.size() -1 != n_seg.size()){
		throw_value_error("Inconsistency between planet sequence length and number of n_seg provided ");
	}
	//check consistency of n_seg
	for (unsigned int i=0;i<n_seg.size();i++){
		if (n_seg[i] <=1) throw_value_error("number of segments must be >=1");
	}

	//init spacecraft in all legs and allocate memory for throttles
	for (unsigned int i=0;i<n_seg.size();i++){
		legs[i].set_spacecraft(spacecraft(mass_,thrust_,isp_));
		legs[i].set_throttles_size(n_seg[i]);
		legs[i].set_mu(sequence[0]->get_mu_central_body());
	}

	//init total_n_seg
	total_n_seg = std::accumulate(n_seg.begin(),n_seg.end(), 0);

}

/// Constructor.
/**
 * Constructs a fb_traj object from a sequence of planets (planet), the number of segments, and a spacecraft.
 * It also preallocate all necessary memory to store constraints violations and the "flight-plan"
 *
 * \param[in] sequence An STL vector containing the sequence of planet from the departure
 * to the trajectory final goal. (e.g. Earth-Earth-Jupiter)
 * \param[in] n_seg The number of segments equal for each leg (e.g 10)
 * \param[in] sc_ A spacecraft object containing details on the system engineering of the spacecraft.
 */

fb_traj::fb_traj(const std::vector<planet_ptr>& sequence, const unsigned int &n_seg, const spacecraft &sc_) : legs(sequence.size()-1),planets(sequence){

	//check consistency of n_seg
	if (n_seg < 1) throw_value_error("number of segments must be >=1");

	//init spacecraft in all legs and allocate memory for throttles
	for (size_t i=0;i<sequence.size()-1;i++){
		legs[i].set_spacecraft(sc_);
		legs[i].set_throttles_size(n_seg);
		legs[i].set_mu(sequence[i]->get_mu_central_body()) ;
	}

	//init total_n_seg
	total_n_seg = n_seg*(sequence.size()-1);

}

/// Constructor.
/**
 * Constructs a fb_traj object from a sequence of planets (planet), the number of segments, and a spacecraft.
 * It also preallocate all necessary memory to store constraints violations and the "flight-plan"
 *
 * \param[in] sequence An STL vector containing the sequence of planet from the departure
 * to the trajectory final goal. (e.g. Earth-Earth-Jupiter)
 * \param[in] n_seg The number of segments equal for each leg (e.g 10)
 * \param[in] mass_ The starting spacecraft wet_mass (kg)
 * \param[in] thrust_ The spacecraft maximum thrust (N)
 * \param[in] is_ The propulsion system specific impulse (sec)
 */

fb_traj::fb_traj(const std::vector<planet_ptr>& sequence, const unsigned int &n_seg, const double &mass_, const double &thrust_, const double &isp_) : legs(sequence.size()-1),planets(sequence){

	//check consistency of n_seg
	if (n_seg <1) throw_value_error("number of segments must be >=1");

	//init spacecraft in all legs and allocate memory for throttles
	for (size_t i=0;i<sequence.size()-1;i++){
		legs[i].set_spacecraft(spacecraft(mass_,thrust_,isp_));
		legs[i].set_throttles_size(n_seg);
		legs[i].set_mu(sequence[0]->get_mu_central_body());
	}

	//init total_n_seg
	total_n_seg = n_seg*(sequence.size()-1);

}


/// Copy constructor.
/**
 * Will deep-copy all the data members.
 *
 * @param[in] other source of the copy
 */
fb_traj::fb_traj(const fb_traj &other):legs(other.legs),total_n_seg(other.total_n_seg)
{
	for (std::vector<planet_ptr>::size_type i = 0; i < other.planets.size(); ++i) {
		planets.push_back(other.planets[i]->clone());
	}
}

/// Assignment operator.
/**
 * Will deep-copy all the data members.
 *
 * @param[in] other source of the assignment
 */
fb_traj &fb_traj::operator=(const fb_traj &other)
{
	// Protect against self-assignment.
	if (this == &other) {
		return *this;
	}
	legs = other.legs;
	planets.clear();
	for (std::vector<planet_ptr>::size_type i = 0; i < other.planets.size(); ++i) {
		planets.push_back(other.planets[i]->clone());
	}
	return *this;
}

/// Overload the stream operator for kep_toolbox::sims_flanagan::leg
/**
 * Streams out the leg object in a human readable format
 *
 * \param[in] s stream to which the output will be sent
 * \param[in] in fb_traj to be sent to stream
 *
 * \return reference to s
 *
 */

std::ostream &operator<<(std::ostream &s, const fb_traj &in ){
	(void)in;
	/*	    s << "Full Vector: " << std::endl;
	    for (size_t i=0;i<in.full_vector.size();i++) s << in.full_vector[i] << " ";
	    s << std::endl <<"Mismatches con: " << std::endl;
	    for (size_t i=0;i<in.mismatches_con.size();i++) s << in.mismatches_con[i] << " ";
	    s << std::endl <<"Throttles Magnitudes con (inequality): " << std::endl;
	    for (size_t i=0;i<in.throttles_con.size();i++) s << in.throttles_con[i] << " ";
	    s << std::endl <<"Same Relative Volicity con: " << std::endl;
	    for (size_t i=0;i<in.fb_relvel_con.size();i++) s << in.fb_relvel_con[i] << " ";
	    s << std::endl <<"Minimum altitude con (inequality): " << std::endl;
	    for (size_t i=0;i<in.fb_altitude_con.size();i++) s << in.fb_altitude_con[i] << " ";
	    s << std::endl <<"Same mass con: " << std::endl;
	    for (size_t i=0;i<in.fb_mass_con.size();i++) s << in.fb_mass_con[i] << " ";
	    s << std::endl <<"Same epoch con: " << std::endl;
	    for (size_t i=0;i<in.fb_epoch_con.size();i++) s << in.fb_epoch_con[i] << " ";*/
	return s;
}


}} //namespaces
