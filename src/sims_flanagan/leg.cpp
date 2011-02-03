
#include <vector>
#include <numeric>

#include "leg.h"
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
};


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
	s << "High-fidelity propagation: " << in.get_high_fidelity() << std::endl;
	s << "Number of segments: " << in.throttles.size() << std::endl << std::endl;
	s << in.get_spacecraft() << std::endl;
	s << "Central body gravitational parameter: " << in.get_mu() << std::endl << std::endl;
	s << "Departure date: " << in.get_t_i() << ", mjd2000: " << in.get_t_i().mjd2000() << std::endl;
	s << "Arrival date: " << in.get_t_f() << ", mjd2000: " << in.get_t_f().mjd2000() << std::endl;
	s << "Initial mass: " << in.get_x_i().get_mass() << " kg" << std::endl;
	s << "Final mass: " << in.get_x_f().get_mass() << " kg" << std::endl;
	s << "State at departure: " << in.get_x_i() << std::endl;
	s << "State at arrival: " << in.get_x_f() << std::endl;

	s << std::endl << "Throttles values: " << std::endl;
	for (size_t i=0; i<in.get_throttles_size(); i++) {
		s << "\t\t\t" << in.throttles[i].get_value()[0] << " " << in.throttles[i].get_value()[1] << " " << in.throttles[i].get_value()[2] << std::endl;
	}

	std::vector<double> temp(in.get_throttles_size());
	in.get_throttles_con(temp.begin(), temp.end());
	sc_state mism;
	try
	{
		in.get_mismatch_con(mism);
		s << std::endl << "Mismatch at the midpoint: ";
		s << mism.get_position() << " " << mism.get_velocity() << " " << mism.get_mass() << std::endl;
	}
	catch (...)
	{
		s << std::endl << "Mismatch at the midpoint: NUMERICAL ERROR!! COULD NOT CALCULATE THE STATE MISMATCH, CHECK YOUR DATA" << std::endl;
	}

	s << "Throttle magnitude constraints (if negative satisfied): [";
	for (size_t i=0;i< in.get_throttles_size();i++) s << temp[i] << " ";
	s << "]";
	return s;
}

}} //namespaces

