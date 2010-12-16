#include"spacecraft.h"

namespace kep_toolbox{ namespace sims_flanagan{


std::string spacecraft::human_readable() const {
	std::ostringstream s;
	s << "NEP spacecraft:" << std::endl << std::endl;
	s << "mass: " << get_mass() << std::endl;
	s << "thrust: " << get_thrust() << std::endl;
	s << "isp: " << get_isp() << std::endl;
	return s.str();
}
std::ostream &operator<<(std::ostream &s, const spacecraft &in ) {
	s << "Spacecraft mass: " << in.get_mass() << std::endl;
	s << "Spacecraft thrust: " << in.get_thrust() << std::endl;
	s << "Spacecraft isp: " << in.get_isp() << std::endl;
	return s;
};

}} //end of namespaces
