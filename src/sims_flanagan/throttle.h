#ifndef THROTTLE_H
#define THROTTLE_H

#include <numeric>
#include<iostream>

// Serialization code
#include "../serialization.h"
// Serialization code (END)

#include "../astro_constants.h"
#include "../epoch.h"
#include "../config.h"

namespace kep_toolbox { namespace sims_flanagan{

/// A single throttle
/**
 * This class models a single throttle in the Sims-Flanagan model. It essentialy contains the cartesian
 * components of one throttle (non dimensional impulse)
 *impulse
 * @author David di Lorenzo
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE throttle {
public:
	throttle():m_start(),m_end() {
		m_value[0] = 0;
		m_value[1] = 0;
		m_value[2] = 0;
	}

	throttle(const epoch &_start, const epoch &_end, const array3D& _value)
		: m_start(_start), m_end(_end), m_value(_value) {}

	const epoch& get_start() const {
		return m_start;
	}

	const epoch& get_end() const {
		return m_end;
	}

	const array3D& get_value() const {
		return m_value;
	}

	void set_start(const epoch &e) {m_start=e;}
	void set_end(const epoch &e) {m_end=e;}
	void set_value(const array3D& e) {m_value=e;}

	double get_norm() const {
		return std::sqrt(std::inner_product(m_value.begin(), m_value.end(), m_value.begin(), 0.));
	}
	std::string human_readable() const {
		std::ostringstream s;
		s << "start = " << m_start << std::endl;
		s << "value = " << m_value << std::endl;
		s << "end = " << m_end;
		return s.str();
	}
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & m_start;
		ar & m_end;
		ar & m_value;
	}
// Serialization code (END)
	epoch m_start;
	epoch m_end;
	array3D m_value;
};

}} //Namespaces

#endif
