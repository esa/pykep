/*****************************************************************************
 *   Copyright (C) 2004-2015 The pykep development team,                     *
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

#ifndef KEP_TOOLBOX_PLANET_SPICE_H
#define KEP_TOOLBOX_PLANET_SPICE_H

#include <string>

#include "base.h"
#include "../util/spice_utils.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox{ namespace planet {

/// A planet using the SPICE Toolbox
/**
 * This class allows to instantiate a planet that uses SPICE toolbox
 * to compute the ephemerides.
 *
 * NOTE: The class does not check upon construction that the required kernels are loaded 
 * in memory. Its only when the ephemerides are actually called that an exception is thrown
 * in case the required kernels are not loaded
 *
 * @see http://naif.jpl.nasa.gov/naif/toolkit.html
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE spice : public base
{
public:
	spice(const std::string & = "CHURYUMOV-GERASIMENKO", 
		const std::string & = "SUN", 
		const std::string & = "ECLIPJ2000", 
		const std::string & = "NONE",
		double = 0, // mu_central_body
		double = 0, // mu_self
		double = 0, // radius
		double = 0  // safe_radius
	);
	planet_ptr clone() const;
	std::string human_readable_extra() const;

private:
	void eph_impl(double mjd2000, array3D &r, array3D &v) const;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<std::string& >(m_target);
		ar & const_cast<std::string& >(m_observer);
		ar & const_cast<std::string& >(m_reference_frame);
		ar & const_cast<std::string& >(m_aberrations);
	}

	const std::string m_target;
	const std::string m_observer;
	const std::string m_reference_frame;
	const std::string m_aberrations;

	// Dummy variables that store intermidiate values to transfer to and from SPICE 
	mutable SpiceDouble m_state[6];
	mutable SpiceDouble m_lt;

};


}} /// End of namespace kep_toolbox

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::spice)

#endif // KEP_TOOLBOX_PLANET_SPICE_H
