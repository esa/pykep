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

#ifndef KEP_TOOLBOX_PLANET_SS_H
#define KEP_TOOLBOX_PLANET_SS_H

#include "planet.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox{

/// Solar System Planet (keplerian)
/**
 * This class derives from the planet class and allow to instantiate planets of
 * the solar system by referring to their common names. The ephemeris used
 * are low_precision ephemeris taken from http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
 * valid in the timeframe 1800AD - 2050 AD
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE planet_ss : public planet
{
public:
	/**
	 * Construct a planet from its common name (e.g. VENUS)
	 * \param[in] name a string describing a planet
	 */
	planet_ss(const std::string & = "earth");
	planet_ptr clone() const;
	/// Computes the planet/system position and velocity w.r.t the Sun
	/**
		* \param[in] when Epoch in which ephemerides are required
		* \param[out] r Planet position at epoch (SI units)
		* \param[out] v Planet velocity at epoch (SI units)
		*/
	void get_eph(const epoch& when, array3D &r, array3D &v) const;
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
		ar & jpl_elements;
		ar & jpl_elements_dot;
	}
// Serialization code (END)

	array6D jpl_elements;
	array6D jpl_elements_dot;
};


} /// End of namespace kep_toolbox

// Serialization code
BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet_ss)
// Serialization code (END)

#endif // KEP_TOOLBOX_PLANET_SS_H
