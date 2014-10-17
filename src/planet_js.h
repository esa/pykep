/*****************************************************************************
 *   Copyright (C) 2004-2014 The PyKEP development team,                     *
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

#ifndef PLANET_JS_H
#define PLANET_JS_H

#include "serialization.h"
#include "planet.h"
#include "config.h"

namespace kep_toolbox{

/// Solar System Planet (keplerian)
/**
 * This class derives from the planet class and allow to instantiate moons of
 * the Jupiter system by referring to their common names. Ephemerides are those
 * used during the GTOC6 competition
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE planet_js : public planet
{
public:
	/**
	 * Construct a Jupiter moon from its common name
	 * \param[in] name a string describing a planet
	 */
	planet_js(const std::string & = "io");
	planet_ptr clone() const;
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
	}
// Serialization code (END)
};


} /// End of namespace kep_toolbox

// Serialization code
BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet_js);
// Serialization code (END)

#endif // PLANET_JS_H
