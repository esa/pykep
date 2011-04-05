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

#ifndef ASTEROID_GTOC2_H
#define ASTEROID_GTOC2_H

// Serialization code
#include "serialization.h"
// Serialization code (END)

#include "planet.h"
#include "config.h"

namespace kep_toolbox{

/// A GTOC2 asteroid
/**
 * This class derives from the planet class and allow to instantiate asteroids
 * from the Global Trajectory Optimization Competition (GTOC) 2nd edition
 *
 * @see http://www.esa.int/gsp/ACT/mad/op/GTOC/index.htm
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */

class __KEP_TOOL_VISIBLE asteroid_gtoc2 : public planet
{
public:
	/// Constructor
	/**
	 * Construct from a consecutive id from 0 to 910 (Earth). The order is that of the original
	 * data file from JPL
	 * Group 1:   0 - 95
	 * Group 2:  96 - 271
	 * Group 3: 272 - 571
	 * Group 4: 572 - 909
	 * Earth:   910
	 * \param[in] name a string describing a planet
	 */
	asteroid_gtoc2(const int & = 0);

	/// Getter
	/**
	 * Gets the group id of the asteroid as defined in the original JPL data file
	 *
	 */
	int get_group() const;
	planet_ptr clone() const;
private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
		ar & m_group;
	}
// Serialization code (END)
	int m_group;
};
} // Namespaces

// Serialization code
BOOST_CLASS_EXPORT(kep_toolbox::asteroid_gtoc2);
// Serialization code (END)

#endif //ASTEROID_GTOC2_H
