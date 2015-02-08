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

#ifndef ASTEROID_GTOC7_H
#define ASTEROID_GTOC7_H

// Serialization code
#include "serialization.h"
// Serialization code (END)

#include "planet.h"
#include "config.h"

namespace kep_toolbox{

/// A GTOC7 asteroid
/**
 * This class derives from the planet class and allow to instantiate asteroids
 * from the Global Trajectory Optimization Competition (GTOC) 7th edition
 *
 * @see http://sophia.estec.esa.int/gtoc_portal/
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE asteroid_gtoc7 : public planet
{
public:
	/// Constructor
	/**
	 * Construct from a consecutive id from 0 (Earth) to 16256 . The order is
	 * that of the original data file from Turin
	 * Earth: 0
	 * Asteroid: 1 - 16256

	 * \param[in] ast_id asteroid id
	 */
	asteroid_gtoc7(int = 0);

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
BOOST_CLASS_EXPORT_KEY(kep_toolbox::asteroid_gtoc7);
// Serialization code (END)

#endif //ASTEROID_GTOC7_H
