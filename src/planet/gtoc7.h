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

#ifndef KEP_TOOLBOX_PLANET_GTOC7_H
#define KEP_TOOLBOX_PLANET_GTOC7_H

#include "keplerian.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox{ namespace planet {

/// A GTOC7 asteroid
/**
 * This class allows to instantiate main belt asteroids
 * from the Global Trajectory Optimization Competition (GTOC) 7th edition
 *
 * @see http://sophia.estec.esa.int/gtoc_portal/
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE gtoc7 : public keplerian
{
public:
	gtoc7(int = 0);
	planet_ptr clone() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<keplerian>(*this);
	}
};


}} /// namespaces

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::gtoc7);


#endif // KEP_TOOLBOX_PLANET_GTOC7_H
