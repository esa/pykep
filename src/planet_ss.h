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

#ifndef PLANET_SS_H
#define PLANET_SS_H


#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
#include "serialization.h"
#endif

#include"planet.h"

namespace kep_toolbox{

/// Solar System Planet (keplerian)
/**
 * This class derives from the planet class and allow to instantiate planets of
 * the solar system by referring to their common names.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class planet_ss : public planet
{
public:
	/**
	 * Construct a planet from its common name (e.g. VENUS)
	 * \param[in] name a string describing a planet
	 */
	planet_ss(const std::string & = "earth");
	planet_ptr clone() const;
private:
#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
	}
#endif
};


} /// End of namespace kep_toolbox

#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT(kep_toolbox::planet_ss);
#endif

#endif // PLANET_SS_H
