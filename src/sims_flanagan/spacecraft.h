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

#ifndef SPACECRAFT_H
#define SPACECRAFT_H

#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
#include "../serialization.h"
#endif

namespace kep_toolbox {

/// Spacecraft
/**
 * A container for system design parameters of a spacecraft.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class spacecraft
{
public:
	spacecraft():mass(0),thrust(0),isp(0) {}
	spacecraft(const double &mass_, const double &thrust_, const double &isp_) : mass(mass_),thrust(thrust_),isp(isp_) {}
	const double& get_mass() const {return mass;}
	const double& get_thrust() const {return thrust;}
	const double& get_isp() const {return isp;}
private:
#ifdef KEP_TOOLBOX_ENABLE_SERIALIZATION
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & mass;
		ar & thrust;
		ar & isp;
	}
#endif
	double mass;
	double thrust;
	double isp;
};

} //Namespaces

#endif // SPACECRAFT_H
