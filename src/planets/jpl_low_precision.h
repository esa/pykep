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

#ifndef KEP_TOOLBOX_PLANET_JPL_LP_H
#define KEP_TOOLBOX_PLANET_JPL_LP_H

#include "base.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox{ namespace planets {

/// Solar System Planet (jpl simplified ephemerides)
/**
 * This class allows to instantiate planets of
 * the solar system by referring to their common names. The ephemeris used
 * are low_precision ephemeris taken from http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
 * valid in the timeframe 1800AD - 2050 AD
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE jpl_lp : public base
{
public:

	jpl_lp(const std::string & = "earth");
	planet_ptr clone() const;
	std::string human_readable_extra() const;

private:
	void eph_impl(double mjd2000, array3D &r, array3D &v) const;

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & jpl_elements;
		ar & jpl_elements_dot;
		ar & const_cast<double&>(ref_mjd2000);
	}

	array6D jpl_elements;
	array6D jpl_elements_dot;
	const double ref_mjd2000;
};


}} /// End of namespaces

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planets::jpl_lp)

#endif // KEP_TOOLBOX_PLANET_JPL_LP_H
