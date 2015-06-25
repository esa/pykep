/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
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

#ifndef KEP_TOOLBOX_PLANET_MPCORB_H
#define KEP_TOOLBOX_PLANET_MPCORB_H

#include <string>

#include "keplerian.h"
#include "../serialization.h"
#include "../config.h"


namespace kep_toolbox { namespace planet {

/// Minor Planet (keplerian)
/**
 * This class allows to instantiate keplerian planets from the MPCORB database.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE mpcorb : public keplerian
{
public:
	mpcorb(const std::string & = "00001    3.34  0.12 K107N 113.41048   72.58976   80.39321   10.58682  0.0791382  0.21432817   2.7653485  0 MPO110568  6063  94 1802-2006 0.61 M-v 30h MPCW       0000      (1) Ceres              20061025");
	planet_ptr clone() const;

	static epoch packed_date2epoch(std::string);
	double get_H() const {return m_H;};
	unsigned int get_n_observations() const {return m_n_observations;};
	unsigned int get_n_oppositions() const {return m_n_oppositions;};
	unsigned int get_year_of_discovery() const {return m_year_of_discovery;};

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<keplerian>(*this);
		ar & m_H;
		ar & m_n_observations;
		ar & m_n_oppositions;
		ar & m_year_of_discovery;
	}

	static int packed_date2number(char c);
	// Absolute Magnitude
	double m_H;
	// Number of observations
	unsigned int m_n_observations;
	// Number of oppositions
	unsigned int m_n_oppositions;
	// Year the asteroid was first discovered
	unsigned int m_year_of_discovery;
};


}} /// End of namespace kep_toolbox

BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet::mpcorb);

#endif // KEP_TOOLBOX_PLANET_MPCORB_H
