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

#ifndef PLANET_MPCORB_H
#define PLANET_MPCORB_H

#include <string>

// Serialization code
#include "serialization.h"
// Serialization code (END)

#include "planet.h"
#include "config.h"


namespace kep_toolbox{

/// Minor Planet (keplerian)
/**
 * This class derives from the planet class and allow to instantiate planets of
 * from the MPCORB database using their names or row id. The file MPCORB.DAT is searched
 * in the current directory.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class __KEP_TOOL_VISIBLE planet_mpcorb : public planet
{
public:
	/**
	 * Construct a minor planet from a line of the MPCORB.DAT file. Default value is the MPCORB.DAT line
	 * for the dwarf planet Ceres.
	 * \param[in] name a string containing one line of MPCORB.DAT
	 */
	planet_mpcorb(const std::string & = "00001    3.34  0.12 K107N 113.41048   72.58976   80.39321   10.58682  0.0791382  0.21432817   2.7653485  0 MPO110568  6063  94 1802-2006 0.61 M-v 30h MPCW       0000      (1) Ceres              20061025");
	planet_ptr clone() const;
	static epoch packed_date2epoch(std::string);
	double get_H() const {return m_H;};
	unsigned int get_n_observations() const {return m_n_observations;};
	unsigned int get_n_oppositions() const {return m_n_oppositions;};
	unsigned int get_year_of_discovery() const {return m_year_of_discovery;};

private:
// Serialization code
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
		ar & m_H;
		ar & m_n_observations;
		ar & m_n_oppositions;
		ar & m_year_of_discovery;
	}
// Serialization code (END)

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


} /// End of namespace kep_toolbox

// Serialization code
BOOST_CLASS_EXPORT_KEY(kep_toolbox::planet_mpcorb);
// Serialization code (END)

#endif // PLANET_MPCORB_H
