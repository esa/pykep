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

#include <fstream>
#include <boost/algorithm/string.hpp>

#include "planet_mpcorb.h"
#include "../exceptions.h"

static const int mpcorb_format[12][2] =
{
	{92,11},	// a (AU)
	{70,9},		// e
	{59,9},		// i (deg)
	{48,9},		// Omega (deg)
	{37,9},		// omega (deg)
	{26,9},		// M (deg)
	{20,5},		// Epoch (shitty format)
	{166,28},	// Asteroid readable name
	{8,5},		// Absolute Magnitude
	{117,5},	// Number of observations
	{123,3},	// Number of oppositions
	{127,4},	// Year of First Observation (only if the number of oppositions is larger than 1)
};

namespace kep_toolbox{

planet_mpcorb::planet_mpcorb(const std::string& line)
{
	std::string linecopy(line);
	boost::algorithm::to_lower(linecopy);
	array6D elem;
	std::string tmp;
	//read keplerian elements from MPCORB.DAT
	for (int i = 0; i < 6; ++i) {
		tmp.clear();
		tmp.append(&linecopy[mpcorb_format[i][0]],mpcorb_format[i][1]);
		boost::algorithm::trim(tmp);
		elem[i] = boost::lexical_cast<double>(tmp);
	}
	// Converting orbital elements to the dictatorial PaGMO units.
	elem[0] *= ASTRO_AU;
	for (int i = 2; i < 6; ++i) {
		elem[i] *= ASTRO_DEG2RAD;
	}

	// Deal with MPCORB data format
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[6][0]],mpcorb_format[6][1]);
	boost::algorithm::trim(tmp);
	kep_toolbox::epoch epoch(packed_date2epoch(tmp));

	// Extract absolute magnitude
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[8][0]],mpcorb_format[8][1]);
	boost::algorithm::trim(tmp);
	m_H = boost::lexical_cast<double>(tmp);

	// Extract number of observations
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[9][0]],mpcorb_format[9][1]);
	boost::algorithm::trim(tmp);

	m_n_observations = boost::lexical_cast<unsigned int>(tmp);

	// Extract number of oppositions
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[10][0]],mpcorb_format[10][1]);
	boost::algorithm::trim(tmp);

	m_n_oppositions = boost::lexical_cast<unsigned int>(tmp);

	// Extract the year of discovery (if only one observation is made this field will, instead, contain the Arc length in Days)
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[11][0]],mpcorb_format[11][1]);
	boost::algorithm::trim(tmp);

	m_year_of_discovery = boost::lexical_cast<unsigned int>(tmp);

	//Now we estimate the asteroid radius, safe_radius and gravity parametes with hyper simplified assumptions
	double radius = 1329000 * std::pow(10,-m_H * 0.2); // This is assuming an albedo of 0.25 (www.physics.sfasu.edu/astro/asteroids/sizemagnitude.html)
	double mu_planet = 4./3. * M_PI * std::pow(radius,3) * 2800 * ASTRO_CAVENDISH;

	// Record asteroid name.
	tmp.clear();
	tmp.append(&linecopy[mpcorb_format[7][0]],mpcorb_format[7][1]);
	boost::algorithm::trim(tmp);

	build_planet(epoch,elem,ASTRO_MU_SUN,mu_planet,radius,radius*1.1,tmp);
}


epoch planet_mpcorb::packed_date2epoch(std::string in) {
	if (in.size()!=5) {
		throw_value_error("mpcorb data format requires 5 characters.");
	}
	boost::algorithm::to_lower(in);
	boost::gregorian::greg_year anno = packed_date2number(in[0]) * 100 + boost::lexical_cast<int>(std::string(&in[1],&in[3]));
	boost::gregorian::greg_month mese = packed_date2number(in[3]);
	boost::gregorian::greg_day giorno = packed_date2number(in[4]);
	return epoch(anno,mese,giorno);
}
// Convert mpcorb packed dates convention into number. (lower case assumed)
// TODO: check locale ASCII.
int planet_mpcorb::packed_date2number(char c)
{
	return static_cast<int>(c) - (boost::algorithm::is_alpha()(c) ? 87 : 48);
}

planet_ptr planet_mpcorb::clone() const
{
	return planet_ptr(new planet_mpcorb(*this));
}

} //namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet_mpcorb);
// Serialization code (END)
