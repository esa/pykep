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

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/none.hpp>

#include "keplerian.h"
#include "../core_functions/ic2par.h"
#include "../core_functions/par2ic.h"
#include "../core_functions/convert_anomalies.h"
#include "../exceptions.h"

namespace kep_toolbox{ namespace planet {

const array6D keplerian::default_elements = {{1.0,0.1,0.1,0.1,0.1,0.1}};

/// Constructor
/**
* Constructs a planet from its elements and its phyisical parameters
* \param[in] ref_epoch epoch to which the elements are referred to
* \param[in] elem A STL vector containing the keplerian parameters (a,e,i,Om,om,M). (SI units)
* \param[in] mu_central_body The gravitational parameter of the attracting body (SI units)
* \param[in] mu_self The gravitational parameter of the planet (SI units)
* \param[in] radius radius of the planet (SI units)
* \param[in] safe_radius mimimual radius that is safe during a fly-by of the planet (SI units)
* \param[in] name C++ string containing the planet name. Default value is "Unknown"
*/
keplerian::keplerian(
	const epoch& ref_epoch, 
	const array6D& keplerian_elements, 
	double mu_central_body, 
	double mu_self, 
	double radius, 
	double safe_radius, 
	const std::string &name) 
	: base(mu_central_body, mu_self, radius, safe_radius, name), m_keplerian_elements(keplerian_elements), m_ref_mjd2000(ref_epoch.mjd2000())
{
	if (keplerian_elements[0] <=0) {
		throw_value_error("The planet semi-major axis needs to a positive number");
	}
	if (keplerian_elements[1] < 0 || keplerian_elements[1] >=1) {
		throw_value_error("The planet eccentricity needs to be in [0,1)");
	}
	m_mean_motion = sqrt(mu_central_body / pow(keplerian_elements[0],3));
}

/// Constructor
/**
* Constructs a planet from its position at epoch and its physical parameters
* \param[in] ref_epoch epoch to which the elements are referred to
* \param[in] r0 A STL vector containing the planet position
* \param[in] v0 A STL vector containing the planet velociy
* \param[in] mu_central_body The gravitational parameter of the attracting body (SI units)
* \param[in] mu_self The gravitational parameter of the planet (SI units)
* \param[in] radius radius of the planet (SI units)
* \param[in] safe_radius mimimual radius that is safe during a fly-by of the planet (SI units)
* \param[in] name C++ string containing the planet name. Default value is "Unknown"
*/

keplerian::keplerian(
	const epoch& ref_epoch, 
	const array3D& r0, 
	const array3D& v0, 
	double mu_central_body, 
	double mu_self, 
	double radius, 
	double safe_radius, 
	const std::string &name) : base(mu_central_body, mu_self, radius, safe_radius, name), m_ref_mjd2000(ref_epoch.mjd2000())
{
	ic2par(r0,v0, get_mu_central_body(), m_keplerian_elements);
	m_keplerian_elements[5] = e2m(m_keplerian_elements[5],m_keplerian_elements[1]);
	m_mean_motion = sqrt(get_mu_central_body() / pow(m_keplerian_elements[0],3));
}

/// Polymorphic copy constructor.
planet_ptr keplerian::clone() const
{
	return planet_ptr(new keplerian(*this));
}

void keplerian::eph_impl(double mjd2000, array3D &r, array3D &v) const {
	double elements[6];
	std::copy(m_keplerian_elements.begin(), m_keplerian_elements.end(), elements);
	double dt = (mjd2000 - m_ref_mjd2000) * ASTRO_DAY2SEC;
	elements[5] += m_mean_motion * dt;
	elements[5] = m2e(elements[5],elements[1]);
	par2ic(elements, get_mu_central_body(), r, v);
}

/// Returns the keplerian elements defining the planet
array6D keplerian::get_elements() const {return m_keplerian_elements;}

/// Returns the reference epoch
kep_toolbox::epoch keplerian::get_ref_epoch() const {return kep_toolbox::epoch(m_ref_mjd2000);}

/// Returns the mean motion
double keplerian::get_mean_motion() const {return m_mean_motion;}

/// Sets the keplerian elements (and computes the mean motion)
void keplerian::set_elements(const array6D& el) {
	m_keplerian_elements = el;
	m_mean_motion = sqrt(get_mu_central_body() / pow(m_keplerian_elements[0],3));
}

/// Sets the reference epoch
void keplerian::set_ref_epoch(const kep_toolbox::epoch& when) {
	m_ref_mjd2000 = when.mjd2000();
}

/// Extra informations streamed in humar readable format
std::string keplerian::human_readable_extra() const {
	std::ostringstream s;
	s << "Keplerian planet elements: "<<std::endl;
	s << "Semi major axis (AU): " << boost::lexical_cast<std::string>(m_keplerian_elements[0] / ASTRO_AU) << std::endl;
	s << "Eccentricity: " << boost::lexical_cast<std::string>(m_keplerian_elements[1]) << std::endl;
	s << "Inclination (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[2] * ASTRO_RAD2DEG) << std::endl;
	s << "Big Omega (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[3] * ASTRO_RAD2DEG) << std::endl;
	s << "Small omega (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[4] * ASTRO_RAD2DEG) << std::endl;
	s << "Mean anomaly (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[5] * ASTRO_RAD2DEG) << std::endl;
	s << "Elements reference epoch: " << epoch(m_ref_mjd2000) << std::endl;
	s << "Ephemerides type: Keplerian" << std::endl;
	return s.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::keplerian);
