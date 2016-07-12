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

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/none.hpp>

#include "j2.h"
#include "../core_functions/ic2par.h"
#include "../core_functions/par2ic.h"
#include "../core_functions/propagate_lagrangian.h"
#include "../core_functions/convert_anomalies.h"
#include "../exceptions.h"

namespace kep_toolbox{ namespace planet {

const array6D j2::default_elements = {{1.0,0.1,0.1,0.1,0.1,0.1}};

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
* \param[in] J2RG2 the product between \f$J_2\f$ and \f$R_G^2\f$ where \f$R_G\f$ is the radius of the oblate primary considered
*/
j2::j2(
	const epoch& ref_epoch,
	const array6D& keplerian_elements,
	double mu_central_body,
	double mu_self,
	double radius,
	double safe_radius,
    double J2RG2,
	const std::string &name
	)
	: base(mu_central_body, mu_self, radius, safe_radius, name), m_keplerian_elements(keplerian_elements), m_ref_mjd2000(ref_epoch.mjd2000()), m_J2RG2(J2RG2)
{
	if (keplerian_elements[0] <=0) {
		throw_value_error("The planet semi-major axis needs to a positive number");
	}
	if (keplerian_elements[1] < 0 || keplerian_elements[1] >=1) {
		throw_value_error("The planet eccentricity needs to be in [0,1)");
	}
	m_mean_motion = sqrt(mu_central_body / pow(keplerian_elements[0],3));
	par2ic(m_keplerian_elements, get_mu_central_body(), m_r, m_v);
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

j2::j2(
	const epoch& ref_epoch,
	const array3D& r0,
	const array3D& v0,
	double mu_central_body,
	double mu_self,
	double radius,
	double safe_radius,
    double J2RG2,
	const std::string &name
    ) : base(mu_central_body, mu_self, radius, safe_radius, name), m_r(r0), m_v(v0), m_ref_mjd2000(ref_epoch.mjd2000()), m_J2RG2(J2RG2)
{
	// This line is  singular (small e and small i) in which case the orbital elements are (simply) not defined
	ic2par(r0,v0, get_mu_central_body(), m_keplerian_elements);
	m_keplerian_elements[5] = e2m(m_keplerian_elements[5],m_keplerian_elements[1]);
	m_mean_motion = sqrt(get_mu_central_body() / pow(m_keplerian_elements[0],3));
}

/// Polymorphic copy constructor.
planet_ptr j2::clone() const
{
	return planet_ptr(new j2(*this));
}

void j2::eph_impl(double mjd2000, array3D &r, array3D &v) const {
	double dt = (mjd2000 - m_ref_mjd2000) * ASTRO_DAY2SEC;
	if (m_keplerian_elements[1] > 1e-5 && m_keplerian_elements[2] > 1e-3)
	{
		double elements[6];
		std::copy(m_keplerian_elements.begin(), m_keplerian_elements.end(), elements);
		elements[5] += m_mean_motion * dt;
		elements[5] = m2e(elements[5],elements[1]);
        double n = std::sqrt(get_mu_central_body() / std::pow(elements[0],3)); // sqrt(mu/a^3)
        double p = elements[0] * (1 - elements[1] * elements[1]); // a(1-e^2)
        double cosi = std::cos(elements[2]);
        double dW = - 3./2. * m_J2RG2 / p / p * n * cosi;
        double dw = 3./4. * m_J2RG2 / p / p * n * (5*cosi*cosi - 1.);
		elements[3] += dW * dt;
		elements[4] += dw * dt;
		par2ic(elements, get_mu_central_body(), r, v);
	} else { // Small inclinations and eccentricities (including nans), we throw directly
		throw_value_error("The planet inclination or eccentricity is too low ... no quick eph computation is avaliable");
	}
}

/// Returns the keplerian elements defining the planet
array6D j2::get_elements() const {return m_keplerian_elements;}

/// Returns the reference epoch
kep_toolbox::epoch j2::get_ref_epoch() const {return kep_toolbox::epoch(m_ref_mjd2000);}

/// Returns the reference epoch in mjd2000
double j2::get_ref_mjd2000() const {return m_ref_mjd2000;}

/// Returns the mean motion
double j2::get_mean_motion() const {return m_mean_motion;}

/// Sets the keplerian elements (and computes the mean motion)
void j2::set_elements(const array6D& el) {
	m_keplerian_elements = el;
	m_mean_motion = sqrt(get_mu_central_body() / pow(m_keplerian_elements[0],3));
}

/// Sets the reference epoch
void j2::set_ref_epoch(const kep_toolbox::epoch& when) {
	m_ref_mjd2000 = when.mjd2000();
}

/// Sets the reference mjd2000
void j2::set_ref_mjd2000(const double &when) {
	m_ref_mjd2000 = when;
}

/// Extra informations streamed in humar readable format
std::string j2::human_readable_extra() const {
	std::ostringstream s;
    s << "Ephemerides type: J2" << "\n\n";
	s << "Orbital elements at epoch: "<<std::endl;
	s << "Semi major axis (AU): " << boost::lexical_cast<std::string>(m_keplerian_elements[0] / ASTRO_AU) << std::endl;
	s << "Eccentricity: " << boost::lexical_cast<std::string>(m_keplerian_elements[1]) << std::endl;
	s << "Inclination (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[2] * ASTRO_RAD2DEG) << std::endl;
	s << "Big Omega (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[3] * ASTRO_RAD2DEG) << std::endl;
	s << "Small omega (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[4] * ASTRO_RAD2DEG) << std::endl;
	s << "Mean anomaly (deg.): " << boost::lexical_cast<std::string>(m_keplerian_elements[5] * ASTRO_RAD2DEG) << std::endl;
	s << "Elements reference epoch: " << epoch(m_ref_mjd2000) << "\n\n";
	s << "J2 RG^2: " << boost::lexical_cast<std::string>(m_J2RG2) << std::endl;
	s << "m_r" << m_r << std::endl;
	s << "m_v" << m_v << std::endl;
	return s.str();
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::j2);
