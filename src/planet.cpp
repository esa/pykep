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

#include "planet.h"
#include"core_functions/ic2par.h"
#include"core_functions/par2ic.h"
#include"core_functions/convert_anomalies.h"
#include"exceptions.h"

namespace kep_toolbox{

planet::planet(const epoch& ref_epoch, const array6D& orbital_elements_, const double & mu_central_body_, const double &mu_self_, const double & radius_, const double & safe_radius_, const std::string &name_)
{
	if (orbital_elements_[0] <=0) {
		throw_value_error("The planet semi-major axis needs to a positive number");
	}
	if (orbital_elements_[1] < 0 || orbital_elements_[1] >=1) {
		throw_value_error("The planet eccentricity needs to be in [0,1)");
	}
	if (radius_ <= 0) {
		throw_value_error("The planet radius needs to be strictly positive");
	}
	if (mu_central_body_ <= 0) {
		throw_value_error("The central body gravitational parameter needs to be strictly positive");
	}
	if (mu_self_ <= 0) {
		throw_value_error("The gravitational parameter of the planet needs to be strictly positive");
	}
	for (int i = 0; i < 3; ++i) {
		cached_r[i] = 0;
		cached_v[i] = 0;
	}
	planet::build_planet(ref_epoch, orbital_elements_, mu_central_body_, mu_self_, radius_, safe_radius_, name_);
}

planet::planet(const epoch& ref_epoch, const array3D& r0, const array3D& v0, const double & mu_central_body_, const double &mu_self_, const double & radius_, const double & safe_radius_, const std::string &name_)
{
	if (radius_ <= 0) {
		throw_value_error("The planet radius needs to be strictly positive");
	}
	if (mu_central_body_ <= 0) {
		throw_value_error("The central body gravitational parameter needs to be strictly positive");
	}
	if (mu_self_ <= 0) {
		throw_value_error("The gravitational parameter of the planet needs to be strictly positive");
	}
	for (int i = 0; i < 3; ++i) {
		cached_r[i] = 0;
		cached_v[i] = 0;
	}
	array6D orbital_elements;
	ic2par(r0,v0, mu_central_body_, orbital_elements);
	orbital_elements[5] = e2m(orbital_elements[5],orbital_elements[1]);
	planet::build_planet(ref_epoch, orbital_elements, mu_central_body_, mu_self_, radius_, safe_radius_, name_);
}


planet_ptr planet::clone() const
{
	return planet_ptr(new planet(*this));
}

void planet::build_planet(const epoch& ref_epoch, const array6D& orbital_elements_, const double & mu_central_body_, const double &mu_self_, const double & radius_, const double & safe_radius_, const std::string &name_)
{
	keplerian_elements = orbital_elements_;
	ref_mjd2000 = ref_epoch.mjd2000();
	radius = radius_;
	safe_radius = safe_radius_;
	mu_self = mu_self_;
	mu_central_body = mu_central_body_;
	m_name = name_;
	mean_motion = sqrt(mu_central_body / pow(keplerian_elements[0],3));
}

void planet::get_eph(const epoch& when, array3D &r, array3D &v) const{
	if(when.mjd2000() != cached_epoch.mjd2000() || cached_epoch.mjd2000() == 0) {
		double elements[6];
		std::copy(keplerian_elements.begin(), keplerian_elements.end(), elements);
		double dt = (when.mjd2000() - ref_mjd2000) * ASTRO_DAY2SEC;
		elements[5] += mean_motion * dt;
		elements[5] = m2e(elements[5],elements[1]);
		par2ic(elements, mu_central_body, cached_r, cached_v);
		cached_epoch = when;
	}
	r = cached_r;
	v = cached_v;
}

array3D planet::get_position(const epoch& when) const {
	array3D r, v;
	get_eph(when, r, v);
	return (r);
}

array3D planet::get_velocity(const epoch& when) const{
	array3D r, v;
	get_eph(when, r, v);
	return (v);
}

array6D planet::get_elements(const epoch& when) const{
	array6D elements(keplerian_elements);
	double dt = (when.mjd2000() - ref_mjd2000) * ASTRO_DAY2SEC;
	elements[5] += mean_motion * dt;
	elements[5] = fmod(elements[5],2*M_PI);
	if (elements[5] < 0) elements[5] = 2*M_PI + elements[5];
	return ( elements );
}

array6D planet::get_elements() const{
	return ( keplerian_elements );
}

/// Destructor.
/**
 * No side effects.
 */
planet::~planet() {}

std::string planet::get_name() const {return m_name;}

double planet::compute_period() const {
	return 2* boost::math::constants::pi<double>() * std::sqrt(std::pow(get_elements()[0],3) / get_mu_central_body());

}

}

/// Overload the stream operator for kep_toolbox::planet
/**
 * Streams out the planet object in a human readable format
 *
 * \param[in] s stream to which the planet will be sent
 * \param[in] body planet to be sent to stream
 *
 * \return reference to s
 *
 */
std::ostream &kep_toolbox::operator<<(std::ostream &s, const kep_toolbox::planet &body) {
	s << "Planet Name: " << body.m_name << std::endl;
	s << "Own gravity parameter: " << boost::lexical_cast<std::string>(body.mu_self) << std::endl;
	s << "Central body gravity parameter: " << boost::lexical_cast<std::string>(body.mu_central_body) << std::endl;
	s << "Planet radius: " << boost::lexical_cast<std::string>(body.radius) << std::endl;
	s << "Planet safe radius: " << boost::lexical_cast<std::string>(body.safe_radius) << std::endl;
	s << "Planet keplerian elements: "<<std::endl;
	array6D elem = body.get_elements(epoch(body.ref_mjd2000));
	s << "Semi major axis (AU): " << boost::lexical_cast<std::string>(elem[0] / ASTRO_AU) << std::endl;
	s << "Eccentricity: " << boost::lexical_cast<std::string>(elem[1]) << std::endl;
	s << "Inclination (deg.): " << boost::lexical_cast<std::string>(elem[2] * ASTRO_RAD2DEG) << std::endl;
	s << "Big Omega (deg.): " << boost::lexical_cast<std::string>(elem[3] * ASTRO_RAD2DEG) << std::endl;
	s << "Small omega (deg.): " << boost::lexical_cast<std::string>(elem[4] * ASTRO_RAD2DEG) << std::endl;
	s << "Mean anomaly (deg.): " << boost::lexical_cast<std::string>(elem[5] * ASTRO_RAD2DEG) << std::endl;
	s << "Elements reference epoch: " << epoch(body.ref_mjd2000) << std::endl;
	return s;
}
