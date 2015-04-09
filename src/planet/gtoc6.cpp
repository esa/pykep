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

#include "gtoc6.h"
#include "../exceptions.h"

namespace kep_toolbox{ namespace planet {

/**
 * Construct a Jupiter moon from its common name
 * \param[in] name a string describing a planet
 */
gtoc6::gtoc6(const std::string& name)
{
	std::map<std::string, int> mapped_planets;
	mapped_planets["io"] = 1; mapped_planets["europa"] = 2; mapped_planets["ganymede"] = 3;
	mapped_planets["callisto"] = 4;
	double mjd = 58849.0;
	array6D keplerian_elements_;
	double mu_central_body_;
	double mu_self_;
	double radius_;
	double safe_radius_;
	const double mu_jupiter =  126686534921800000.0; //m^3/s^2
	std::string lower_case_name = name;
	boost::algorithm::to_lower(lower_case_name);
	switch ( mapped_planets[lower_case_name] ) {
	case (1): {
			double E[6] = {422029687.14001,  4.308524661773E-03,  40.11548686966E-03 * ASTRO_DEG2RAD,   -79.640061742992 * ASTRO_DEG2RAD,    37.991267683987 * ASTRO_DEG2RAD,  286.85240405645 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 1826500;
			safe_radius_ = radius_ + 50000.0;
			mu_self_ = 5959916000000.0;
			mu_central_body_ = mu_jupiter;
		}
		break;
	case (2): {
			double E[6] = {671224237.12681, 9.384699662601E-03, 0.46530284284480 * ASTRO_DEG2RAD, -132.15817268686 * ASTRO_DEG2RAD, -79.571640035051 * ASTRO_DEG2RAD, 318.00776678240 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 1561000.0;
			safe_radius_ = radius_ + 50000.0;
			mu_self_ = 3202739000000.0;
			mu_central_body_ = mu_jupiter;
		}
		break;
	case (3): {
			double E[6] = {1070587469.2374000, 1.953365822716E-03, 0.13543966756582 * ASTRO_DEG2RAD, -50.793372416917 * ASTRO_DEG2RAD, -42.876495018307 * ASTRO_DEG2RAD, 220.59841030407 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 2634000.0;
			safe_radius_ = radius_+ 50000.0;
			mu_self_ = 9887834000000.0;
			mu_central_body_ = mu_jupiter;
		}
		break;
	case (4): {
			double E[6] = { 1883136616.7305, 7.337063799028E-03, 0.25354332731555 * ASTRO_DEG2RAD, 86.723916616548 * ASTRO_DEG2RAD, -160.76003434076 * ASTRO_DEG2RAD, 321.07650614246 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 2408000;
			safe_radius_ = radius_ + 50000.0;
			mu_self_ = 7179289000000.0;
			mu_central_body_ = mu_jupiter;
		}
		break;

	default : {
		throw_value_error(std::string("unknown planet name") + name);
		}
	}
	set_mu_central_body(mu_central_body_);
	set_mu_self(mu_self_);
	set_radius(radius_);
	set_safe_radius(safe_radius_ / radius_);
	set_name(lower_case_name);
	set_elements(keplerian_elements_);
	set_ref_epoch(epoch(mjd,epoch::MJD));
}

planet_ptr gtoc6::clone() const
{
	return planet_ptr(new gtoc6(*this));
}

}} //namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::gtoc6)
// Serialization code (END)
