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

#include"planet_ss.h"
#include"exceptions.h"

namespace kep_toolbox{

planet_ss::planet_ss(const std::string& name)
{
	std::map<std::string, int> mapped_planets;
	mapped_planets["mercury"] = 1; mapped_planets["venus"] = 2; mapped_planets["earth"] = 3;
	mapped_planets["mars"] = 4; mapped_planets["jupiter"] = 5; mapped_planets["saturn"] = 6;
	mapped_planets["uranus"] = 7; mapped_planets["neptune"] = 8;

	double mjd2000 = 0;
	array6D keplerian_elements_;
	double mu_central_body_;
	double mu_self_;
	double radius_;
	double safe_radius_;
	std::string lower_case_name = name;
	boost::algorithm::to_lower(lower_case_name);
	switch ( mapped_planets[lower_case_name] ) {
	case (1): {
			double E[6] = {3.8717591e-01 * ASTRO_AU, 2.0563012e-01,7.0057713e+00 * ASTRO_DEG2RAD,4.8323635e+01 * ASTRO_DEG2RAD, 2.9131238e+01 * ASTRO_DEG2RAD, 1.7274971e+02 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 2440000;
			safe_radius_ = radius_ * 1.1;
			mu_self_ = 22032e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (2): {
			double E[6] = {7.2347372e-01 * ASTRO_AU, 6.757291e-03, 3.3948519e+00 * ASTRO_DEG2RAD, 7.6659746e+01 * ASTRO_DEG2RAD, 5.5219953e+01 * ASTRO_DEG2RAD, 4.9298376e+01 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 6052000;
			safe_radius_ = radius_ * 1.1;
			mu_self_ = 324859e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (3): {
			double E[6] = {9.9998805e-01 * ASTRO_AU, 1.6716812e-02, 8.8543531e-04 * ASTRO_DEG2RAD, 1.7540648e+02 * ASTRO_DEG2RAD, 2.8761578e+02 * ASTRO_DEG2RAD, 2.5760684e+02 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			mjd2000 = epoch(54000.0,epoch::MJD).mjd2000();
			radius_ = 6378000;
			safe_radius_ = radius_*1.1;
			mu_self_ = 398600.4418e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (4): {
			double E[6] = { 1.5239844e+00 * ASTRO_AU, 9.3314935e-02, 1.8506136e+00 * ASTRO_DEG2RAD, 4.9535248e+01 * ASTRO_DEG2RAD, 2.865642e+02 * ASTRO_DEG2RAD, 1.909443e+01  * ASTRO_DEG2RAD };
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 3397000;
			safe_radius_ = radius_*1.1;
			mu_self_ = 42828e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (5): {
			double E[6] = {5.2107645e+00 * ASTRO_AU, 4.9715759e-02, 1.3044197e+00 * ASTRO_DEG2RAD, 1.0044249e+02 * ASTRO_DEG2RAD,2.7550414e+02 * ASTRO_DEG2RAD,1.8387562e+01 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 71492000;
			safe_radius_ = radius_ * 9;
			mu_self_ = 126686534e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (6): {
			double E[6] = {9.5869202e+00 * ASTRO_AU, 5.5944004e-02, 2.4847848e+00 * ASTRO_DEG2RAD, 1.1361884e+02 * ASTRO_DEG2RAD, 3.3583259e+02 * ASTRO_DEG2RAD,-3.9463561e+01 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 60330000;
			safe_radius_ = radius_*1.1;
			mu_self_ = 37931187e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (7): {
			double E[6] = {1.9234105e+01 * ASTRO_AU, 4.4369076e-02, 7.7287008e-01 * ASTRO_DEG2RAD, 7.3908932e+01 * ASTRO_DEG2RAD, 9.6656163e+01 * ASTRO_DEG2RAD, 1.4291587e+02 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 25559000;
			safe_radius_ = radius_*1.1;
			mu_self_ = 5793939e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	case (8): {
			double E[6] = {3.0111359e+01 * ASTRO_AU, 1.1211871e-02, 1.7672166e+00 * ASTRO_DEG2RAD, 1.3176686e+02 * ASTRO_DEG2RAD, 2.6539295e+02 * ASTRO_DEG2RAD, -9.1954294e+01 * ASTRO_DEG2RAD};
			std::copy(E, E + 6, keplerian_elements_.begin());
			radius_ = 24764000;
			safe_radius_ = radius_*1.1;
			mu_self_ = 6836528e9;
			mu_central_body_ = ASTRO_MU_SUN;
		}
		break;
	default : {
		throw_value_error(std::string("unknown planet name") + name);
		}
	}
	build_planet(epoch(mjd2000,epoch::MJD2000),keplerian_elements_,mu_central_body_,mu_self_,radius_,safe_radius_,lower_case_name);
}

planet_ptr planet_ss::clone() const
{
	return planet_ptr(new planet_ss(*this));
}

} //namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet_ss);
// Serialization code (END)
