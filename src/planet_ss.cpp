/*****************************************************************************
 *   Copyright (C) 2004-2012 The PyKEP development team,                     *
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
            double E[6] = {3.772942436743051E-01 * ASTRO_AU, 2.200014582454420E-01, 7.072030953380287E+00 * ASTRO_DEG2RAD, 4.754470860361798E+01 * ASTRO_DEG2RAD, 2.767991252328971E+01 * ASTRO_DEG2RAD, 1.899912269338703E+02 * ASTRO_DEG2RAD};
            std::copy(E, E + 6, keplerian_elements_.begin());
            mjd2000 = epoch(2458849.5,epoch::JD).mjd2000();
            radius_ = 2440000;
            safe_radius_ = radius_ * 1.1;
            mu_self_ = 22032e9;
            mu_central_body_ = ASTRO_MU_SUN;
        }
        break;
    case (2): {
            double E[6] = {7.159637980420038E-01 * ASTRO_AU, 1.103010682131625E-02, 3.416850348494068E+00 * ASTRO_DEG2RAD, 7.675615379502699E+01 * ASTRO_DEG2RAD, 1.337130135915192E+02 * ASTRO_DEG2RAD, 1.537146531316461E+02 * ASTRO_DEG2RAD};
            std::copy(E, E + 6, keplerian_elements_.begin());
            mjd2000 = epoch(2458849.5,epoch::JD).mjd2000();
            radius_ = 6052000;
            safe_radius_ = radius_ * 1.1;
            mu_self_ = 324859e9;
            mu_central_body_ = ASTRO_MU_SUN;
        }
        break;
    case (3): {
            double E[6] = {1.016281672453605E+00 * ASTRO_AU, 2.467020792109220E-02, 1.812249943317730E-03 * ASTRO_DEG2RAD, 1.420066138465537E+02 * ASTRO_DEG2RAD, 3.142255714626496E+02 * ASTRO_DEG2RAD, 3.472751672370176E+00 * ASTRO_DEG2RAD};
            std::copy(E, E + 6, keplerian_elements_.begin());
            mjd2000 = epoch(2458849.5,epoch::JD).mjd2000();
            radius_ = 6378000;
            safe_radius_ = radius_*1.1;
            mu_self_ = 398600.4418e9;
            mu_central_body_ = ASTRO_MU_SUN;
        }
        break;
    case (4): {
            double E[6] = { 1.519223936334168E+00 * ASTRO_AU, 9.868757921143217E-02, 1.846247867207550E+00 * ASTRO_DEG2RAD, 4.925441005680170E+01 * ASTRO_DEG2RAD, 2.868638400110784E+02 * ASTRO_DEG2RAD, 2.473483982901654E+02  * ASTRO_DEG2RAD };
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
