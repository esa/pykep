/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
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

#include "jpl_low_precision.h"
#include <keplerian_toolbox/core_functions/convert_anomalies.hpp>
#include "../core_functions/par2ic.h"
#include <keplerian_toolbox/epoch.hpp>
#include <keplerian_toolbox/exceptions.hpp>

namespace kep_toolbox
{
namespace planet
{

// Data from http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt

static const double mercury_el[6] = {0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593};
static const double mercury_el_dot[6] = {0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081};
static const double venus_el[6] = {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255};
static const double venus_el_dot[6] = {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418};
static const double earth_moon_el[6] = {1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0};
static const double earth_moon_el_dot[6] = {0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0};
static const double mars_el[6] = {1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891};
static const double mars_el_dot[6] = {0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343};
static const double jupiter_el[6] = {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909};
static const double jupiter_el_dot[6] = {-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106};
static const double saturn_el[6] = {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448};
static const double saturn_el_dot[6] = {-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794};
static const double uranus_el[6] = {19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503};
static const double uranus_el_dot[6] = {-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589};
static const double neptune_el[6] = {30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574};
static const double neptune_el_dot[6] = {0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664};
static const double pluto_el[6] = {39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684};
static const double pluto_el_dot[6] = {-0.00031596, 0.00005170, 0.00004818, 145.20780515, -0.04062942, -0.01183482};

/**
 * Construct a planet from its common name (e.g. VENUS)
 *
 * \param[in] name a string describing a planet
 */
jpl_lp::jpl_lp(const std::string &name) : ref_mjd2000(epoch(2451545.0, epoch::JD).mjd2000())
{
    std::map<std::string, int> mapped_planets;
    mapped_planets["mercury"] = 1;
    mapped_planets["venus"] = 2;
    mapped_planets["earth"] = 3;
    mapped_planets["mars"] = 4;
    mapped_planets["jupiter"] = 5;
    mapped_planets["saturn"] = 6;
    mapped_planets["uranus"] = 7;
    mapped_planets["neptune"] = 8;
    mapped_planets["pluto"] = 9;

    double mu_central_body;
    double mu_self;
    double radius;
    double safe_radius;
    std::string lower_case_name = name;
    boost::algorithm::to_lower(lower_case_name);
    switch (mapped_planets[lower_case_name]) {
        case (1): {
            std::copy(mercury_el, mercury_el + 6, &jpl_elements[0]);
            std::copy(mercury_el_dot, mercury_el_dot + 6, &jpl_elements_dot[0]);
            radius = 2440000.;
            safe_radius = 1.1;
            mu_self = 22032e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (2): {
            std::copy(venus_el, venus_el + 6, &jpl_elements[0]);
            std::copy(venus_el_dot, venus_el_dot + 6, &jpl_elements_dot[0]);
            radius = 6052000.;
            safe_radius = 1.1;
            mu_self = 324859e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (3): {
            std::copy(earth_moon_el, earth_moon_el + 6, &jpl_elements[0]);
            std::copy(earth_moon_el_dot, earth_moon_el_dot + 6, &jpl_elements_dot[0]);
            radius = 6378000.;
            safe_radius = 1.1;
            mu_self = 398600.4418e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (4): {
            std::copy(mars_el, mars_el + 6, &jpl_elements[0]);
            std::copy(mars_el_dot, mars_el_dot + 6, &jpl_elements_dot[0]);
            radius = 3397000.;
            safe_radius = 1.1;
            mu_self = 42828e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (5): {
            std::copy(jupiter_el, jupiter_el + 6, &jpl_elements[0]);
            std::copy(jupiter_el_dot, jupiter_el_dot + 6, &jpl_elements_dot[0]);
            radius = 71492000.;
            safe_radius = 9.;
            mu_self = 126686534e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (6): {
            std::copy(saturn_el, saturn_el + 6, &jpl_elements[0]);
            std::copy(saturn_el_dot, saturn_el_dot + 6, &jpl_elements_dot[0]);
            radius = 60330000.;
            safe_radius = 1.1;
            mu_self = 37931187e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (7): {
            std::copy(uranus_el, uranus_el + 6, &jpl_elements[0]);
            std::copy(uranus_el_dot, uranus_el_dot + 6, &jpl_elements_dot[0]);
            radius = 25362000.;
            safe_radius = 1.1;
            mu_self = 5793939e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (8): {
            std::copy(neptune_el, neptune_el + 6, &jpl_elements[0]);
            std::copy(neptune_el_dot, neptune_el_dot + 6, &jpl_elements_dot[0]);
            radius = 24622000.;
            safe_radius = 1.1;
            mu_self = 6836529e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        case (9): {
            std::copy(pluto_el, pluto_el + 6, &jpl_elements[0]);
            std::copy(pluto_el_dot, pluto_el_dot + 6, &jpl_elements_dot[0]);
            radius = 1153000.;
            safe_radius = 1.1;
            mu_self = 871e9;
            mu_central_body = ASTRO_MU_SUN;
        } break;
        default: {
            throw_value_error(std::string("unknown planet name: ") + name);
        }
    }
    set_mu_central_body(mu_central_body);
    set_mu_self(mu_self);
    set_radius(radius);
    set_safe_radius(safe_radius);
    set_name(lower_case_name);
}

/// Polymorphic copy constructor.
planet_ptr jpl_lp::clone() const
{
    return planet_ptr(new jpl_lp(*this));
}

/// Computes the low-precision ephemerides
void jpl_lp::eph_impl(double mjd2000, array3D &r, array3D &v) const
{
    if (mjd2000 <= -73048.0 || mjd2000 >= 18263.0) {
        throw_value_error("Ephemeris are out of range [1800-2050]");
    }
    // algorithm from http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt downloded 2013
    array6D elements, elements2;
    double dt = (mjd2000 - ref_mjd2000) / 36525.0; // Number of centuries passed since J2000.0
    for (unsigned int i = 0; i < 6; ++i) {
        elements[i] = (jpl_elements[i] + jpl_elements_dot[i] * dt);
    }
    elements2[0] = elements[0] * ASTRO_AU;
    elements2[1] = elements[1];
    elements2[2] = elements[2] * ASTRO_DEG2RAD;
    elements2[3] = elements[5] * ASTRO_DEG2RAD;
    elements2[4] = (elements[4] - elements[5]) * ASTRO_DEG2RAD;
    elements2[5] = (elements[3] - elements[4]) * ASTRO_DEG2RAD;
    elements2[5] = m2e(elements2[5], elements2[1]);
    par2ic(elements2, get_mu_central_body(), r, v);
}

/// Extra informations streamed in human readable format
std::string jpl_lp::human_readable_extra() const
{
    std::ostringstream s;
    s << "Ephemerides type: JPL low-precision" << std::endl;
    return s.str();
}
}
} // namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::jpl_lp)
// Serialization code (END)
