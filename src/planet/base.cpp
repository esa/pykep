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

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>

#include <keplerian_toolbox/core_functions/convert_anomalies.hpp>
#include <keplerian_toolbox/core_functions/ic2par.hpp>
#include <keplerian_toolbox/exceptions.hpp>
#include <keplerian_toolbox/planet/base.hpp>

namespace kep_toolbox
{
namespace planet
{

/// Constructor
/** \param[in] mu_central_body gravitational parameter of the central attracting body [m^3/sec^2]
 * \param[in] mu_self gravitational parameter of the body [m^3/sec^2]
 * \param[in] radius body radius [m]
 * \param[in] safe_radius body safe radius [m]
 * \param[in] name body name
*/
base::base(double mu_central_body, double mu_self, double radius, double safe_radius, const std::string &name)
    : m_mu_central_body(mu_central_body), m_mu_self(mu_self), m_radius(radius), m_safe_radius(safe_radius), m_name(name)
{
    if (radius < 0) {
        throw_value_error("The planet radius needs to be positive");
    }
    if (mu_central_body < 0) {
        throw_value_error("The central body gravitational parameter needs to be positive");
    }
    if (mu_self < 0) {
        throw_value_error("The gravitational parameter of the planet needs to be positive");
    }
    if (radius > safe_radius) {
        throw_value_error("Safe radius must be larger than radius");
    }
}

/// Computes the osculating keplerian elements at epoch
/**
 * /param[in] when Epoch at which the osculating elements are computed
 *
 * @return The six orbital elements (a, e, i, W, w, M) in SI units
 */
array6D base::compute_elements(const epoch &when) const
{
    array3D r, v;
    array6D elements;
    this->eph_impl(when.mjd2000(), r, v);
    kep_toolbox::ic2par(r, v, m_mu_central_body, elements);
    elements[5] = kep_toolbox::e2m(elements[5], elements[1]);
    return (elements);
}

/// Gets the planet position and velocity from epoch
/**
* \param[in] when Epoch in which ephemerides are required
* \param[out] r Planet position at epoch (SI units)
* \param[out] v Planet velocity at epoch (SI units)
*/
void base::eph(const epoch &when, array3D &r, array3D &v) const
{
    this->eph_impl(when.mjd2000(), r, v);
}

/// Gets the planet position and velocity from mjd2000
/**
* \param[in]  when mjd2000 in which ephemerides are required
* \param[out] r Planet position at epoch (SI units)
* \param[out] v Planet velocity at epoch (SI units)
*/
void base::eph(const double mjd2000, array3D &r, array3D &v) const
{
    this->eph_impl(mjd2000, r, v);
}

/// Computes the orbital period of the planet at epoch
/**
* \param[in]  when mjd2000 in which ephemerides are required
* \returns The orbital period is seconds
*/
double base::compute_period(const epoch &when) const
{
    return 2 * boost::math::constants::pi<double>()
           * std::sqrt(std::pow(compute_elements(when)[0], 3) / get_mu_central_body());
}

/// Getter for the central body gravitational parameter
/**
 * Gets the gravitational parameter of the central body
 *
 * @return mu_central_body (SI Units)
 */
double base::get_mu_central_body() const
{
    return m_mu_central_body;
};

/// Getter for the planet gravitational parameter
/**
 * Gets the gravitational parameter of the planet
 *
 * @return mu_self (SI Units)
 */
double base::get_mu_self() const
{
    return m_mu_self;
};

/// Getter for the planet radius
/**
 * Gets the radius of the planet
 *
 * @return const reference to radius (SI Units)
 */
double base::get_radius() const
{
    return m_radius;
};

/// Getter for the planet safe-radius
/**
 * Gets the safe-radius of the planet. This is intended to be the minimum distance
 * from the planet center that is safe ... It may be used, for example,  during fly-bys as a constarint
 * on the spacecraft trajectory
 *
 * @return const reference to safe_radius (SI Units)
 */
double base::get_safe_radius() const
{
    return m_safe_radius;
};

/// Returns the planet name
std::string base::get_name() const
{
    return m_name;
};

/// Setter for the planet safe-radius
/**
 * Sets the safe-radius of the planet. This is intended to be the minimum distance
 * from the planet center that is safe ... It is used, for example,  during fly-bys as a constarint
 * on the spacecraft trajectory
 *
 * \param[in] safe_radius Minimum allowed planetary distance (in planetary radius units)
 * \throws value_error if safe_radius in < 1
 */
void base::set_safe_radius(double sr)
{
    if (sr < 1) {
        throw_value_error("Trying to set a safe_radius that is smaller than the planetary radius");
    }
    m_safe_radius = sr * get_radius();
};

/// Setter for the central body gravity parameter
/**
 * Sets the central body gravity parameter
 *
 * \param[in] mu central body gravity parameter
 * \throws value_error if mu is < 0
 */
void base::set_mu_central_body(double mu)
{
    if (mu < 0) {
        throw_value_error("Gravity parameter must be larger than zero");
    }
    m_mu_central_body = mu;
}

/// Setter for the body gravity parameter
/**
 * Sets the  body gravity parameter
 *
 * \param[in] mu body gravity parameter [m^3/sec^2]
 * \throws value_error if mu is < 0
 */
void base::set_mu_self(double mu)
{
    if (mu < 0) {
        throw_value_error("Gravity parameter must be larger than zero");
    }
    m_mu_self = mu;
}

/// Setter for the planet radius
/**
 * Sets the planet radius
 *
 * \param[in] radius planet radius [m]
 * \throws value_error if radius is < 0
 */
void base::set_radius(double radius)
{
    if (radius < 0) {
        throw_value_error("Radius must be larger than zero");
    }
    m_radius = radius;
}

/// Setter for the planet name
/**
 * Sets the planet name
 *
 * \param[in] name planet name
 */
void base::set_name(const std::string &radius)
{
    m_name = radius;
}

/// Human readable output for the planet
std::string base::human_readable() const
{
    std::ostringstream s;
    s << "Planet Name: " << m_name << std::endl;
    s << "Own gravity parameter: " << boost::lexical_cast<std::string>(m_mu_self) << std::endl;
    s << "Central body gravity parameter: " << boost::lexical_cast<std::string>(m_mu_central_body) << std::endl;
    s << "Planet radius: " << boost::lexical_cast<std::string>(m_radius) << std::endl;
    s << "Planet safe radius: " << boost::lexical_cast<std::string>(m_safe_radius) << std::endl;
    s << human_readable_extra();
    return s.str();
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
std::ostream &operator<<(std::ostream &s, const base &body)
{
    s << body.human_readable();
    return s;
}
}
} // Namespaces
