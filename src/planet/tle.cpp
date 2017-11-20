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

#include <cmath>
#include <string>

#include "../core_functions/convert_anomalies.h"
#include "../core_functions/par2ic.h"
#include "../exceptions.h"
#include "../third_party/libsgp4/Eci.h"
#include "../third_party/libsgp4/Globals.h"
#include "../third_party/libsgp4/SGP4.h"
#include "../third_party/libsgp4/SatelliteException.h"
#include "../third_party/libsgp4/Tle.h"
#include "../third_party/libsgp4/TleException.h"
#include "tle.h"

namespace kep_toolbox
{
namespace planet
{

/**
 * Construct a planet_tle from two strings containing the two line elements
 * \param[in] line1 first line
 * \param[in] line2 second line
 */
tle::tle(const std::string &line1, const std::string &line2) try : base(),
                                                                   m_line1(line1),
                                                                   m_line2(line2),
                                                                   m_tle(Tle("TLE satellite", line1, line2)),
                                                                   m_sgp4_propagator(SGP4(m_tle)) {
    // We read the osculating elements of the satellite
    array6D keplerian_elements;
    double mu_central_body = kMU * 1E09;                                                     // (m^3/s^2)
    double mean_motion = m_tle.MeanMotion() * 2 * kPI / kSECONDS_PER_DAY;                    // [rad/s]
    keplerian_elements[0] = std::pow(mu_central_body / (mean_motion * mean_motion), 1. / 3); // a [m]
    keplerian_elements[1] = m_tle.Eccentricity();                                            // e
    keplerian_elements[2] = m_tle.Inclination(false);                                        // i [rad]
    keplerian_elements[3] = m_tle.RightAscendingNode(false);                                 // Om [rad]
    keplerian_elements[4] = m_tle.ArgumentPerigee(false);                                    // om [rad]
    keplerian_elements[5] = m_tle.MeanAnomaly(false);                                        // M [rad]

    std::string year_str = m_tle.IntDesignator().substr(0, 2);
    int year_int = std::stoi(year_str);
    std::string rest = m_tle.IntDesignator().substr(2);
    int prefix = (year_int > 56) ? (19) : (20);

    std::string object_name(std::to_string(prefix) + year_str + std::string("-") + rest);
    set_mu_central_body(mu_central_body);
    set_name(object_name);
    m_ref_mjd2000 = epoch(m_tle.Epoch().ToJulian(), epoch::JD).mjd2000();

} catch (TleException &e) {
    // std::cout << "TleException cought in planet_tle constructor" << std::endl;
    throw_value_error(e.what());
} catch (SatelliteException &e) {
    // std::cout << "SatelliteException cought in planet_tle constructor" << std::endl;
    throw_value_error(e.what());
}

/// Polymorphic copy constructor.le::clone() const
planet_ptr tle::clone() const
{
    return planet_ptr(new tle(*this));
}

void tle::eph_impl(double mjd2000, array3D &r, array3D &v) const
{
    Vector position;
    Vector velocity;
    double minutes_since = (mjd2000 - m_ref_mjd2000) * 24 * 60;

    try {
        Eci eci = m_sgp4_propagator.FindPosition(minutes_since);
        position = eci.Position();
        velocity = eci.Velocity();
        r[0] = position.x * 1000;
        r[1] = position.y * 1000;
        r[2] = position.z * 1000;
        v[0] = velocity.x * 1000;
        v[1] = velocity.y * 1000;
        v[2] = velocity.z * 1000;
    } catch (SatelliteException &e) {
        // std::cout << "SatelliteException caught while computing ephemerides" << std::endl;
        throw_value_error(e.what());
    } catch (DecayedException &e) {
        // std::cout << "DecayedException caught while computing ephemerides" << std::endl;
        throw_value_error(e.what());
    }
}

/// Getter for the reference mjd2000
double tle::get_ref_mjd2000() const
{
    return m_ref_mjd2000;
}

/// Getter for TLE line1
std::string tle::get_line1() const
{
    return m_line1;
}

/// Getter for TLE line2
std::string tle::get_line2() const
{
    return m_line2;
}

/// Setter for the epoch of the TLE (workaround for 2056/2057 bug)
void tle::set_epoch(const unsigned int year, const double day)
{
    m_tle.setEpoch(year, day);
    m_sgp4_propagator.SetTle(m_tle);
    m_ref_mjd2000 = epoch(m_tle.Epoch().ToJulian(), epoch::JD).mjd2000();
}

/// Extra informations streamed in human readable format
std::string tle::human_readable_extra() const
{
    std::ostringstream s;
    s << "Ephemerides type: SGP4 propagator" << std::endl;
    s << "TLE epoch: " << epoch(m_ref_mjd2000, epoch::MJD2000) << std::endl;
    s << "TLE 1: " << m_line1 << std::endl;
    s << "TLE 2: " << m_line2 << std::endl;
    return s.str();
}
}
} // namespace

BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::tle)
