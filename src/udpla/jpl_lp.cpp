// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <stdexcept>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/jpl_lp.hpp>

namespace kep3::udpla
{

// Data from: https://ssd.jpl.nasa.gov/planets/approx_pos.html
// a,e,i,L,W,w
// clang-format off
static const std::array<double, 6> mercury_el = {0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593};
static const std::array<double, 6> mercury_el_dot = {0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081};
static const std::array<double, 6> venus_el = {0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255};
static const std::array<double, 6> venus_el_dot = {0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418};
static const std::array<double, 6> earth_moon_el = {1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0};
static const std::array<double, 6> earth_moon_el_dot = {0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0};
static const std::array<double, 6> mars_el = {1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891};
static const std::array<double, 6> mars_el_dot = {0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343};
static const std::array<double, 6> jupiter_el = {5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909};
static const std::array<double, 6> jupiter_el_dot = {-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106};
static const std::array<double, 6> saturn_el = {9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448};
static const std::array<double, 6> saturn_el_dot = {-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794};
static const std::array<double, 6> uranus_el = {19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503};
static const std::array<double, 6> uranus_el_dot = {-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589};
static const std::array<double, 6> neptune_el = {30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574};
static const std::array<double, 6> neptune_el_dot = {0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664};
// clang-format on

// NOLINTNEXTLINE(cert-err58-cpp)
const std::unordered_map<std::string, int> mapped_planets
    = {{"mercury", 1}, {"venus", 2},  {"earth", 3},   {"mars", 4}, {"jupiter", 5},
       {"saturn", 6},  {"uranus", 7}, {"neptune", 8}, {"pluto", 9}};

jpl_lp::jpl_lp(std::string name)
    : m_elements(), m_elements_dot(), m_name(std::move(name)), m_mu_central_body(kep3::MU_SUN)
{

    boost::algorithm::to_lower(m_name);
    switch (mapped_planets.at(m_name)) {
        case (1): {
            m_elements = mercury_el;
            m_elements_dot = mercury_el_dot;
            m_radius = 2440000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 22032e9;
        } break;
        case (2): {
            m_elements = venus_el;
            m_elements_dot = venus_el_dot;
            m_radius = 6052000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 324859e9;
        } break;
        case (3): {
            m_elements = earth_moon_el;
            m_elements_dot = earth_moon_el_dot;
            m_radius = 6378000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 398600.4418e9;
        } break;
        case (4): {
            m_elements = mars_el;
            m_elements_dot = mars_el_dot;
            m_radius = 3397000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 42828e9;
        } break;
        case (5): {
            m_elements = jupiter_el;
            m_elements_dot = jupiter_el_dot;
            m_radius = 71492000.;
            m_safe_radius = 9. * m_radius;
            m_mu_self = 126686534e9;
        } break;
        case (6): {
            m_elements = saturn_el;
            m_elements_dot = saturn_el_dot;
            m_radius = 60330000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 37931187e9;
        } break;
        case (7): {
            m_elements = uranus_el;
            m_elements_dot = uranus_el_dot;
            m_radius = 25362000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 5793939e9;
        } break;
        case (8): {
            m_elements = neptune_el;
            m_elements_dot = neptune_el_dot;
            m_radius = 24622000.;
            m_safe_radius = 1.1 * m_radius;
            m_mu_self = 6836529e9;
        } break;
            // LCOV_EXCL_START
        default: {
            throw std::logic_error("unknown planet name: ");
            // LCOV_EXCL_END
        }
    }
    m_name = m_name + "(jpl_lp)";
}

// Computes the kep3::KEP_F elements (osculating with true anomaly) at epoch.
std::array<double, 6> jpl_lp::_f_elements(double mjd2000) const
{
    if (mjd2000 <= -73048.0 || mjd2000 >= 18263.0) {
        throw std::domain_error("Low precision Ephemeris are only valid in the "
                                "range range [1800-2050]");
    }
    // algorithm from https://ssd.jpl.nasa.gov/planets/approx_pos.html accessed
    // 2023.
    std::array<double, 6> elements_updated{}, elements_f{};
    double dt = (mjd2000 - 0.5) / 36525.; // Number of centuries passed since J2000. T = (T_eph - 2451545.0)/36525
    for (unsigned int i = 0; i < 6; ++i) {
        elements_updated[i] = (m_elements[i] + m_elements_dot[i] * dt);
    }
    elements_f[0] = elements_updated[0] * kep3::AU;
    elements_f[1] = elements_updated[1];
    elements_f[2] = elements_updated[2] * kep3::DEG2RAD;
    elements_f[3] = elements_updated[5] * kep3::DEG2RAD;
    elements_f[4] = (elements_updated[4] - elements_updated[5]) * kep3::DEG2RAD;
    elements_f[5] = (elements_updated[3] - elements_updated[4]) * kep3::DEG2RAD;
    elements_f[5] = kep3::m2f(elements_f[5], elements_f[1]);
    return elements_f;
}

std::array<std::array<double, 3>, 2> jpl_lp::eph(double mjd2000) const
{
    auto elements_f = _f_elements(mjd2000);
    return par2ic(elements_f, get_mu_central_body());
}

std::array<double, 6> jpl_lp::elements(double mjd2000, kep3::elements_type el_type) const
{
    auto elements = _f_elements(mjd2000);
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            break;
        case kep3::elements_type::KEP_M:
            elements[5] = kep3::f2m(elements[5], elements[1]);
            break;
        case kep3::elements_type::MEE:
            elements = kep3::par2mee(elements, false);
            break;
        case kep3::elements_type::MEE_R:
            elements = kep3::par2mee(elements, true);
            break;
        default:
            throw std::logic_error("You should not go here!");
    }
    return elements;
}

std::string jpl_lp::get_name() const
{
    return m_name;
}

double jpl_lp::get_mu_central_body() const
{
    return m_mu_central_body;
}

double jpl_lp::get_mu_self() const
{
    return m_mu_self;
}

double jpl_lp::get_radius() const
{
    return m_radius;
}

double jpl_lp::get_safe_radius() const
{
    return m_safe_radius;
}

std::string jpl_lp::get_extra_info() const
{
    auto par = elements(0.);
    auto pos_vel = eph(0.);

    std::string retval = fmt::format("\nLow-precision planet elements: \n")
                         + fmt::format("Semi major axis (AU): {}\n", par[0] / kep3::AU)
                         + fmt::format("Eccentricity: {}\n", par[1])
                         + fmt::format("Inclination (deg.): {}\n", par[2] * kep3::RAD2DEG)
                         + fmt::format("Big Omega (deg.): {}\n", par[3] * kep3::RAD2DEG)
                         + fmt::format("Small omega (deg.): {}\n", par[4] * kep3::RAD2DEG)
                         + fmt::format("True anomly (deg.): {}\n", par[5] * kep3::RAD2DEG);
    retval += fmt::format("Mean anomly (deg.): {}\n", kep3::f2m(par[5], par[1]) * kep3::RAD2DEG);
    retval += fmt::format("Elements reference epoch (MJD2000): {}\n", 0.)
              + fmt::format("Elements reference epoch (date): {}\n", kep3::epoch(0.))
              + fmt::format("r at ref. = {}\n", pos_vel[0]) + fmt::format("v at ref. = {}\n", pos_vel[1]);
    return retval;
}

void jpl_lp::set_safe_radius(double safe_radius)
{
    if (safe_radius < m_radius) {
        throw std::domain_error("The udpla safe radius must be at least the planet's radius ("
                                + std::to_string(m_radius) + "), while it is: " + std::to_string(safe_radius));
    }
    m_safe_radius = safe_radius;
}

std::ostream &operator<<(std::ostream &os, const kep3::udpla::jpl_lp &udpla)
{
    os << udpla.get_extra_info() << "\n";
    return os;
}

} // namespace kep3::udpla

// NOLINTNEXTLINE
KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(kep3::udpla::jpl_lp, kep3::detail::planet_iface)