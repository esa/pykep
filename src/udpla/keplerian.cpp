// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cmath>
#include <limits>

#include <chrono>
#include <stdexcept>
#include <string>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>
#include <kep3/udpla/keplerian.hpp>

namespace kep3::udpla
{

keplerian::keplerian(const epoch &ref_epoch, const std::array<std::array<double, 3>, 2> &pos_vel,
                     double mu_central_body, std::string name, std::array<double, 3> added_params)
    : m_ref_epoch(ref_epoch), m_name(std::move(name)), m_mu_central_body(mu_central_body), m_mu_self(added_params[0]),
      m_radius(added_params[1]), m_safe_radius(added_params[2]), m_period(), m_ellipse(), m_pos_vel_0(pos_vel),
      m_kep_f_elements()
{
    const double R = std::sqrt(pos_vel[0][0] * pos_vel[0][0] + pos_vel[0][1] * pos_vel[0][1] + pos_vel[0][2] * pos_vel[0][2]);
    const double en = (pos_vel[1][0] * pos_vel[1][0] + pos_vel[1][1] * pos_vel[1][1] + pos_vel[1][2] * pos_vel[1][2]) / 2.
                - mu_central_body / R;
    (en > 0) ? m_ellipse = false : m_ellipse = true;
    const double a = -m_mu_central_body / 2. / en;
    if (m_ellipse) {
        m_period = kep3::pi * 2. * std::sqrt(a * a * a / m_mu_central_body);
    } else {
        m_period = std::numeric_limits<double>::quiet_NaN();
    }
    m_kep_f_elements = kep3::ic2par(pos_vel, mu_central_body);
}

keplerian::keplerian(const epoch &ref_epoch, const std::array<double, 6> &elem, double mu_central_body,
                     std::string name, std::array<double, 3> added_params, kep3::elements_type el_type)
    : m_ref_epoch(ref_epoch), m_name(std::move(name)), m_mu_central_body(mu_central_body), m_mu_self(added_params[0]),
      m_radius(added_params[1]), m_safe_radius(added_params[2]), m_period(), m_ellipse(), m_pos_vel_0(),
      m_kep_f_elements(elem)
{
    // We check if the input is consistent with the pykep convention of a,e
    if (m_kep_f_elements[0] * (1 - m_kep_f_elements[1]) <= 0) {
        throw std::domain_error("A Keplerian planet constructor was called with with non compatible "
                                "a,e:"
                                "The following must hold: a<0 -> e>1 [a>0 -> e<1].");
    }

    // If so, we convert to the kep_f elements according to the chosen type for the input elem
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            break;
        case kep3::elements_type::KEP_M:
            if (m_kep_f_elements[0] < 0) {
                throw std::logic_error("Mean anomaly is only available for ellipses.");
            }
            m_kep_f_elements[5] = kep3::m2f(m_kep_f_elements[5], m_kep_f_elements[1]);
            break;
        case kep3::elements_type::MEE:
            m_kep_f_elements = kep3::mee2par(m_kep_f_elements, false);
            break;
        case kep3::elements_type::MEE_R:
            m_kep_f_elements = kep3::mee2par(m_kep_f_elements, true);
            break;
        default:
            throw std::logic_error("You should not go here!");
    }

    m_pos_vel_0 = kep3::par2ic(m_kep_f_elements, mu_central_body);

    (m_kep_f_elements[0] < 0) ? m_ellipse = false : m_ellipse = true;
    if (m_ellipse) {
        m_period = kep3::pi * 2.
                   * std::sqrt(m_kep_f_elements[0] * m_kep_f_elements[0] * m_kep_f_elements[0] / m_mu_central_body);
    } else {
        m_period = std::numeric_limits<double>::quiet_NaN();
    }
}

// This is the mandatory method for the planet interface
std::array<std::array<double, 3>, 2> keplerian::eph(double mjd2000) const
{
    // 1 - We compute the dt
    auto dt = (mjd2000 - m_ref_epoch.mjd2000()) * kep3::DAY2SEC;
    // 2 - We propagate
    auto retval = kep3::propagate_lagrangian(m_pos_vel_0, dt, m_mu_central_body).first;
    return retval;
}

std::string keplerian::get_name() const
{
    return m_name;
}

double keplerian::get_mu_central_body() const
{
    return m_mu_central_body;
}

double keplerian::get_mu_self() const
{
    return m_mu_self;
}

double keplerian::get_radius() const
{
    return m_radius;
}

double keplerian::get_safe_radius() const
{
    return m_safe_radius;
}

double keplerian::period(double) const
{
    return m_period;
}

kep3::epoch keplerian::get_ref_epoch() const
{
    return m_ref_epoch;
}

//std::array<double, 6> keplerian::elements(double, kep3::elements_type el_type) const
//{
//    std::array<double, 6> retval{};
//    switch (el_type) {
//        case kep3::elements_type::KEP_F:
//            retval = m_kep_f_elements;
//            break;
//        case kep3::elements_type::KEP_M:
//            if (!m_ellipse) {
//                throw std::logic_error("Mean anomaly is only available for ellipses.");
//            }
//            retval = m_kep_f_elements;
//            retval[5] = kep3::f2m(retval[5], retval[1]);
//            break;
//        case kep3::elements_type::MEE:
//            retval = m_kep_f_elements;
//            retval = kep3::par2mee(retval, false);
//            break;
//        case kep3::elements_type::MEE_R:
//            retval = m_kep_f_elements;
//            retval = kep3::par2mee(retval, true);
//            break;
//        default:
//            throw std::logic_error("You should not go here!");
//    }
//    return retval;
//}

std::string keplerian::get_extra_info() const
{
    std::string retval = fmt::format("Keplerian planet elements: \n")
                         + fmt::format("Semi major axis: {}\n", m_kep_f_elements[0])
                         + fmt::format("Semi major axis (AU): {}\n", m_kep_f_elements[0] / kep3::AU)
                         + fmt::format("Eccentricity: {}\n", m_kep_f_elements[1])
                         + fmt::format("Inclination (deg.): {}\n", m_kep_f_elements[2] * kep3::RAD2DEG)
                         + fmt::format("Big Omega (deg.): {}\n", m_kep_f_elements[3] * kep3::RAD2DEG)
                         + fmt::format("Small omega (deg.): {}\n", m_kep_f_elements[4] * kep3::RAD2DEG)
                         + fmt::format("True anomly (deg.): {}\n", m_kep_f_elements[5] * kep3::RAD2DEG);
    if (m_ellipse) {
        retval += fmt::format("Mean anomly (deg.): {}\n",
                              kep3::f2m(m_kep_f_elements[5], m_kep_f_elements[1]) * kep3::RAD2DEG);
    }
    retval += fmt::format("Elements reference epoch (MJD2000): {}\n", m_ref_epoch.mjd2000())
              + fmt::format("Elements reference epoch (UTC): {}\n", m_ref_epoch)
              + fmt::format("r at ref. = {}\n", m_pos_vel_0[0]) + fmt::format("v at ref. = {}\n", m_pos_vel_0[1]);
    return retval;
}

std::ostream &operator<<(std::ostream &os, const kep3::udpla::keplerian &udpla)
{
    os << udpla.get_extra_info() << "\n";
    return os;
}

} // namespace kep3::udpla

// NOLINTNEXTLINE
KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(kep3::udpla::keplerian, kep3::detail::planet_iface)
