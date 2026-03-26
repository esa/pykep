// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>
#include <limits>

#include <boost/core/demangle.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/epoch.hpp>
#include <kep3/exceptions.hpp>
#include <kep3/planet.hpp>

namespace kep3::detail
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double period_from_energy(const std::array<double, 3> &r, const std::array<double, 3> &v, double mu)
{
    double const R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double const v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double const en = v2 / 2. - mu / R;
    if (en > 0) {
        // If the energy is positive we have an hyperbolae and we return nan
        return std::numeric_limits<double>::quiet_NaN();
    } else {
        double const a = -mu / 2. / en;
        return kep3::pi * 2. * std::sqrt(a * a * a / mu);
    }
}

std::array<double, 6> elements_from_posvel(const std::array<std::array<double, 3>, 2> &pos_vel, double mu,
                                           kep3::elements_type el_type)
{
    std::array<double, 6> retval;
    switch (el_type) {
        case kep3::elements_type::KEP_F:
            retval = kep3::ic2par(pos_vel, mu);
            break;
        case kep3::elements_type::KEP_M:
            retval = kep3::ic2par(pos_vel, mu);
            if (retval[0] < 1) {
                throw std::logic_error("Mean anomaly is only available for ellipses.");
            }
            retval[5] = kep3::f2m(retval[5], retval[1]);
            break;
        case kep3::elements_type::MEE:
            retval = kep3::ic2mee(pos_vel, mu, false);
            break;
        case kep3::elements_type::MEE_R:
            retval = kep3::ic2mee(pos_vel, mu, true);
            break;
        // LCOV_EXCL_START
        default:
            throw std::logic_error("You should not go here!");
            // LCOV_EXCL_END
    }
    return retval;
}

std::array<std::array<double, 3>, 2> null_udpla::eph(double)
{
    std::array<double, 3> const pos = {1., 0., 0.};
    std::array<double, 3> const vel = {0., 1., 0.};
    return {pos, vel};
}

planet_iface::~planet_iface() = default;

std::array<std::array<double, 3>, 2> planet_iface::eph(const epoch &ep) const
{
    return eph(ep.mjd2000());
}

std::array<double, 3> planet_iface::acc(const epoch &ep) const
{
    return acc(ep.mjd2000());
}

std::array<double, 6> planet_iface::elements(const epoch &ep, kep3::elements_type el_type) const
{
    return elements(ep.mjd2000(), el_type);
}

double planet_iface::period(const epoch &ep) const
{
    return period(ep.mjd2000());
}

} // namespace kep3::detail

namespace kep3
{
std::ostream &operator<<(std::ostream &os, const planet &p)
{
    os << "Planet name: " << p.get_name();
    os << "\nC++ class name: " << boost::core::demangle(p.get_type_index().name()) << '\n';
    os << fmt::format("\nmu central body (-1 if not defined): {}\n", p.get_mu_central_body());
    os << fmt::format("mu body (-1 if not defined): {}\n", p.get_mu_self());
    os << fmt::format("radius body (-1 if not defined): {}\n", p.get_radius());
    os << fmt::format("safe body radius (-1 if not defined): {}\n", p.get_safe_radius());

    const auto extra_str = p.get_extra_info();
    if (!extra_str.empty()) {
        os << "\nExtra info:\n" << extra_str;
    }
    return os;
}
} // namespace kep3

// NOLINTNEXTLINE
KEP3_S11N_EXPORT_IMPLEMENT_AND_INSTANTIATE(kep3::detail::null_udpla, kep3::detail::planet_iface)