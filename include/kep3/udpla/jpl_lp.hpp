// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_UDPLA_JPL_LP_H
#define kep3_UDPLA_JPL_LP_H

#include <array>

#include <fmt/ostream.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/planet.hpp>

namespace kep3::udpla
{

/// Solar System Planet (jpl simplified ephemerides)
/**
 * This class allows to instantiate planets of
 * the solar system by referring to their common names. The ephemeris used
 * are low_precision ephemeris defined in the JPL pages:
 * https://ssd.jpl.nasa.gov/planets/approx_pos.html and valid for the timeframe
 * 1800AD - 2050 AD
 */

class kep3_DLL_PUBLIC jpl_lp
{

    std::array<double, 6> m_elements;
    std::array<double, 6> m_elements_dot;
    std::string m_name;
    double m_mu_central_body;
    double m_mu_self;
    double m_radius;
    double m_safe_radius;

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar & m_elements;
        ar & m_elements_dot;
        ar & m_name;
        ar & m_mu_central_body;
        ar & m_mu_self;
        ar & m_radius;
        ar & m_safe_radius;
    }

public:
    // Constructor
    explicit jpl_lp(std::string = "earth");
    // Mandatory UDPLA methods
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double) const;

    // Optional UDPLA methods
    [[nodiscard]] std::string get_name() const;
    [[nodiscard]] double get_mu_central_body() const;
    [[nodiscard]] double get_mu_self() const;
    [[nodiscard]] double get_radius() const;
    [[nodiscard]] double get_safe_radius() const;
    [[nodiscard]] std::string get_extra_info() const;
    [[nodiscard]] std::array<double, 6> elements(double = 0., kep3::elements_type = kep3::elements_type::KEP_F) const;

    // JPL does not really define a safe radius for the planets, so we allow the user to change it if needed.
    // This value is mostly used to constrain fly-bys at a safe distance from the planet.
    void set_safe_radius(double);

private:
    [[nodiscard]] std::array<double, 6> _f_elements(double = 0.) const;
};
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const kep3::udpla::jpl_lp &);
} // namespace kep3::udpla

// fmt formatter redirecting to the stream operator
template <>
struct fmt::formatter<kep3::udpla::jpl_lp> : ostream_formatter {
};

// necessary for serialization
KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(kep3::udpla::jpl_lp, kep3::detail::planet_iface)

#endif // KEP_TOOLBOX_PLANET_JPL_LP_H
