// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_UDPLA_KEPLERIAN_H
#define kep3_UDPLA_KEPLERIAN_H

#include <array>

#include <fmt/ostream.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/planet.hpp>

namespace kep3::udpla
{

class kep3_DLL_PUBLIC keplerian
{

    kep3::epoch m_ref_epoch;
    std::string m_name;
    double m_mu_central_body;
    double m_mu_self;
    double m_radius;
    double m_safe_radius;
    double m_period;
    bool m_ellipse;
    std::array<std::array<double, 3>, 2> m_pos_vel_0;
    std::array<double, 6> m_kep_f_elements;

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_ref_epoch;
        ar &m_name;
        ar &m_mu_central_body;
        ar &m_mu_self;
        ar &m_radius;
        ar &m_safe_radius;
        ar &m_period;
        ar &m_ellipse;
        ar &m_safe_radius;
        ar &m_pos_vel_0;
        ar &m_kep_f_elements;
    }

public:
    // NOTE: added_param contains mu_self, radius and safe_radius
    explicit keplerian(const epoch &ref_epoch, const std::array<double, 6> &elem, double mu_central_body = 1.,
                       std::string name = "Unknown UDPLA", std::array<double, 3> added_params = {-1., -1., -1.},
                       kep3::elements_type el_t = kep3::elements_type::KEP_F);
    // Constructor from pos_vel
    explicit keplerian(const epoch &ref_epoch = kep3::epoch(0),
                       const std::array<std::array<double, 3>, 2> &pos_vel = {{{1.0, 0.0, 0.0}, {0., 1.0, 0.0}}},
                       double mu_central_body = 1., std::string name = "Unknown UDPLA",
                       std::array<double, 3> added_params = {-1., -1., -1.});
    // Mandatory UDPLA methods
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double) const;

    // Optional UDPLA methods
    [[nodiscard]] std::string get_name() const;
    [[nodiscard]] double get_mu_central_body() const;
    [[nodiscard]] double get_mu_self() const;
    [[nodiscard]] double get_radius() const;
    [[nodiscard]] double get_safe_radius() const;
    [[nodiscard]] std::string get_extra_info() const;
    [[nodiscard]] double period(double =0.) const;

    // Other methods
    [[nodiscard]] kep3::epoch get_ref_epoch() const;
    //[[nodiscard]] std::array<double, 6> elements(double = 0., kep3::elements_type = kep3::elements_type::KEP_F) const;
};
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const kep3::udpla::keplerian &);
} // namespace kep3::udpla

// fmt formatter redirecting to the stream operator
template <>
struct fmt::formatter<kep3::udpla::keplerian> : ostream_formatter {
};

KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(kep3::udpla::keplerian, kep3::detail::planet_iface)

#endif // kep3_UDPLA_KEPLERIAN_H