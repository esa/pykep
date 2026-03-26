// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_SIMS_FLANAGAN_H
#define kep3_LEG_SIMS_FLANAGAN_H

#include <array>
#include <tuple>
#include <vector>

#include <fmt/ostream.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>

namespace kep3::leg
{
/// The Sims-Flanagan leg model
/**
 * This class represents, generically, a low-thrust leg as a sequence of successive
 * impulses of magnitude compatible with the low-thrust propulsion system of a spacecraft.
 * The leg achieves to transfer a given spacecraft from an initial to a final state in the
 * time given (and can be considered as feasible) whenever the method evaluate_mismatch
 * returns all zeros and the method get_throttles_con returns all values less than zero.
 */
class kep3_DLL_PUBLIC sims_flanagan
{
public:
    // Default Constructor.
    sims_flanagan() = default;
    // Constructors
    sims_flanagan(const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
                  const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double veff,
                  double mu, double cut = 0.5);

    // Setters
    void set_tof(double tof);
    void set_rvs(const std::array<std::array<double, 3>, 2> &rv);
    void set_ms(double mass);
    void set_throttles(const std::vector<double> &throttles);
    void set_throttles(const std::vector<double>::const_iterator &it1, const std::vector<double>::const_iterator &it2);
    void set_rvf(const std::array<std::array<double, 3>, 2> &rv);
    void set_mf(double mass);
    void set_max_thrust(double max_thrust);
    void set_veff(double veff);
    void set_mu(double mu);
    void set_cut(double cut);
    void set(const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double veff,
             double mu, double cut = 0.5);

    // Getters
    [[nodiscard]] double get_tof() const;
    [[nodiscard]] const std::array<std::array<double, 3>, 2> &get_rvs() const;
    [[nodiscard]] double get_ms() const;
    [[nodiscard]] const std::vector<double> &get_throttles() const;
    [[nodiscard]] const std::array<std::array<double, 3>, 2> &get_rvf() const;
    [[nodiscard]] double get_mf() const;
    [[nodiscard]] double get_max_thrust() const;
    [[nodiscard]] double get_veff() const;
    [[nodiscard]] double get_mu() const;
    [[nodiscard]] double get_cut() const;
    [[nodiscard]] unsigned get_nseg() const;
    [[nodiscard]] unsigned get_nseg_fwd() const;
    [[nodiscard]] unsigned get_nseg_bck() const;

    // Compute constraints
    [[nodiscard]] std::array<double, 7> compute_mismatch_constraints() const;
    [[nodiscard]] std::vector<double> compute_throttle_constraints() const;

    // Compute mismatch constraint gradients (w.r.t. rvm state and w.r.t. throttles, tof)
    [[nodiscard]] std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
    compute_mc_grad() const;

    // Compute throttle constraint gradients
    [[nodiscard]] std::vector<double> compute_tc_grad() const;

private:
    [[nodiscard]] std::pair<std::array<double, 49>, std::vector<double>>
    gradients_multiple_impulses(std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
                                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                const std::array<std::array<double, 3>, 2> &rvs, double ms, double c, double a,
                                double dt) const;

    [[nodiscard]] std::pair<std::array<double, 49>, std::vector<double>>
    gradients_fwd(std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
                  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                  const std::array<std::array<double, 3>, 2> &rvs, double ms, double c, double a, double dt) const;

    [[nodiscard]] std::pair<std::array<double, 49>, std::vector<double>>
    gradients_bck(std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
                  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                  const std::array<std::array<double, 3>, 2> &rvf_orig, double mf, double c, double a, double dt) const;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & m_rvs;
        ar & m_ms;
        ar & m_throttles;
        ar & m_tof;
        ar & m_rvf;
        ar & m_mf;
        ar & m_max_thrust;
        ar & m_veff;
        ar & m_mu;
        ar & m_cut;
        ar & m_nseg;
        ar & m_nseg_fwd;
        ar & m_nseg_bck;
    }

    // Initial spacecraft state.
    std::array<std::array<double, 3>, 2> m_rvs{{{1., 0., 0.}, {0., 1., 0.}}};
    double m_ms = 1.;
    // Sequence of throttles.
    std::vector<double> m_throttles{0., 0., 0., 0., 0., 0.};
    // Final spacecraft state.
    std::array<std::array<double, 3>, 2> m_rvf{{{0., 1., 0.}, {-1., 0., 0.}}};
    double m_mf = 1.;
    // Time of flight (defaults to 1/4 of the period)
    double m_tof = kep3::pi / 2;
    // Spacecraft propulsion system maximum thrust.
    double m_max_thrust{1.};
    // Spacecraft propulsion system effective velocity.
    double m_veff{1.};
    // Spacecraft gravitational parameter.
    double m_mu{1.};
    // The cut parameter
    double m_cut = 0.5;
    // Segment sizes
    unsigned m_nseg = 2u;
    unsigned m_nseg_fwd = 1u;
    unsigned m_nseg_bck = 1u;
};

// Streaming operator for the class kep3::leg::sims_flanagan.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const sims_flanagan &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::sims_flanagan> : fmt::ostream_formatter {
};

#endif