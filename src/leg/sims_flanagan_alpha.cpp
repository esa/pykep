// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <fmt/base.h>
#include <iterator>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan_alpha.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg
{

using kep3::linalg::mat13;
using kep3::linalg::mat61;
using kep3::linalg::mat63;
using kep3::linalg::mat66;

// Constructors
sims_flanagan_alpha::sims_flanagan_alpha(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::vector<double> &throttles, const std::vector<double> &talphas,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                             double veff, double mu, double cut)
    : m_rvs(rvs), m_ms(ms), m_throttles(throttles),m_talphas(talphas), m_rvf(rvf), m_mf(mf), m_tof(tof),
      m_max_thrust(max_thrust), m_veff(veff), m_mu(mu), m_cut(cut),
      m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    kep3::leg::_sanity_checks_alpha(m_throttles, m_talphas, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_nseg, m_nseg_fwd, m_nseg_bck);
}

// Setters
void sims_flanagan_alpha::set_tof(double tof)
{
    kep3::leg::_check_tof(tof);
    m_tof = tof;
}
void sims_flanagan_alpha::set_rvs(const std::array<std::array<double, 3>, 2> &rv)
{
    m_rvs = rv;
}
void sims_flanagan_alpha::set_ms(double mass)
{
    m_ms = mass;
}
void sims_flanagan_alpha::set_throttles(const std::vector<double> &throttles)
{
    set_throttles(throttles.begin(), throttles.end());
}
void sims_flanagan_alpha::set_throttles(const std::vector<double>::const_iterator &it1,
                                  const std::vector<double>::const_iterator &it2)
{
    if (((std::distance(it1, it2) % 3) != 0) || std::distance(it1, it2) <= 0) {
        throw std::logic_error("The throttles of a sims_flanagan_alpha leg are being set with invalid iterators.");
    }
    m_throttles.resize(static_cast<size_t>(std::distance(it1, it2)));
    std::copy(it1, it2, m_throttles.begin());
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan_alpha::set_rvf(const std::array<std::array<double, 3>, 2> &rv)
{
    m_rvf = rv;
}
void sims_flanagan_alpha::set_mf(double mass)
{
    m_mf = mass;
}
void sims_flanagan_alpha::set_max_thrust(double max_thrust)
{
    kep3::leg::_check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan_alpha::set_veff(double veff)
{
    kep3::leg::_check_veff(veff);
    m_veff = veff;
}
void sims_flanagan_alpha::set_mu(double mu)
{
    kep3::leg::_check_mu(mu);
    m_mu = mu;
}
void sims_flanagan_alpha::set_cut(double cut)
{
    kep3::leg::_check_cut(cut);
    m_cut = cut;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
// Setter for m_talphas
void sims_flanagan_alpha::set_talphas(const std::vector<double> &talphas)
{
    m_talphas = talphas;
}
void sims_flanagan_alpha::set_talphas(const std::vector<double>::const_iterator &it1,
    const std::vector<double>::const_iterator &it2)
{
    if ( std::distance(it1, it2) <= 0) {
        throw std::logic_error("The talphas of a sims_flanagan_alpha leg are being set with invalid iterators.");
    }
    m_talphas.resize(static_cast<size_t>(std::distance(it1, it2)));
    std::copy(it1, it2, m_talphas.begin());
}


void sims_flanagan_alpha::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                        const std::vector<double> &throttles, const std::vector<double> &talphas,
                        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                        double veff, double mu, double cut)
{
    kep3::leg::_sanity_checks_alpha(m_throttles, m_talphas, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_nseg, m_nseg_fwd, m_nseg_bck);
    m_rvs = rvs;
    m_ms = ms;
    m_throttles = throttles;
    m_talphas = talphas;  // Set the talphas vector
    m_rvf = rvf;
    m_mf = mf;
    m_tof = tof;
    m_max_thrust = max_thrust;
    m_veff = veff;
    m_mu = mu;
    m_cut = cut;
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}

// Getters
double sims_flanagan_alpha::get_tof() const
{
    return m_tof;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan_alpha::get_rvs() const
{
    return m_rvs;
}
double sims_flanagan_alpha::get_ms() const
{
    return m_ms;
}
const std::vector<double> &sims_flanagan_alpha::get_throttles() const
{
    return m_throttles;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan_alpha::get_rvf() const
{
    return m_rvf;
}
double sims_flanagan_alpha::get_mf() const
{
    return m_mf;
}
double sims_flanagan_alpha::get_max_thrust() const
{
    return m_max_thrust;
}
double sims_flanagan_alpha::get_veff() const
{
    return m_veff;
}
double sims_flanagan_alpha::get_mu() const
{
    return m_mu;
}
double sims_flanagan_alpha::get_cut() const
{
    return m_cut;
}
unsigned sims_flanagan_alpha::get_nseg() const
{
    return m_nseg;
}
unsigned sims_flanagan_alpha::get_nseg_fwd() const
{
    return m_nseg_fwd;
}
unsigned sims_flanagan_alpha::get_nseg_bck() const
{
    return m_nseg_bck;
}
// Getter for m_talphas
const std::vector<double> &sims_flanagan_alpha::get_talphas() const
{
    return m_talphas;
}

// The core routines
std::array<double, 7> sims_flanagan_alpha::compute_mismatch_constraints() const
{
    // We introduce some convenience variables
    std::array<double, 3> dv{};
    const double veff = m_veff;
    // const double dtorig = m_tof / static_cast<double>(m_nseg);
    // const double c = m_max_thrust * dt;

    // Find initial and final legs
    double dt0 = (*(m_talphas.begin()));
    double dtf = (*(m_talphas.end() - 1l));
    
    // Initialise intermediate
    double dti = 0.0;
    double dti1 = 0.0;
    double c = 0.0;

    // Forward pass
    // Initial state
    std::array<std::array<double, 3>, 2> rv_fwd(get_rvs());
    double mass_fwd = get_ms();
    // We propagate for a first dt/2 (only if there is at least one forward segment)
    if (m_nseg_fwd > 0) {
        rv_fwd = propagate_lagrangian(rv_fwd, dt0 / 2, m_mu, false).first;
    }
    // We now loop through the forward segments and 1) add a dv + 2) propagate for dt (except on the last segment, where
    // we propagate for dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg_fwd; ++i) {
    // for (auto i = 0u; i < m_nseg_fwd; ++i) {

        // Compute time interval and c
        dti = m_talphas[i];
        dti1 = (i == m_nseg_fwd - 1) ? 0 : m_talphas[i+1];
        c = m_max_thrust * dti;

        // We compute the the dv
        dv[0] = c / mass_fwd * m_throttles[3 * i];
        dv[1] = c / mass_fwd * m_throttles[3 * i + 1];
        dv[2] = c / mass_fwd * m_throttles[3 * i + 2];
        // Add it to the current spacecraft velocity
        rv_fwd[1][0] += dv[0];
        rv_fwd[1][1] += dv[1];
        rv_fwd[1][2] += dv[2];
        // Update the mass accordingly
        const double norm_dv = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
        mass_fwd *= std::exp(-norm_dv / veff);
        // Perform the propagation
        const double prop_duration = (i == m_nseg_fwd - 1) ? dti / 2 : dti/2 + dti1/2;

        // fmt::print("F i {} prop_duration {} dti {} dti1 {} dtorig {}\n",i, prop_duration, dti, dti1, dtorig);
        rv_fwd = propagate_lagrangian(rv_fwd, prop_duration, m_mu, false).first;
    }

    // Backward pass
    // Final state
    std::array<std::array<double, 3>, 2> rv_bck(get_rvf());
    double mass_bck = get_mf();
    // We propagate for a first dt/2 (only if there is at least one backward segment)
    if (m_nseg_bck > 0) {
        rv_bck = propagate_lagrangian(rv_bck, -dtf / 2, m_mu, false).first;
    }
    // We now loop through the backward segments and 1) add a dv + 2) propagate for -dt (except on the last segment,
    // where we propagate for -dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg_bck; ++i) {
    // for (auto i = 0u; i < m_nseg_bck; ++i) {
        
        // Compute time interval and c
        dti = m_talphas[m_talphas.size()-(i+1)];
        dti1 = (i == m_nseg_bck - 1) ? 0 : m_talphas[m_talphas.size()-(i+2)];
        c = m_max_thrust * dti;
        
        // We compute the the dv
        dv[0] = c / mass_bck * m_throttles[m_throttles.size() - 1 - 3 * i - 2];
        dv[1] = c / mass_bck * m_throttles[m_throttles.size() - 1 - 3 * i - 1];
        dv[2] = c / mass_bck * m_throttles[m_throttles.size() - 1 - 3 * i];
        // Subtract it (remember we are going backward) to the current spacecraft velocity
        rv_bck[1][0] -= dv[0];
        rv_bck[1][1] -= dv[1];
        rv_bck[1][2] -= dv[2];
        // Update the mass accordingly (will increase as we go backward)
        double norm_dv = std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
        mass_bck *= std::exp(norm_dv / veff);
        // Perform the propagation
        double prop_duration = (i == m_nseg_bck - 1) ? -dti / 2 : -dti/2 - dti1/2;
        
        // fmt::print("B i {} prop_duration {} dti {} dti1 {} dtorig {}\n",i, prop_duration, dti, dti1, dtorig);
        rv_bck = propagate_lagrangian(rv_bck, prop_duration, m_mu, false).first;
    }

    return {rv_fwd[0][0] - rv_bck[0][0], rv_fwd[0][1] - rv_bck[0][1], rv_fwd[0][2] - rv_bck[0][2],
            rv_fwd[1][0] - rv_bck[1][0], rv_fwd[1][1] - rv_bck[1][1], rv_fwd[1][2] - rv_bck[1][2],
            mass_fwd - mass_bck};
}

std::vector<double> sims_flanagan_alpha::compute_throttle_constraints() const
{
    std::vector<double> retval(m_nseg);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan_alpha &sf)
{
    s << fmt::format("Number of segments: {}\n", sf.get_nseg());
    s << fmt::format("Number of fwd segments: {}\n", sf.get_nseg_fwd());
    s << fmt::format("Number of bck segments: {}\n", sf.get_nseg_bck());
    s << fmt::format("Maximum thrust: {}\n", sf.get_max_thrust());
    s << fmt::format("Central body gravitational parameter: {}\n", sf.get_mu());
    s << fmt::format("Specific impulse: {}\n\n", sf.get_veff());
    s << fmt::format("Time of flight: {}\n", sf.get_tof());
    s << fmt::format("Initial mass: {}\n", sf.get_ms());
    s << fmt::format("Final mass: {}\n", sf.get_mf());
    s << fmt::format("State at departure: {}\n", sf.get_rvs());
    s << fmt::format("State at arrival: {}\n", sf.get_rvf());
    s << fmt::format("Throttles values: {}\n\n", sf.get_throttles());
    s << fmt::format("Mismatch constraints: {}\n", sf.compute_mismatch_constraints());
    s << fmt::format("Throttle constraints: {}\n\n", sf.compute_throttle_constraints());
    return s;
}

} // namespace kep3::leg