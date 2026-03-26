// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cstddef>
#include <iterator>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/leg/sf_checks.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg
{

using kep3::linalg::_dot;
using kep3::linalg::mat13;
using kep3::linalg::mat61;
using kep3::linalg::mat63;
using kep3::linalg::mat66;

// Constructors
sims_flanagan::sims_flanagan(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                             const std::vector<double> &throttles,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                             double veff, double mu, double cut)
    : m_rvs(rvs), m_ms(ms), m_throttles(throttles), m_rvf(rvf), m_mf(mf), m_tof(tof), m_max_thrust(max_thrust),
      m_veff(veff), m_mu(mu), m_cut(cut), m_nseg(static_cast<unsigned>(m_throttles.size()) / 3u),
      m_nseg_fwd(static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut)), m_nseg_bck(m_nseg - m_nseg_fwd)
{
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_nseg, m_nseg_fwd, m_nseg_bck);
}

// Setters
void sims_flanagan::set_tof(double tof)
{
    kep3::leg::_check_tof(tof);
    m_tof = tof;
}
void sims_flanagan::set_rvs(const std::array<std::array<double, 3>, 2> &rv)
{
    m_rvs = rv;
}
void sims_flanagan::set_ms(double mass)
{
    m_ms = mass;
}

void sims_flanagan::set_throttles(const std::vector<double>::const_iterator &it1,
                                  const std::vector<double>::const_iterator &it2)
{
    if (((std::distance(it1, it2) % 3) != 0) || std::distance(it1, it2) <= 0) {
        throw std::logic_error("The throttles of a sims_flanagan leg are being set with invalid iterators.");
    }
    m_throttles.resize(static_cast<size_t>(std::distance(it1, it2)));
    std::copy(it1, it2, m_throttles.begin());
    m_nseg = static_cast<unsigned>(m_throttles.size()) / 3u;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan::set_throttles(const std::vector<double> &throttles)
{
    set_throttles(throttles.begin(), throttles.end());
}

void sims_flanagan::set_rvf(const std::array<std::array<double, 3>, 2> &rv)
{
    m_rvf = rv;
}
void sims_flanagan::set_mf(double mass)
{
    m_mf = mass;
}
void sims_flanagan::set_max_thrust(double max_thrust)
{
    kep3::leg::_check_max_thrust(max_thrust);
    m_max_thrust = max_thrust;
}
void sims_flanagan::set_veff(double veff)
{
    kep3::leg::_check_veff(veff);
    m_veff = veff;
}
void sims_flanagan::set_mu(double mu)
{
    kep3::leg::_check_mu(mu);
    m_mu = mu;
}
void sims_flanagan::set_cut(double cut)
{
    kep3::leg::_check_cut(cut);
    m_cut = cut;
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}
void sims_flanagan::set(const std::array<std::array<double, 3>, 2> &rvs, double ms,
                        const std::vector<double> &throttles,
                        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                        const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust,
                        double veff, double mu, double cut)
{
    kep3::leg::_sanity_checks(m_throttles, m_tof, m_max_thrust, m_veff, m_mu, m_cut, m_nseg, m_nseg_fwd, m_nseg_bck);
    m_rvs = rvs;
    m_ms = ms;
    m_throttles = throttles;
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
double sims_flanagan::get_tof() const
{
    return m_tof;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan::get_rvs() const
{
    return m_rvs;
}
double sims_flanagan::get_ms() const
{
    return m_ms;
}
const std::vector<double> &sims_flanagan::get_throttles() const
{
    return m_throttles;
}
const std::array<std::array<double, 3>, 2> &sims_flanagan::get_rvf() const
{
    return m_rvf;
}
double sims_flanagan::get_mf() const
{
    return m_mf;
}
double sims_flanagan::get_max_thrust() const
{
    return m_max_thrust;
}
double sims_flanagan::get_veff() const
{
    return m_veff;
}
double sims_flanagan::get_mu() const
{
    return m_mu;
}
double sims_flanagan::get_cut() const
{
    return m_cut;
}
unsigned sims_flanagan::get_nseg() const
{
    return m_nseg;
}
unsigned sims_flanagan::get_nseg_fwd() const
{
    return m_nseg_fwd;
}
unsigned sims_flanagan::get_nseg_bck() const
{
    return m_nseg_bck;
}

// The core routines
std::array<double, 7> sims_flanagan::compute_mismatch_constraints() const
{
    // We introduce some convenience variables
    std::array<double, 3> dv{};
    const double veff = m_veff;
    const double dt = m_tof / static_cast<double>(m_nseg);
    const double c = m_max_thrust * dt;
    // Forward pass
    // Initial state
    std::array<std::array<double, 3>, 2> rv_fwd(get_rvs());
    double mass_fwd = get_ms();
    // We propagate for a first dt/2 (only if there is at least one forward segment)
    if (m_nseg_fwd > 0) {
        rv_fwd = propagate_lagrangian(rv_fwd, dt / 2, m_mu, false).first;
    }
    // We now loop through the forward segments and 1) add a dv + 2) propagate for dt (except on the last segment, where
    // we propagate for dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg_fwd; ++i) {
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
        const double prop_duration = (i == m_nseg_fwd - 1) ? dt / 2 : dt;
        rv_fwd = propagate_lagrangian(rv_fwd, prop_duration, m_mu, false).first;
    }

    // Backward pass
    // Final state
    std::array<std::array<double, 3>, 2> rv_bck(get_rvf());
    double mass_bck = get_mf();
    // We propagate for a first dt/2 (only if there is at least one backward segment)
    if (m_nseg_bck > 0) {
        rv_bck = propagate_lagrangian(rv_bck, -dt / 2, m_mu, false).first;
    }
    // We now loop through the backward segments and 1) add a dv + 2) propagate for -dt (except on the last segment,
    // where we propagate for -dt/2).
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg_bck; ++i) {
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
        double prop_duration = (i == m_nseg_bck - 1) ? -dt / 2 : -dt;
        rv_bck = propagate_lagrangian(rv_bck, prop_duration, m_mu, false).first;
    }

    return {rv_fwd[0][0] - rv_bck[0][0], rv_fwd[0][1] - rv_bck[0][1], rv_fwd[0][2] - rv_bck[0][2],
            rv_fwd[1][0] - rv_bck[1][0], rv_fwd[1][1] - rv_bck[1][1], rv_fwd[1][2] - rv_bck[1][2],
            mass_fwd - mass_bck};
}

std::vector<double> sims_flanagan::compute_throttle_constraints() const
{
    std::vector<double> retval(m_nseg);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i] = m_throttles[3 * i] * m_throttles[3 * i] + m_throttles[3 * i + 1] * m_throttles[3 * i + 1]
                    + m_throttles[3 * i + 2] * m_throttles[3 * i + 2] - 1.;
    }
    return retval;
}

mat61 _dyn(std::array<std::array<double, 3>, 2> rv, double mu)
{
    mat61 retval;
    auto R3 = std::pow(rv[0][0] * rv[0][0] + rv[0][1] * rv[0][1] + rv[0][2] * rv[0][2], 1.5);
    retval(0, 0) = rv[1][0];
    retval(1, 0) = rv[1][1];
    retval(2, 0) = rv[1][2];
    retval(3, 0) = -mu / R3 * rv[0][0];
    retval(4, 0) = -mu / R3 * rv[0][1];
    retval(5, 0) = -mu / R3 * rv[0][2];
    return retval;
}

// Performs the state updates for nseg sarting from rvs, ms. Computes all gradient information
std::pair<std::array<double, 49>, std::vector<double>> sims_flanagan::gradients_multiple_impulses(
    std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    const std::array<std::array<double, 3>, 2> &rvs, double ms, double c, double a, double dt) const
{
    assert(std::distance(th1, th2) % 3 == 0u);
    auto nseg = static_cast<unsigned>(std::distance(th1, th2) / 3u);

    // Corner case: nseg is zero
    if (nseg == 0) {
        std::array<double, 49> grad_rvm{}; // The mismatch constraints gradient w.r.t. extended state r,v,m
        auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
        xgrad_rvm = xt::eye(7);
        std::vector<double> grad(7, 0.); // The mismatch constraints gradient w.r.t. throttles (0 in this case) and tof
        return std::make_pair(grad_rvm, std::move(grad));
    }
    // Allocate memory.
    std::vector<mat13> u(nseg);
    std::vector<xt::xarray<double>> du(nseg, xt::zeros<double>({3u, nseg * 3u + 2u}));
    std::vector<double> m(nseg + 1, 0.);
    std::vector<xt::xarray<double>> dm(nseg + 1u, xt::zeros<double>({1u, nseg * 3u + 2u}));
    xt::xarray<double> dtof = xt::zeros<double>({1u, nseg * 3u + 2u});
    std::vector<mat13> Dv(nseg);
    std::vector<xt::xarray<double>> dDv(nseg, xt::zeros<double>({3u, nseg * 3u + 2u}));
    std::vector<mat66> M(nseg + 1);  // The STMs
    std::vector<mat66> Mc(nseg + 1); // Mc will contain [Mn@..@M0,Mn@..@M1, Mn]
    std::vector<mat61> f(nseg + 1, xt::zeros<double>({6u, 1u}));
    // Initialize values
    m[0] = ms;
    unsigned i_tmp = 0u;
    for (auto it = th1; it != th2; it += 3) {
        u[i_tmp](0, 0) = *it;
        u[i_tmp](0, 1) = *(it + 1);
        u[i_tmp](0, 2) = *(it + 2);
        du[i_tmp](0, 3 * i_tmp) = 1.;
        du[i_tmp](1, 3 * i_tmp + 1) = 1.;
        du[i_tmp](2, 3 * i_tmp + 2) = 1.;
        i_tmp++;
    }
    dm[0](0, nseg * 3u) = 1.;
    dtof(0, nseg * 3u + 1) = 1.;
    // 1 - We compute the mass schedule and related gradients
    for (decltype(nseg) i = 0; i < nseg; ++i) {
        Dv[i] = c / m[i] * u[i];
        double un = std::sqrt(u[i](0, 0) * u[i](0, 0) + u[i](0, 1) * u[i](0, 1) + u[i](0, 2) * u[i](0, 2));
        double Dvn = c / m[i] * un;
        dDv[i] = c / m[i] * du[i] - c / m[i] / m[i] * xt::linalg::dot(xt::transpose(u[i]), dm[i])
                 + m_max_thrust / m[i] * xt::linalg::dot(xt::transpose(u[i]), dtof) / nseg;
        auto dDvn = c / m[i] / un * xt::linalg::dot(u[i], du[i]) - c / m[i] / m[i] * un * dm[i]
                    + m_max_thrust / m[i] * un * dtof / nseg;
        m[i + 1] = m[i] * std::exp(-Dvn * a);
        dm[i + 1] = -m[i + 1] * a * dDvn + std::exp(-Dvn * a) * dm[i];
    }
    // 2 - We compute the various STMs
    std::array<std::array<double, 3>, 2> rv_it(rvs);
    std::optional<std::array<double, 36>> M_it;
    for (decltype(nseg) i = 0; i < nseg + 1; ++i) {
        auto dur = dt;
        if (i == 0 || i == nseg) {
            dur = dt / 2;
        }

        std::tie(rv_it, M_it) = kep3::propagate_lagrangian(rv_it, dur, m_mu, true);
        // Now we have the STM in M_it, but its a vector, we must operate on an xtensor object instead.
        assert(M_it);
        // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
        M[i] = xt::adapt(*M_it, {6, 6});
        f[i] = _dyn(rv_it, m_mu);
        // And add the impulse if needed
        if (i < nseg) {
            rv_it[1][0] += Dv[i](0, 0);
            rv_it[1][1] += Dv[i](0, 1);
            rv_it[1][2] += Dv[i](0, 2);
        }
    }

    // 3 - We now need to apply the chain rule to assemble the gradients we want (i.e. not w.r.t DV but w.r.t. u etc...)
    mat63 Iv = xt::zeros<double>({6u, 3u}); // This is the gradient of x (rv) w.r.t. v
    Iv(3, 0) = 1.;
    Iv(4, 1) = 1.;
    Iv(5, 2) = 1.;
    Mc[nseg] = M[nseg]; // Mc will contain [Mn@..@M0,Mn@..@M1, Mn]
    for (decltype(nseg) i = 1; i < nseg + 1; ++i) {
        Mc[nseg - i] = _dot(Mc[nseg - i + 1], M[nseg - i]);
    }
    // grad_tof./
    // First the d/dtof term - example: (0.5 * f3 + M3 @ f2 + M3 @ M2 @ f1 + 0.5 * M3 @ M2 @ M1 @ f0) / N

    mat61 grad_tof = 0.5 * f[nseg];
    for (decltype(nseg) i = 0; i + 1 < nseg; ++i) { // i+1 < nseg avoids overflow
        grad_tof += _dot(Mc[i + 2], f[i + 1]);
    }

    grad_tof += 0.5 * _dot(Mc[1], f[0]);
    grad_tof /= nseg;
    // Then we add the d/Dvi * dDvi/dtof - example: M3 @ Iv @ dDv2 + M3 @ M2 @ Iv @ dDv1 + M3 @ M2 @ M1 @ Iv @ dDv0
    for (decltype(nseg) i = 0; i < nseg; ++i) {
        grad_tof += xt::linalg::dot(_dot(Mc[i + 1], Iv),
                                    xt::eval(xt::view(dDv[i], xt::all(), xt::range(nseg * 3 + 1, nseg * 3 + 2))));
    }
    // grad_u
    xt::xarray<double> grad_u = xt::zeros<double>({6u, nseg * 3u});
    for (decltype(nseg) i = 0u; i < nseg; ++i) {
        grad_u += xt::linalg::dot(_dot(Mc[i + 1], Iv), xt::eval(xt::view(dDv[i], xt::all(), xt::range(0, nseg * 3))));
    }
    // grad_ms
    xt::xarray<double> grad_ms = xt::zeros<double>({6u, 1u});
    for (decltype(nseg) i = 0u; i < nseg; ++i) {
        grad_ms += xt::linalg::dot(_dot(Mc[i + 1], Iv),
                                   xt::eval(xt::view(dDv[i], xt::all(), xt::range(nseg * 3, nseg * 3 + 1))));
    }
    // grad_xs
    mat66 grad_xs = Mc[0];

    // Allocate the return values
    std::array<double, 49> grad_rvm{}; // The mismatch constraints gradient w.r.t. extended state r,v,m
    std::vector<double> grad((nseg * 3lu + 1) * 7,
                             0.); // The mismatch constraints gradient w.r.t. throttles and tof
    // Copying in the computed derivatives
    // a) xgrad (the xtensor gradient w.r.t. throttles and tof)
    auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
    auto xgrad = xt::adapt(grad, {7u, nseg * 3 + 1u});
    xt::view(xgrad, xt::range(0u, 6u), xt::range(0u, nseg * 3u)) = grad_u;
    xt::view(xgrad, xt::range(0u, 6u), xt::range(nseg * 3, nseg * 3 + 1)) = grad_tof;
    xt::view(xgrad, xt::range(6u, 7u), xt::all()) = xt::view(dm[nseg], xt::all(), xt::range(0u, nseg * 3 + 1));
    // At this point since the variable order is u,m,tof we have put dmf/dms in rather than dms/dtof. So we fix this.
    xgrad(6u, nseg * 3) = dm[nseg](0, nseg * 3 + 1);
    // b) xgrad_rvm (the xtensor gradient w.r.t. the initial conditions)
    xt::view(xgrad_rvm, xt::range(0, 6), xt::range(0, 6)) = grad_xs;
    xt::view(xgrad_rvm, xt::range(0, 6), xt::range(6, 7)) = grad_ms;
    xgrad_rvm(6, 6) = dm[nseg](0, nseg * 3);
    return std::make_pair(grad_rvm, std::move(grad));
}

std::pair<std::array<double, 49>, std::vector<double>>
sims_flanagan::gradients_fwd(std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<std::array<double, 3>, 2> &rvs, double ms, double c, double a,
                             double dt) const
{
    return gradients_multiple_impulses(th1, th2, rvs, ms, c, a, dt);
}

std::pair<std::array<double, 49>, std::vector<double>>
sims_flanagan::gradients_bck(std::vector<double>::const_iterator th1, std::vector<double>::const_iterator th2,
                             // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                             const std::array<std::array<double, 3>, 2> &rvf_orig, double mf, double c, double a,
                             double dt) const
{
    // 1) we invert the starting velocity.
    auto rvf = rvf_orig;
    rvf[1][0] = -rvf[1][0];
    rvf[1][1] = -rvf[1][1];
    rvf[1][2] = -rvf[1][2];

    // 2) we reverse the throttles ([1,2,3,4,5,6] -> [4,5,6,1,2,3])
    auto size = static_cast<unsigned>(std::distance(th1, th2));
    // Create a new vector to store the reversed values three by three.
    // Here we allocate a vector. Might be not necessary using the C++ range library?
    std::vector<double> reversed_throttles(size);
    // Iterate in reverse order with a step of three
    for (decltype(size) i = 0u, j = size - 1; i < size; i += 3, j -= 3) {
        // Copy three elements at a time in reverse order
        reversed_throttles[j - 2] = *(th1 + i);
        reversed_throttles[j - 1] = *(th1 + i + 1);
        reversed_throttles[j] = *(th1 + i + 2);
    }

    // 3) We reverse the veff, hence veff (a = 1/veff)
    a = -a;

    // 4) We then compute gradients as if this was a forward leg
    auto [grad_rvm, grad]
        = gradients_multiple_impulses(reversed_throttles.begin(), reversed_throttles.end(), rvf, mf, c, a, dt);
    // 5) We have computed dxf/dxs and dxf/dus, but the initial and final velocites (and us) had their sign
    // inverted! We thus need to account for that and change sign once again of the relevant entries.
    // We also must account for changes in the mass equation (now -a)
    auto xgrad_rvm = xt::adapt(grad_rvm, {7u, 7u});
    xt::view(xgrad_rvm, xt::range(3, 6), xt::all()) *= -1; // dvf/dall
    xt::view(xgrad_rvm, xt::all(), xt::range(0, 3)) *= -1; // dmc/drs
    xt::view(xgrad_rvm, xt::all(), xt::range(6, 7)) *= -1; // dmc/dmf

    auto xgrad = xt::adapt(grad, {7u, size + 1u});
    xt::view(xgrad, xt::range(3, 6), xt::all()) *= -1;    // dvf/dall
    xt::view(xgrad, xt::all(), xt::range(0, size)) *= -1; // dmc/dus

    // 6) Note that the throttles in xgrad are ordered in reverse. Before returning we must restore the forward order
    xt::view(xgrad, xt::all(), xt::range(0, size)) = xt::flip(xt::view(xgrad, xt::all(), xt::range(0, size)), 1);
    for (decltype(size) i = 0u; i < size / 3; ++i) {
        xt::view(xgrad, xt::all(), xt::range(3 * i, 3 * i + 3))
            = xt::flip(xt::view(xgrad, xt::all(), xt::range(3 * i, 3 * i + 3)), 1);
    }
    // And finally return.
    return std::make_pair(grad_rvm, std::move(grad));
}

// Computes the gradient of the mismatch constraints w.r.t. xs, xf and [throttles, tof]
std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>> sims_flanagan::compute_mc_grad() const
{
    // Preliminaries
    const auto dt = m_tof / static_cast<double>(m_nseg); // dt
    const auto c = m_max_thrust * dt;                    // T*tof/nseg
    const auto a = 1. / m_veff;                // 1/veff

    // We compute for the forward half-leg: dxf/dxs and dxf/dxu (the gradients w.r.t. initial state ant throttles )
    auto [grad_rvm, grad_fwd]
        = gradients_fwd(m_throttles.begin(), m_throttles.begin() + static_cast<unsigned>(3 * m_nseg_fwd), get_rvs(),
                        get_ms(), c, a, dt);
    // We compute for the backward half-leg: dxf/dxs and dxf/dxu (the gradients w.r.t. final state and throttles )
    auto [grad_rvm_bck, grad_bck] = gradients_bck(m_throttles.begin() + static_cast<unsigned>(3 * m_nseg_fwd),
                                                  m_throttles.end(), get_rvf(), get_mf(), c, a, dt);

    // We assemble the final results
    std::vector<double> grad_final(static_cast<size_t>(7) * (m_nseg * 3u + 1u), 0.);
    auto xgrad_final = xt::adapt(grad_final, {7u, static_cast<unsigned>(m_nseg) * 3u + 1u});
    auto xgrad_fwd = xt::adapt(grad_fwd, {7u, static_cast<unsigned>(m_nseg_fwd) * 3u + 1u});
    auto xgrad_bck = xt::adapt(grad_bck, {7u, static_cast<unsigned>(m_nseg - m_nseg_fwd) * 3u + 1u});

    // Copy the gradient w.r.t. the forward throttles as is
    xt::view(xgrad_final, xt::all(), xt::range(0, m_nseg_fwd * 3))
        = xt::view(xgrad_fwd, xt::all(), xt::range(0, m_nseg_fwd * 3));

    // Copy the gradient w.r.t. the backward throttles as is
    xt::view(xgrad_final, xt::all(), xt::range(m_nseg_fwd * 3, m_nseg * 3))
        = xt::view(xgrad_bck, xt::all(), xt::range(0, (m_nseg - m_nseg_fwd) * 3));

    // Copy the gradient w.r.t. tof as fwd-bck
    xt::view(xgrad_final, xt::all(), xt::range(m_nseg * 3, m_nseg * 3 + 1))
        = xt::view(xgrad_fwd, xt::all(), xt::range(m_nseg_fwd * 3, m_nseg_fwd * 3 + 1)) / m_nseg * m_nseg_fwd
          - xt::view(xgrad_bck, xt::all(), xt::range((m_nseg - m_nseg_fwd) * 3, (m_nseg - m_nseg_fwd) * 3 + 1)) / m_nseg
                * (m_nseg - m_nseg_fwd);
    return {grad_rvm, grad_rvm_bck, std::move(grad_final)};
}

std::vector<double> sims_flanagan::compute_tc_grad() const
{
    std::vector<double> retval(static_cast<size_t>(m_nseg) * m_nseg * 3u, 0);
    for (decltype(m_throttles.size()) i = 0u; i < m_nseg; ++i) {
        retval[i * m_nseg * 3 + 3 * i] = 2 * m_throttles[3 * i];
        retval[i * m_nseg * 3 + 3 * i + 1] = 2 * m_throttles[3 * i + 1];
        retval[i * m_nseg * 3 + 3 * i + 2] = 2 * m_throttles[3 * i + 2];
    }
    return retval;
}

std::ostream &operator<<(std::ostream &s, const sims_flanagan &sf)
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