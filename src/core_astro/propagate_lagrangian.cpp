// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

#include <boost/math/tools/roots.hpp>
#include <fmt/core.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/kepler_equations.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/core_astro/special_functions.hpp>

namespace kep3
{

/// Lagrangian propagation
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion. Lagrange coefficients are used as basic
 * numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 */
std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_lagrangian(const std::array<std::array<double, 3>, 2> &pos_vel0, const double tof, const double mu, bool stm)
{
    const auto &[r0, v0] = pos_vel0;
    auto pos_velf = pos_vel0;
    auto &[rf, vf] = pos_velf;
    double R0 = std::sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    double Rf = 0.;
    double const V02 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
    double DX = 0.;
    double const energy = (V02 / 2 - mu / R0);
    double a = -mu / 2.0 / energy; // will be negative for hyperbolae
    double sqrta = 0.;
    double F = 0., G = 0., Ft = 0., Gt = 0.;
    double s0 = 0., c0 = 0.;
    double sigma0 = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / std::sqrt(mu);

    if (a > 0) { // Solve Kepler's equation in DE, elliptical case
        sqrta = std::sqrt(a);
        double DM = std::sqrt(mu / std::pow(a, 3)) * tof;
        double const sinDM = std::sin(DM);
        double const cosDM = std::cos(DM);
        // Here we use the atan2 to recover the mean anomaly difference in the
        // [0,2pi] range. This makes sure that for high value of M no catastrophic
        // cancellation occurs, as would be the case using std::fmod(DM, 2pi)
        double DM_cropped = std::atan2(sinDM, cosDM);
        if (DM_cropped < 0) {
            DM_cropped += 2 * kep3::pi;
        }
        s0 = sigma0 / sqrta;
        c0 = (1 - R0 / a);
        // This initial guess was developed applying Lagrange expansion theorem to
        // the Kepler's equation in DE. We stopped at 3rd order.
        double const IG = DM_cropped + c0 * sinDM - s0 * (1 - cosDM)
                          + (c0 * cosDM - s0 * sinDM) * (c0 * sinDM + s0 * cosDM - s0)
                          + 0.5 * (c0 * sinDM + s0 * cosDM - s0)
                                * (2 * std::pow(c0 * cosDM - s0 * sinDM, 2)
                                   - (c0 * sinDM + s0 * cosDM - s0) * (c0 * sinDM + s0 * cosDM));

        // Solve Kepler Equation for ellipses in DE (eccentric anomaly difference)
        const int digits = std::numeric_limits<double>::digits;
        std::uintmax_t max_iter = 100u;
        // NOTE: Halley iterates may result into instabilities (specially with a
        // poor IG)

        double DE = boost::math::tools::newton_raphson_iterate(
            [DM_cropped, sigma0, sqrta, a, R0](double DE) {
                return std::make_tuple(kepDE(DE, DM_cropped, sigma0, sqrta, a, R0), d_kepDE(DE, sigma0, sqrta, a, R0));
            },
            IG, IG - pi, IG + pi, digits, max_iter);
        // LCOV_EXCL_START
        if (max_iter == 100u) {
            throw std::domain_error(fmt::format("Maximum number of iterations exceeded when solving Kepler's "
                                                "equation for the eccentric anomaly in propagate_lagrangian.\n"
                                                "DM={}\nsigma0={}\nsqrta={}\na={}\nR={}\nDE={}",
                                                DM, sigma0, sqrta, a, R0, DE));
        }
        // LCOV_EXCL_STOP
        Rf = a + (R0 - a) * std::cos(DE) + sigma0 * sqrta * std::sin(DE);

        // Lagrange coefficients
        F = 1 - a / R0 * (1 - std::cos(DE));
        G = a * sigma0 / std::sqrt(mu) * (1 - std::cos(DE)) + R0 * std::sqrt(a / mu) * std::sin(DE);
        Ft = -std::sqrt(mu * a) / (Rf * R0) * std::sin(DE);
        Gt = 1 - a / Rf * (1 - std::cos(DE));
        DX = DE;
    } else { // Solve Kepler's equation in DH, hyperbolic case
        sqrta = std::sqrt(-a);
        double DN = std::sqrt(-mu / a / a / a) * tof;
        double IG = 0.;
        s0 = sigma0 / sqrta;
        c0 = (1 - R0 / a);
        tof > 0. ? IG = 1. : IG = -1.; // TODO(darioizzo): find a better initial guess.
                                       // I tried with 0 and DN (both have numercial
                                       // problems and result in exceptions)

        // Solve Kepler Equation for ellipses in DH (hyperbolic anomaly difference)
        const int digits = std::numeric_limits<double>::digits;
        std::uintmax_t max_iter = 100u;
        // NOTE: Halley iterates may result into instabilities (specially with a
        // poor IG)
        double DH = boost::math::tools::newton_raphson_iterate(
            [DN, sigma0, sqrta, a, R0](double DH) {
                return std::make_tuple(kepDH(DH, DN, sigma0, sqrta, a, R0), d_kepDH(DH, sigma0, sqrta, a, R0));
            },
            IG, IG - 50, IG + 50, digits,
            max_iter); // TODO (dario): study this hyperbolic equation in more
                       // details as to provide decent and well proved bounds
        // LCOV_EXCL_START
        if (max_iter == 100u) {
            throw std::domain_error(fmt::format("Maximum number of iterations exceeded when solving Kepler's "
                                                "equation for the hyperbolic anomaly in propagate_lagrangian.\n"
                                                "DN={}\nsigma0={}\nsqrta={}\na={}\nR={}\nDH={}",
                                                DN, sigma0, sqrta, a, R0, DH));
        }
        // LCOV_EXCL_STOP

        // Note: the following equation, according to Battin's (4.63, pag 170), is different. I suspect a typo in Battin
        // book deriving from the confusion on the convention for the semi-majoraxis being negative. The following
        // expression instead works in this context.
        Rf = a + (R0 - a) * std::cosh(DH) + sigma0 * sqrta * std::sinh(DH);

        // Lagrange coefficients
        F = 1. - a / R0 * (1. - std::cosh(DH));
        G = a * sigma0 / std::sqrt(mu) * (1. - std::cosh(DH)) + R0 * std::sqrt(-a / mu) * std::sinh(DH);
        Ft = -std::sqrt(-mu * a) / (Rf * R0) * std::sinh(DH);
        Gt = 1. - a / Rf * (1. - std::cosh(DH));
        DX = DH;
    }

    for (auto i = 0u; i < 3; i++) {
        rf[i] = F * r0[i] + G * v0[i];
        vf[i] = Ft * r0[i] + Gt * v0[i];
    }
    if (stm) {
        auto retval_stm = kep3::stm_lagrangian(pos_vel0, tof, mu, R0, Rf, energy, sigma0, a, s0, c0, DX, F, G, Ft, Gt);
        return {pos_velf, retval_stm};
    } else {
        return {pos_velf, std::nullopt};
    }
}

std::vector<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>>
propagate_lagrangian_grid(const std::array<std::array<double, 3>, 2> &pos_vel, const std::vector<double> &time_grid,
                       double mu, bool stm)
{
    auto n = time_grid.size();
    std::vector<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>> retval(n);
    for (decltype(n) i = 0u; i < n; ++i) {
        retval[i] = kep3::propagate_lagrangian(pos_vel, time_grid[i] - time_grid[0], mu, stm);
    }
    return retval;
}

/// Universial Variables version
/**
 * This function has the same prototype as kep3::propagate_lgrangian, but
 * internally makes use of universal variables formulation for the Lagrange
 * Coefficients. Its slower so not the main choice in kep3.
 */
std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_lagrangian_u(const std::array<std::array<double, 3>, 2> &pos_vel0, const double dt, const double mu, // NOLINT
                       bool)
{
    // If time is negative we need to invert time and velocities. Unlike the other
    // formulation of the propagate lagrangian we cannot rely on negative times to
    // automatically mean back-propagation
    double dt_copy = dt;
    const auto &[r0, v0] = pos_vel0;
    std::array<std::array<double, 3>, 2> pos_velf = pos_vel0;
    auto &[rf, vf] = pos_velf;

    // posvelf is at this point storing v0 (will store vf at the end). We cannot use v0 here as its const.
    if (dt < 0) {
        dt_copy = -dt;
        vf[0] = -vf[0];
        vf[1] = -vf[1];
        vf[2] = -vf[2];
    }

    double F = 0., G = 0., Ft = 0., Gt = 0.;
    double const R0 = std::sqrt(rf[0] * rf[0] + rf[1] * rf[1] + rf[2] * rf[2]);
    double const V02 = vf[0] * vf[0] + vf[1] * vf[1] + vf[2] * vf[2];
    // the reciprocal of the semi-major axis
    double const alpha = 2 / R0 - V02 / mu;
    // initial radial velocity (see we use the inverted vf = +-v0 according to dt)
    double const VR0 = (rf[0] * vf[0] + rf[1] * vf[1] + rf[2] * vf[2]) / R0;

    // solve kepler's equation in the universal anomaly DS
    double IG = 0;
    alpha > 0. ? IG = std::sqrt(mu) * dt_copy * std::abs(alpha)
               : IG = 3.; // TODO(darioizzo): initial guess for the universal
                          // anomaly. For hyperbolas is 3 .... can be better?

    // Solve Kepler Equation in DS (univrsal anomaly difference)
    const int digits = std::numeric_limits<double>::digits;
    std::uintmax_t max_iter = 100u;
    // NOTE: Halley iterates may result into instabilities (specially with a poor
    // IG)
    double const DS = boost::math::tools::newton_raphson_iterate(
        [dt_copy, R0, VR0, alpha, mu](double DS) {
            return std::make_tuple(kepDS(DS, dt_copy, R0, VR0, alpha, mu), d_kepDS(DS, R0, VR0, alpha, mu));
        },
        IG, IG - 2 * pi, IG + 2 * pi, digits,
        max_iter); // limiting the IG error within
                   // only pi will not work.
    if (max_iter == 100u) {
        throw std::domain_error("Maximum number of iterations exceeded when solving Kepler's "
                                "equation for the universal anomaly in propagate_lagrangian_u.");
    }
    // evaluate the lagrangian coefficients F and G
    double const S = stumpff_s(alpha * DS * DS);
    double const C = stumpff_c(alpha * DS * DS);
    //
    double const z = alpha * DS * DS;
    F = 1 - DS * DS / R0 * C;
    G = dt_copy - 1 / std::sqrt(mu) * DS * DS * DS * S;

    // compute the final position (vf = +-v0 according to dt)
    rf[0] = F * r0[0] + G * vf[0];
    rf[1] = F * r0[1] + G * vf[1];
    rf[2] = F * r0[2] + G * vf[2];
    double const RF = std::sqrt(rf[0] * rf[0] + rf[1] * rf[1] + rf[2] * rf[2]);

    // compute the lagrangian coefficients Ft, Gt
    Ft = std::sqrt(mu) / RF / R0 * (z * S - 1) * DS;
    Gt = 1 - DS * DS / RF * C;

    // compute the final velocity (vf = +-v0 according to dt)
    vf[0] = Ft * r0[0] + Gt * vf[0];
    vf[1] = Ft * r0[1] + Gt * vf[1];
    vf[2] = Ft * r0[2] + Gt * vf[2];

    if (dt < 0) {
        vf[0] = -vf[0];
        vf[1] = -vf[1];
        vf[2] = -vf[2];
    }
    // since the stm is here not implemented yet, we return always nullopt.
    return {pos_velf, std::nullopt};
}

/// Keplerian (not using the lagrangian coefficients) propagation
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion. Simple conversions are used to compute
 * M0 then Mt, etc.. It only here for study purposes as its x10 slower (strange
 * such a high factor ..investigate?)
 */
std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_keplerian(const std::array<std::array<double, 3>, 2> &pos_vel0, const double dt, const double mu, // NOLINT
                    bool)
{
    // 1 - Compute the orbital parameters at t0
    auto par = kep3::ic2par(pos_vel0, mu);
    if (par[0] > 0) {
        // 2e - Compute the mean anomalies
        double const n = std::sqrt(mu / par[0] / par[0] / par[0]);
        double const M0 = kep3::f2m(par[5], par[1]);
        double const Mf = M0 + n * dt;
        // 3e - Update elements (here Kepler's equation gets solved)
        par[5] = kep3::m2f(Mf, par[1]);
    } else {
        // 2h - Compute the mean hyperbolic anomalies
        double const n = std::sqrt(-mu / par[0] / par[0] / par[0]);
        double const N0 = kep3::f2n(par[5], par[1]);
        double const Nf = N0 + n * dt;
        // 3h - Update elements (here Kepler's equation gets solved in its
        // hyperbolic version)
        par[5] = kep3::n2f(Nf, par[1]);
    }
    // Update posvel
    auto pos_velf = kep3::par2ic(par, mu);
    return {pos_velf, std::nullopt};
}

} // namespace kep3