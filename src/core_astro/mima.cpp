// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cmath>

#include <boost/math/tools/roots.hpp>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include <kep3/core_astro/mima.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/linalg.hpp>

#include "kep3/core_astro/constants.hpp"

namespace kep3
{

std::pair<double, double> mima(const std::array<double, 3> &dv1, const std::array<double, 3> &dv2, double tof,
                               double Tmax, // NOLINT
                               double veff) // NOLINT
{
    std::array<double, 3> dv{dv1[0] + dv2[0], dv1[1] + dv2[1], dv1[2] + dv2[2]};
    std::array<double, 3> dv_diff{-dv1[0] + dv2[0], -dv1[1] + dv2[1], -dv1[2] + dv2[2]};
    double ab = dv[0] * dv_diff[0] + dv[1] * dv_diff[1] + dv[2] * dv_diff[2];                // ab = dv@dv_diff
    double aa = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];                               // aa = dv@dv
    double bb = dv_diff[0] * dv_diff[0] + dv_diff[1] * dv_diff[1] + dv_diff[2] * dv_diff[2]; // bb = dv_diff@dv_diff
    double ad = std::sqrt(aa + 2 * bb + 2 * std::sqrt((ab * ab + bb * bb))) / tof;
    double mima = 2 * Tmax / ad / (1 + std::exp(-ad * tof / veff));
    return {mima, ad};
}

std::pair<double, double> mima_from_hop(const kep3::planet &pl_s, const kep3::planet &pl_f, const kep3::epoch &when_s,
                                        const kep3::epoch &when_f,
                                        double Tmax, // NOLINT
                                        double veff) // NOLINT
{
    double tof = (when_f.mjd2000() - when_s.mjd2000()) * kep3::DAY2SEC;
    double mu = pl_s.get_mu_central_body();
    const auto &[r_s, v_s] = pl_s.eph(when_s);
    const auto &[r_f, v_f] = pl_f.eph(when_f);
    auto l = kep3::lambert_problem(r_s, r_f, tof, mu, false, 0u);
    std::array<double, 3> dv1 = {l.get_v0()[0][0] - v_s[0], l.get_v0()[0][1] - v_s[1], l.get_v0()[0][2] - v_s[2]};
    std::array<double, 3> dv2 = {-l.get_v1()[0][0] + v_f[0], -l.get_v1()[0][1] + v_f[1], -l.get_v1()[0][2] + v_f[2]};
    return mima(dv1, dv2, tof, Tmax, veff);
}

std::pair<double, double> _mima_compute_transfer(double x, const std::array<std::array<double, 3>, 2> &posvel1,
                                                 double tof, const std::array<double, 3> &dv1_flat,
                                                 const std::array<double, 3> &dv2_flat, double mu)
{
    // Some shortcuts for xtensor functions
    using kep3::linalg::_dot;
    using kep3::linalg::mat31;
    using kep3::linalg::mat33;
    using kep3::linalg::mat36;
    using kep3::linalg::mat61;
    using kep3::linalg::mat66;
    using xt::linalg::inv;

    // Start of the algorithm (see the paper for details)
    double tau = (x / std::sqrt(x * x + 1.) + 1.) / 2.;
    double t1 = tof * tau;
    double t2 = tof * (1. - tau);
    auto res12 = propagate_lagrangian(posvel1, t1 / 2., mu, true);
    auto res22 = propagate_lagrangian(posvel1, tof - t2 / 2., mu, true);
    auto res11 = propagate_lagrangian(posvel1, t1, mu, true);
    auto res21 = propagate_lagrangian(posvel1, tof - t2, mu, true);
    auto res3 = propagate_lagrangian(posvel1, tof, mu, true);

    // Adapting all kep3 flattened vectors and arrays into xtensor objects
    mat31 dv1 = xt::adapt(dv1_flat);
    mat31 dv2 = xt::adapt(dv2_flat);
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    mat66 M12 = xt::adapt(res12.second.value());
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    mat66 M22 = xt::adapt(res22.second.value());
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    mat66 M11 = xt::adapt(res11.second.value());
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    mat66 M21 = xt::adapt(res21.second.value());
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    mat66 M3 = xt::adapt(res3.second.value());

    // Simpson's rule
    mat66 &stm0 = M3; // alias for convenience
    mat66 tmp = inv(M11) + 4 * inv(M12);
    mat66 stm1 = (_dot(M3, tmp) + M3) / 6.;
    mat66 tmp2 = inv(M21) + 4 * inv(M22);
    mat66 stm2 = (_dot(M3, tmp2) + xt::eye<double>(6)) / 6.;
    // b = np.hstack((stm0[0:3, 3:6]@dv1, dv2+stm0[3:6, 3:6]@dv1))
    mat33 b1_tmp = xt::view(stm0, xt::range(0, 3), xt::range(3, 6));
    mat31 b1 = _dot(b1_tmp, dv1);
    mat33 b2_tmp = xt::view(stm0, xt::range(3, 6), xt::range(3, 6));
    mat31 b2 = dv2 + _dot(b2_tmp, dv1);
    mat61 b = xt::concatenate(xt::xtuple(b1, b2));
    // M = np.vstack((np.hstack((stm1[0:3, 3:6], stm2[0:3, 3:6])), np.hstack(
    //    (stm1[3:6, 3:6], stm2[3:6, 3:6]))))
    mat36 M1 = xt::concatenate(
        xt::xtuple(xt::view(stm1, xt::range(0, 3), xt::range(3, 6)), xt::view(stm2, xt::range(0, 3), xt::range(3, 6))),
        1);
    mat36 M2 = xt::concatenate(
        xt::xtuple(xt::view(stm1, xt::range(3, 6), xt::range(3, 6)), xt::view(stm2, xt::range(3, 6), xt::range(3, 6))),
        1);
    mat66 M = xt::concatenate(xt::xtuple(M1, M2));
    // dvs = np.linalg.inv(M)@b
    // a1 = dvs[0:3]/tau/tof
    // a2 = dvs[3:6]/(1-tau)/tof
    // err = tau**2*(1-tau)**2*(a2[0]*a2[0]+a2[1]*a2[1]+a2[2]*a2[2]-a1[0]*a1[0]-a1[1]*a1[1]-a1[2]*a1[2])
    mat66 invM = inv(M);
    mat61 dvs = _dot(invM, b);
    mat31 a1 = xt::view(dvs, xt::range(0, 3)) / tau / tof;
    mat31 a2 = xt::view(dvs, xt::range(3, 6)) / (1. - tau) / tof;

    double err = tau * tau * (1 - tau) * (1 - tau)
                 * (a2(0, 0) * a2(0, 0) + a2(1, 0) * a2(1, 0) + a2(2, 0) * a2(2, 0) - a1(0, 0) * a1(0, 0)
                    - a1(1, 0) * a1(1, 0) - a1(2, 0) * a1(2, 0));
    return {err, std::sqrt(a1(0, 0) * a1(0, 0) + a1(1, 0) * a1(1, 0) + a1(2, 0) * a1(2, 0))};
}

kep3_DLL_PUBLIC std::pair<double, double> mima2(const std::array<std::array<double, 3>, 2> &posvel1,
                                                const std::array<double, 3> &dv1, const std::array<double, 3> &dv2,
                                                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                                                double tof, double Tmax, double veff, double mu)
{
    const boost::uintmax_t maxit = 100u;
    boost::uintmax_t it = maxit;
    unsigned digits = 10u;                                 // No need to compute this approximation precisely
    boost::math::tools::eps_tolerance<double> tol(digits); // Set the tolerance.
    double guess = _mima_compute_transfer(0., posvel1, tof, dv1, dv2, mu).first > 0 ? -0.5 : 0.5;
    auto r = boost::math::tools::bracket_and_solve_root(
        [&](double x) { return _mima_compute_transfer(x, posvel1, tof, dv1, dv2, mu).first; }, guess, 2., true, tol,
        it);

    if (it >= maxit) { //
        throw std::domain_error("Maximum number of iterations exceeded when computing mima2");
    }
    auto root = r.first + (r.second - r.first) / 2; // Midway between brackets.
    auto acc = _mima_compute_transfer(root, posvel1, tof, dv1, dv2, mu).second;
    auto mima2 = 2. * Tmax / acc / (1. + std::exp(-acc * tof / veff));
    return {mima2, acc};
}

std::pair<double, double> mima2_from_hop(const kep3::planet &pl_s, const kep3::planet &pl_f, const kep3::epoch &when_s,
                                         const kep3::epoch &when_f,
                                         double Tmax, // NOLINT
                                         double veff) // NOLINT
{
    double tof = (when_f.mjd2000() - when_s.mjd2000()) * kep3::DAY2SEC;
    double mu = pl_s.get_mu_central_body();
    const auto &[r_s, v_s] = pl_s.eph(when_s);
    const auto &[r_f, v_f] = pl_f.eph(when_f);
    auto l = kep3::lambert_problem(r_s, r_f, tof, mu, false, 0u);
    std::array<double, 3> dv1 = {l.get_v0()[0][0] - v_s[0], l.get_v0()[0][1] - v_s[1], l.get_v0()[0][2] - v_s[2]};
    std::array<double, 3> dv2 = {-l.get_v1()[0][0] + v_f[0], -l.get_v1()[0][1] + v_f[1], -l.get_v1()[0][2] + v_f[2]};
    return mima2({r_s, l.get_v0()[0]}, dv1, dv2, tof, Tmax, veff, mu);
}

} // namespace kep3