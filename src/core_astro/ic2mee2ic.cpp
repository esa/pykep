// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>

#include <fmt/core.h>

#include <fmt/ranges.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>

using xt::linalg::cross;
using xt::linalg::dot;

namespace kep3
{

// Implementation following:
// Cefola: Equinoctial orbit elements - Application to artificial satellite
// orbitsCefola, P., 1972, September. Equinoctial orbit elements-Application to
// artificial satellite orbits. In Astrodynamics Conference (p. 937).

std::array<double, 6> ic2mee(const std::array<std::array<double, 3>, 2> &pos_vel, double mu, bool retrogade)
{
    {
        // Switch between the element types.
        int I = 0;
        if (retrogade) {
            I = -1;
        } else {
            I = 1;
        }
        // 0 - We prepare a few xtensor constants.
        auto r0 = xt::adapt(pos_vel[0]);
        auto v0 = xt::adapt(pos_vel[1]);
        // The equinoctial reference frame
        xt::xtensor_fixed<double, xt::xshape<3>> fv = {0.0, 0.0, 0.0};
        xt::xtensor_fixed<double, xt::xshape<3>> gv = {0.0, 0.0, 0.0};

        // angular momentum
        auto ang = cross(r0, v0);

        // 0 - We compute the semi-major axis
        double R0 = xt::linalg::norm(r0);
        double V0 = xt::linalg::norm(v0);
        double sma = 1. / (2. / R0 - V0 * V0 / mu);

        // 1 - We compute the equinoctial frame
        auto w = cross(r0, v0);
        w = w / xt::linalg::norm(w);

        double k = w(0) / (1. + I * w(2));
        double h = -w(1) / (1. + I * w(2));
        double den = k * k + h * h + 1;
        fv(0) = (1. - k * k + h * h) / den;
        fv(1) = (2. * k * h) / den;
        fv(2) = (-2. * I * k) / den;

        gv(0) = (2. * I * k * h) / den;
        gv(1) = (1. + k * k - h * h) * I / den;
        gv(2) = (2. * h) / den;

        // 2 - We compute evett: the eccentricity vector
        auto evett = cross(v0, ang) / mu - r0 / R0; // e = (v x h)/mu - r0/R0;

        double g = dot(evett, gv)(0);
        double f = dot(evett, fv)(0);
        double ecc = xt::linalg::norm(evett);

        // 3 - We compute the true longitude L
        // This solution is certainly not the most elegant, but it works and will
        // never be singular.

        double det1 = (gv(1) * fv(0) - fv(1) * gv(0)); // xy
        double det2 = (gv(2) * fv(0) - fv(2) * gv(0)); // xz
        double det3 = (gv(2) * fv(1) - fv(2) * gv(1)); // yz

        std::array<double, 3> dets = {std::abs(det1), std::abs(det2), std::abs(det3)};
        std::size_t idx
            = static_cast<std::size_t>(std::distance(dets.begin(), std::max_element(dets.begin(), dets.end())));

        double X = 0., Y = 0.;
        switch (idx) {
            case 0:
                X = (gv(1) * r0(0) - gv(0) * r0(1)) / det1;
                Y = (-fv(1) * r0(0) + fv(0) * r0(1)) / det1;
                break;
            case 1:
                X = (gv(2) * r0(0) - gv(0) * r0(2)) / det2;
                Y = (-fv(2) * r0(0) + fv(0) * r0(2)) / det2;
                break;
            case 2:
                X = (gv(2) * r0(1) - gv(1) * r0(2)) / det3;
                Y = (-fv(2) * r0(1) + fv(1) * r0(2)) / det3;
                break;
        }

        double L = std::atan2(Y, X);

        // 5 - We assign the results
        return {sma * (1. - ecc * ecc), f, g, h, k, L};
    }
}

std::array<std::array<double, 3>, 2> mee2ic(const std::array<double, 6> &eq, double mu, bool retrogade)
{
    std::array<std::array<double, 3>, 2> retval{};
    int I = 0;
    if (retrogade) {
        I = -1;
    } else {
        I = 1;
    }

    // p = a (1-e^2) will be negative for eccentricities > 1, we here need a
    // positive number for the following computations to make sense
    double par = std::abs(eq[0]);
    double f = eq[1];
    double g = eq[2];
    double h = eq[3];
    double k = eq[4];
    double L = eq[5];

    // We compute the equinoctial reference frame
    double den = k * k + h * h + 1;
    double fx = (1 - k * k + h * h) / den;
    double fy = (2 * k * h) / den;
    double fz = (-2 * I * k) / den;

    double gx = (2 * I * k * h) / den;
    double gy = (1 + k * k - h * h) * I / den;
    double gz = (2 * h) / den;

    // Auxiliary
    double radius = par / (1 + g * std::sin(L) + f * std::cos(L));
    // In the equinoctial reference frame
    double X = radius * std::cos(L);
    double Y = radius * std::sin(L);
    double VX = -std::sqrt(mu / par) * (g + std::sin(L));
    double VY = std::sqrt(mu / par) * (f + std::cos(L));

    // Results
    retval[0][0] = X * fx + Y * gx;
    retval[0][1] = X * fy + Y * gy;
    retval[0][2] = X * fz + Y * gz;

    retval[1][0] = VX * fx + VY * gx;
    retval[1][1] = VX * fy + VY * gy;
    retval[1][2] = VX * fz + VY * gz;

    return retval;
}
} // namespace kep3