// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_HELPERS_H
#define kep3_TEST_HELPERS_H

#include <algorithm>
#include <array>
#include <cmath>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xadapt.hpp>

namespace kep3_tests
{
using xt::linalg::cross;
using xt::linalg::dot;
// This is a float test which, while controversial, will test for abs
// differences in small numbers, relative otherwise.
inline double floating_point_error(double a, double b)
{
    return std::abs(a - b) / std::max(1., std::max(std::abs(a), std::abs(b)));
}

// Self explanatory
template <typename T>
double L_infinity_norm(T a, T b)
{
    if (a.size() != b.size()) {
        throw std::domain_error("Subtracting two vectors having unequal size.");
    }
    double retval = 0.;

    for (decltype(a.size()) i=0u; i < a.size() ; ++i) {
        std::abs(a[i]-b[i]) > retval ? retval = std::abs(a[i]-b[i]) : retval;
    }
    return retval;
}

// Takes the floating_point_error to compute the infinity norm ... (thus not an infinity norm)
template <typename T>
double L_infinity_norm_rel(T a, T b)
{
    double retval = 0.;
    for (decltype(a.size()) i=0u; i < a.size() ; ++i) {
        std::abs(a[i]-b[i]) > retval ? retval = floating_point_error(a[i], b[i]) : retval;
    }
    return retval;
}

// This tests how close two 3D vectors are in the euclidean metric. err = |(r2-r1)|
inline double floating_point_error_vector(const std::array<double, 3> &r1, const std::array<double, 3> &r2)
{
    double const R1 = std::sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    std::array<double, 3> r12 = {{r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]}};
    double const R12 = std::sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
    return R12 / std::max(1., R1);
}

// see Battin: "An Introduction to the Mathematics and Methods of
// Astrodynamics, Revised Edition", Introduction.
//
// On Keplerian dynamics the following must hold.
//
// (v1 x r1).(v1 x (r2 - r1)) + mu r2 . (r2/|r2| - r1/{|r1|})
inline double delta_guidance_error(const std::array<double, 3> &r1, const std::array<double, 3> &r2,
                                   const std::array<double, 3> &v1, double mu)
{
    const auto r1_x = xt::adapt(r1);
    const auto r2_x = xt::adapt(r2);
    const auto v1_x = xt::adapt(v1);
    auto F = dot(cross(v1_x, r1_x), cross(v1_x, (r2_x - r1_x)))
             + mu * dot(r2_x, (r2_x / xt::linalg::norm(r2_x) - r1_x / xt::linalg::norm(r1_x)));
    return F(0);
}
} // namespace kep3_tests

#endif // kep3_TEST_HELPERS_H
