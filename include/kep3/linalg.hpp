// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_DETAIL_XTENSOR_HELPERS_HPP
#define kep3_DETAIL_XTENSOR_HELPERS_HPP

#include <xtensor/containers/xarray.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>

namespace kep3::linalg
{
using mat31 = xt::xtensor_fixed<double, xt::xshape<3, 1>>;
using mat13 = xt::xtensor_fixed<double, xt::xshape<1, 3>>;
using mat33 = xt::xtensor_fixed<double, xt::xshape<3, 3>>;
using mat36 = xt::xtensor_fixed<double, xt::xshape<3, 6>>;
using mat66 = xt::xtensor_fixed<double, xt::xshape<6, 6>>;
using mat32 = xt::xtensor_fixed<double, xt::xshape<3, 2>>;
using mat62 = xt::xtensor_fixed<double, xt::xshape<6, 2>>;
using mat61 = xt::xtensor_fixed<double, xt::xshape<6, 1>>;
using mat16 = xt::xtensor_fixed<double, xt::xshape<1, 6>>;
using mat63 = xt::xtensor_fixed<double, xt::xshape<6, 3>>;
using mat77 = xt::xtensor_fixed<double, xt::xshape<7, 7>>;
using mat74 = xt::xtensor_fixed<double, xt::xshape<7, 4>>;
using mat71 = xt::xtensor_fixed<double, xt::xshape<7, 1>>;
using mat17 = xt::xtensor_fixed<double, xt::xshape<1, 7>>;

// -----------------------------------------------------------------------------------------
// Linear algebra helpers to speed up xtensor when small, fixed size matrices and vectors are involved

// Matrix multiplication
// Overload for xt::xtensor_fixed types: compile-time shape enforcement, no need to specify template params explicitly
template <std::size_t m, std::size_t k, std::size_t n>
xt::xtensor_fixed<double, xt::xshape<m, n>>
_dot(const xt::xtensor_fixed<double, xt::xshape<m, k>>& a,
     const xt::xtensor_fixed<double, xt::xshape<k, n>>& b)
{
    auto A_eval = xt::eval(a);
    auto B_eval = xt::eval(b);
    xt::xtensor_fixed<double, xt::xshape<m, n>> C{};
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = 0;
            for (size_t l = 0; l < k; ++l) {
                C(i, j) += A_eval(i, l) * B_eval(l, j);
            }
        }
    }
    return C;
}

// General template for other types (expressions, dynamic shapes, etc.)
template <std::size_t m, std::size_t k, std::size_t n, class A, class B>
xt::xtensor_fixed<double, xt::xshape<m, n>> _dot(const A& a, const B& b) {
    auto A_eval = xt::eval(a);
    auto B_eval = xt::eval(b);
    xt::xtensor_fixed<double, xt::xshape<m, n>> C{};
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C(i, j) = 0;
            for (size_t l = 0; l < k; ++l) {
                C(i, j) += A_eval(i, l) * B_eval(l, j);
            }
        }
    }
    return C;
}

// Small linear algebra operations
mat33 _skew(const mat31 &v);
mat31 _cross(const mat31 &v1, const mat31 &v2);
// ---------------------------------------------------------------------------------------

// Linear algebra helpers for std::array<double, 3> types.
inline void _normalize(std::array<double, 3> &v1)
{
    const double norm = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    v1[0] /= norm;
    v1[1] /= norm;
    v1[2] /= norm;
}

inline std::array<double, 3> operator+(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline std::array<double, 3> operator-(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline std::array<double, 3> operator*(const std::array<double, 3> &a, double s)
{
    return {a[0] * s, a[1] * s, a[2] * s};
}

inline std::array<double, 3> operator*(double s, const std::array<double, 3> &a)
{
    return a * s;
}

inline std::array<double, 3> operator/(const std::array<double, 3> &a, double s)
{
    return {a[0] / s, a[1] / s, a[2] / s};
}

inline double _dot(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double _norm(const std::array<double, 3> &a)
{
    return std::sqrt(_dot(a, a));
}

inline std::array<double, 3> _cross(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}


} // namespace kep3::linalg

#endif
