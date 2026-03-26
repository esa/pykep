// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

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

// Linear algebra helpers to speed up xtensor when, small, fixed size matrices and vectors are involved
// ---------------------------------------------------------------------------------------
template <std::size_t m, std::size_t k, std::size_t n>
xt::xtensor_fixed<double, xt::xshape<m, n>> _dot(const xt::xtensor_fixed<double, xt::xshape<m, k>> &A,
                                                 const xt::xtensor_fixed<double, xt::xshape<k, n>> &B)
{
    xt::xtensor_fixed<double, xt::xshape<m, n>> C{};
    for (decltype(m) i = 0u; i < m; ++i) {
        for (decltype(n) j = 0u; j < n; ++j) {
            C(i, j) = 0;
            for (decltype(k) l = 0u; l < k; ++l) {
                {
                    C(i, j) += A(i, l) * B(l, j);
                }
            }
        }
    }
    return C;
}
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
