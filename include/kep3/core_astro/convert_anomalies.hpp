// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_CONVERT_ANOMALIES_H
#define kep3_CONVERT_ANOMALIES_H

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/kepler_equations.hpp>

namespace kep3
{

// mean to eccentric (only ellipses) e<1. Preserves the sign and integer number
// of revolutions.
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
inline double m2e(double M, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    // We compute the sin and cos of the mean anomaly which are used also later
    // for the initial guess (3rd order expansion of Kepler's equation).
    double sinM = std::sin(M), cosM = std::cos(M);
    // Here we use the atan2 to recover the mean anomaly in the [0,2pi] range.
    // This makes sure that for high value of M no catastrophic cancellation
    // occurs, as would be the case using std::fmod(M, 2pi)
    double M_cropped = std::atan2(sinM, cosM);

    // The Initial guess follows from a third order expansion of Kepler's
    // equation.
    double IG = M_cropped + ecc * sinM + ecc * ecc * sinM * cosM + ecc * ecc * ecc * sinM * (1.5 * cosM * cosM - 0.5);

    const int digits = std::numeric_limits<double>::digits;
    std::uintmax_t max_iter = 100u;

    // Newton-raphson or Halley iterates can be used here. Similar performances,
    // thus we choose the simplest algorithm.
    double sol = boost::math::tools::newton_raphson_iterate(
        [M_cropped, ecc](double E) { return std::make_tuple(kepE(E, M_cropped, ecc), d_kepE(E, ecc)); }, IG,
        IG - kep3::pi, IG + kep3::pi, digits, max_iter);
    if (max_iter == 100u) {
        throw std::domain_error("Maximum number of iterations exceeded when solving Kepler's "
                                "equation for the eccentric anomaly in m2e.");
    }
    return sol;
}
// eccentric to mean (only ellipses) e<1
inline double e2m(double E, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (E - ecc * std::sin(E));
}

// eccentric to true (only ellipses) e<1 (returns in range [-pi,pi])
inline double e2f(double E, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return 2 * std::atan(std::sqrt((1 + ecc) / (1 - ecc)) * std::tan(E / 2));
}

// true to eccentric (only ellipses) e<1 (returns in range [-pi,pi])
inline double f2e(double f, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return 2 * std::atan(std::sqrt((1 - ecc) / (1 + ecc)) * std::tan(f / 2));
}

// mean to true (only ellipses) e<1 (returns in range [-pi,pi])
inline double m2f(double M, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return e2f(m2e(M, ecc), ecc);
}

// true to mean (only ellipses) e<1 (returns in range [-pi,pi])
inline double f2m(double f, double ecc)
{
    if (ecc >= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return e2m(f2e(f, ecc), ecc);
}

// gudermannian to true (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double zeta2f(double f, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return 2 * std::atan(std::sqrt((1 + ecc) / (ecc - 1)) * std::tan(f / 2));
}

// true to gudermannian (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double f2zeta(double zeta, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return 2 * std::atan(std::sqrt((ecc - 1) / (1 + ecc)) * std::tan(zeta / 2));
}

// mean to hyperbolic (only hyperbolas) e>1.
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
inline double n2h(double N, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    // The Initial guess (TODO(darioizo) improve)
    double IG = 1.;

    const int digits = std::numeric_limits<double>::digits;
    std::uintmax_t max_iter = 100u;

    // Newton-raphson iterates.
    double sol = boost::math::tools::newton_raphson_iterate(
        [N, ecc](double H) { return std::make_tuple(kepH(H, N, ecc), d_kepH(H, ecc)); }, IG, IG - 20 * kep3::pi,
        IG + 20 * kep3::pi, digits, max_iter);
    if (max_iter == 100u) {
        throw std::domain_error("Maximum number of iterations exceeded when solving Kepler's "
                                "equation for the hyperbolic anomaly in m2h.");
    }
    return sol;
}

// hyperbolic H to hyperbolic mean N (only hyperbolas) e>1
inline double h2n(double H, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return (ecc * std::sinh(H) - H);
}

// hyperbolic H to true (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double h2f(double H, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return 2 * std::atan(std::sqrt((1 + ecc) / (ecc - 1)) * std::tanh(H / 2));
}

// true to hyperbolic H (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double f2h(double f, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return 2 * std::atanh(std::sqrt((ecc - 1) / (1 + ecc)) * std::tan(f / 2));
}

// mean hyperbolic to true (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double n2f(double N, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return h2f(n2h(N, ecc), ecc);
}

// true to mean hyperbolic (only hyperbolas) e>1 (returns in range [-pi,pi])
inline double f2n(double f, double ecc)
{
    if (ecc <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    return h2n(f2h(f, ecc), ecc);
}

} // namespace kep3
#endif // kep3_TOOLBOX_M2E_H
