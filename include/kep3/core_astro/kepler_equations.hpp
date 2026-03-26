// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_KEPLER_EQUATIONS_H
#define kep3_KEPLER_EQUATIONS_H

#include <cmath>

#include <kep3/core_astro/special_functions.hpp>

namespace kep3
{
// In terms of the eccentric anomaly (E)
// -------------------------------------------
inline double kepE(double E, double M, double ecc)
{
    return (E - ecc * std::sin(E) - M);
}
// Its first derivative
inline double d_kepE(double E, double ecc)
{
    return (1 - ecc * std::cos(E));
}
// And its second derivative
inline double dd_kepE(double E, double ecc)
{
    return ecc * std::sin(E);
}
// -------------------------------------------

// In terms of the hyperbolic anomaly (E)
// -------------------------------------------
inline double kepH(double H, double M, double ecc)
{
    return (ecc * std::sinh(H) - H - M);
}
// Its first derivative
inline double d_kepH(double H, double ecc)
{
    return (ecc * std::cosh(H) - 1);
}
// And its second derivative
inline double dd_kepH(double H, double ecc)
{
    return ecc * std::sinh(H);
}
// -------------------------------------------

// In terms of the eccentric anomaly difference (DE)
// -------------------------------------------
inline double kepDE(double DE, double DM, double sigma0, double sqrta, double a, double R)
{
    return -DM + DE + sigma0 / sqrta * (1 - std::cos(DE)) - (1 - R / a) * std::sin(DE);
}

inline double d_kepDE(double DE, double sigma0, double sqrta, double a, double R)
{
    return 1 + sigma0 / sqrta * std::sin(DE) - (1 - R / a) * std::cos(DE);
}

inline double dd_kepDE(double DE, double sigma0, double sqrta, double a, double R)
{
    return sigma0 / sqrta * std::cos(DE) + (1 - R / a) * std::sin(DE);
}

// In terms of the hyperbolic anomaly difference (DH)
// -------------------------------------------
inline double kepDH(double DH, double DN, double sigma0, double sqrta, double a, double R)
{
    return -DN - DH + sigma0 / sqrta * (std::cosh(DH) - 1) + (1 - R / a) * std::sinh(DH);
}

inline double d_kepDH(double DH, double sigma0, double sqrta, double a, double R)
{
    return -1. + sigma0 / sqrta * std::sinh(DH) + (1 - R / a) * std::cosh(DH);
}

inline double dd_kepDH(double DH, double sigma0, double sqrta, double a, double R)
{
    return sigma0 / sqrta * std::cosh(DH) + (1 - R / a) * std::sinh(DH);
}

// In terms of the universal anomaly difference (DS)
// -------------------------------------------
inline double kepDS(const double &DS, const double &DT, const double &r0, const double &vr0, const double &alpha,
                    const double &mu)
{
    double S = stumpff_s(alpha * DS * DS);
    double C = stumpff_c(alpha * DS * DS);
    double retval
        = -std::sqrt(mu) * DT + r0 * vr0 * DS * DS * C / std::sqrt(mu) + (1 - alpha * r0) * DS * DS * DS * S + r0 * DS;
    return (retval);
}

inline double d_kepDS(const double &DS, const double &r0, const double &vr0, const double &alpha, const double &mu)
{
    double S = stumpff_s(alpha * DS * DS);
    double C = stumpff_c(alpha * DS * DS);
    double retval = r0 * vr0 / std::sqrt(mu) * DS * (1 - alpha * DS * DS * S) + (1 - alpha * r0) * DS * DS * C + r0;
    return (retval);
}

inline double dd_kepDS(const double &DS, const double &r0, const double &vr0, const double &alpha, const double &mu)
{
    double S = stumpff_s(alpha * DS * DS);
    double C = stumpff_c(alpha * DS * DS);
    double retval = r0 * vr0 / std::sqrt(mu) * (1 - alpha * DS * DS * C) + (1 - alpha * r0) * (1 - DS * DS * S);

    return (retval);
}

} // namespace kep3
#endif // kep3_KEPLER_EQUATIONS_H
