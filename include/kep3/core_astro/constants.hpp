// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_CONSTANTS_H
#define kep3_CONSTANTS_H

#include <boost/math/constants/constants.hpp>

namespace kep3
{

enum elements_type {
    KEP_M,  // Keplerian Osculating (with Mean Anomaly)
    KEP_F,  // Keplerian Osculating (with True Anomaly)
    MEE,    // Modified Equinoctial Elements
    MEE_R,  // Modified Equinoctial Elements (retrogade)
    POSVEL, // position and Velocity
};

enum optimality_type {
    MASS,  // Mass Optimality
    TIME,  // Time Optimality
};

inline constexpr double pi = boost::math::constants::pi<double>();
inline constexpr double half_pi = boost::math::constants::half_pi<double>();
inline constexpr double AU = 149597870700.0;                 // Astronomical Unit (m) - IAU 2012 Resolution B1
inline constexpr double CAVENDISH = 6.67430e-11;             // Cavendish constant (m^3/s^2/kg)
inline constexpr double MU_SUN = 1.32712440041279419e20;     // Sun's gravitational parameter (m^3/s^2 kg) - DE440
inline constexpr double MU_EARTH = 398600435507000.0;        // Earth's gravitational parameter (m^3/s^2 kg) - DE440
inline constexpr double MU_MOON =  4902800118000.0;          // Moon's gravitational parameter (m^3/s^2 kg) - DE440
inline constexpr double EARTH_VELOCITY = 29784.69183430911;  // Earth's velocity. (m/s)
inline constexpr double EARTH_J2 = 1.08262668E-03;
inline constexpr double EARTH_RADIUS = 6378137; // Earth's radius (m)
inline constexpr double DEG2RAD = (pi / 180.0);
inline constexpr double RAD2DEG = (180.0 / pi);
inline constexpr double DAY2SEC = 86400.0;
inline constexpr double SEC2DAY = (1. / DAY2SEC);
inline constexpr double YEAR2DAY = (365.25);
inline constexpr double DAY2YEAR = (1. / YEAR2DAY);
inline constexpr double G0 = 9.80665; // Acceleration at Earth's surface (m/s^2)
inline constexpr double CR3BP_MU_EARTH_MOON = MU_MOON / (MU_MOON + MU_EARTH); // M_moon / (M_moon + M_Earth) 
inline constexpr double BCP_MU_EARTH_MOON = MU_MOON / (MU_MOON + MU_EARTH); // M_moon / (M_moon + M_Earth) 
inline constexpr double BCP_MU_S = MU_SUN / (MU_MOON + MU_EARTH); // M_Sun / (M_Earth + M_Moon)
inline constexpr double BCP_RHO_S = 3.88811143E2; // Scaled Sun–(Earth + Moon) distance (from Simo' BCP paper)
inline constexpr double BCP_OMEGA_S = -9.25195985E-01; // Scaled angular velocity of the Sun (from Simo' BCP paper)
} // namespace kep3

#endif // kep3_CONSTANTS_H