// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

/// From cartesian to osculating Keplerian
/**
 * Transforms cartesian coordinates (r,v) to Keplerian elements (a,e,i,W,w,E).
 * Note that we use the eccentric anomaly (or Gudermannian if e > 1)
 */

#ifndef kep3_STM_H
#define kep3_STM_H

#include <array>
#include <optional>
#include <utility>

#include <kep3/detail/visibility.hpp>

namespace kep3
{
// From:
// Reynolds, Reid G. "Direct Solution of the Keplerian State Transition Matrix." Journal of Guidance, Control, and
// Dynamics 45, no. 6 (2022): 1162-1165.
kep3_DLL_PUBLIC std::array<double, 36> stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0,
                                                    const std::array<std::array<double, 3>, 2> &pos_vel, double tof,
                                                    double mu = 1.);

// From:
// Lagrange Coefficients and their (manually done) derivatives -- faster than Reynolds, more difficult to
// implement --
kep3_DLL_PUBLIC std::array<double, 36> stm_lagrangian(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof,
                                                      double mu,                           // NOLINT
                                                      double R0, double Rf, double energy, // NOLINT
                                                      double sigma0,                       // NOLINT
                                                      double a, double s0, double c0,      // NOLINT
                                                      double DX, double F, double G, double Ft, double Gt);

// For consistency we offer the same interface we have for propagate lagrangian to access Reynolds stm.
kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_stm_reynolds(const std::array<std::array<double, 3>, 2> &pos_vel0, double tof, double mu = 1.,
                       bool stm = false);

} // namespace kep3
#endif // kep3_IC2mee2IC_H