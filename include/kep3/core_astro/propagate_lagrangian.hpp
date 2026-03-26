// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_PROPAGATE_LAGRANGIAN_H
#define kep3_PROPAGATE_LAGRANGIAN_H

#include <array>
#include <optional>
#include <utility>
#include <vector>

#include <kep3/core_astro/stm.hpp>
#include <kep3/detail/visibility.hpp>

namespace kep3
{

/// Lagrangian propagation
/**
 * This function propagates an initial Cartesian state for a time t assuming a
 * central body and a keplerian motion. Lagrange coefficients are used as basic
 * numerical technique. All units systems can be used, as long
 * as the input parameters are all expressed in the same system.
 * The information on the state transition matrix can be optionally asked.
 */
kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_lagrangian(const std::array<std::array<double, 3>, 2> &pos_vel, double tof, double mu, bool stm = false);

kep3_DLL_PUBLIC std::vector<std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>>
propagate_lagrangian_grid(const std::array<std::array<double, 3>, 2> &pos_vel, const std::vector<double> &time_grid, double mu,
                       bool stm = false);

// These are backup functions that use a different algorithm to get the same as propagate_lagrangian.
// We offer them with an identical interface even if the stm is not implemented.
kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_lagrangian_u(const std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu, bool = false);

kep3_DLL_PUBLIC std::pair<std::array<std::array<double, 3>, 2>, std::optional<std::array<double, 36>>>
propagate_keplerian(const std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu, bool = false);
} // namespace kep3

#endif // KEP_TOOLBOX_PROPAGATE_LAGRANGIAN_H