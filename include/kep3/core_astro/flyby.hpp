// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_FLYBY_H
#define kep3_FLYBY_H

#include <array>
#include <utility>
#include <optional>

#include <heyoka/expression.hpp>

#include <kep3/detail/visibility.hpp>
#include <kep3/planet.hpp>

namespace kep3
{

// Returns the fly-by constraints [eq_V2, ineq_delta].
// eq_V2 = Vin2 - Vout2 (equality constraint).
// ineq_delta = (1 - 2 / e_min^2) - cos(alpha) with e_min = 1 + safe_radius / mu * Vin2
// (inequality constraint, negative if satisfied).
kep3_DLL_PUBLIC std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in,
                                                 const std::array<double, 3> &v_rel_out, double mu, double safe_radius);

kep3_DLL_PUBLIC std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in,
                                                 const std::array<double, 3> &v_rel_out, const kep3::planet &pl);

// Returns symbolic expressions [eq_V2, ineq_delta] and, optionally, their Jacobian wrt
// [vx_i, vy_i, vz_i, vx_o, vy_o, vz_o].
kep3_DLL_PUBLIC std::pair<std::vector<heyoka::expression>, std::optional<std::vector<heyoka::expression>>> fb_con(bool jacobian = false);


// Returns the dv needed to make a fly-by feasible. (assuming one DV at the out conditions).
kep3_DLL_PUBLIC double fb_dv(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out, double mu,
                             double safe_radius);

kep3_DLL_PUBLIC double fb_dv(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out,
                             const kep3::planet &pl);

// Propagates an absolute incoming velocity through a fly-by.
kep3_DLL_PUBLIC std::array<double, 3> fb_vout(const std::array<double, 3> &v_in, const std::array<double, 3> &v_pla,
                                              double rp, double beta, double mu);

} // namespace kep3
#endif // kep3_FLYBY_H