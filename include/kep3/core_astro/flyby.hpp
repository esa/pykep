// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_FLYBY_H
#define kep3_FLYBY_H

#include <array>
#include <utility>

#include <kep3/detail/visibility.hpp>
#include <kep3/planet.hpp>

namespace kep3
{

// Returns the constraints [v2, alpha] of a fly-by. The first is an equality constraint, the second an inequality
// (negative if satisfied).
kep3_DLL_PUBLIC std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in,
                                                 const std::array<double, 3> &v_rel_out, double mu, double safe_radius);

kep3_DLL_PUBLIC std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in,
                                                 const std::array<double, 3> &v_rel_out, const kep3::planet &pl);

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