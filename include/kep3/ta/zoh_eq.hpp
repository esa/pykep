// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_TA_zoh_eq_H
#define kep3_TA_zoh_eq_H

#include <vector>

#include <kep3/detail/visibility.hpp>

#include <heyoka/expression.hpp>
#include <heyoka/taylor.hpp>

namespace kep3::ta
{
// Returns the ZOH equinoctial low-thrust dynamics.
// 7 states, 5 parameters: [T, i_r, i_t, i_n, c].
kep3_DLL_PUBLIC std::vector<std::pair<heyoka::expression, heyoka::expression>> zoh_eq_dyn();

// These return const references to function level static variables of type heyoka::taylor_adaptive<double>.
// NOTE: The object returned are expected to be copied to then be modified.
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &get_ta_zoh_eq(double tol);
kep3_DLL_PUBLIC const heyoka::taylor_adaptive<double> &
get_ta_zoh_eq_var(double tol); // variational wrt (p,f,g,h,k,L,m,T,ir,it,in), first order.

// Methods to access the cache dimensions.
kep3_DLL_PUBLIC size_t get_ta_zoh_eq_cache_dim();
kep3_DLL_PUBLIC size_t get_ta_zoh_eq_var_cache_dim();
} // namespace kep3::ta

#endif
