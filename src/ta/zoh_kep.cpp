// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <heyoka/kw.hpp>
#include <mutex>
#include <unordered_map>
#include <vector>

#include <heyoka/expression.hpp>
#include <heyoka/math/exp.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/ta/zoh_kep.hpp>

using heyoka::exp;
using heyoka::expression;
using heyoka::make_vars;
using heyoka::par;
using heyoka::pow;
using heyoka::prime;
using heyoka::sum;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
std::vector<std::pair<expression, expression>> zoh_kep_dyn()
{
    // The symbolic variables.
    auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");

    // Renaming parameters.
    const auto &thrust = par[0];
    const auto &ix = par[1];
    const auto &iy = par[2];
    const auto &iz = par[3];
    const auto &c = par[4];

    // The square of the radius.
    const auto r2 = sum({pow(x, 2.), pow(y, 2.), pow(z, 2.)});

    // The equations of motion.
    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot = -pow(r2, -3. / 2.) * x + thrust * ix / m;
    const auto vydot = -pow(r2, -3. / 2.) * y + thrust * iy / m;
    const auto vzdot = -pow(r2, -3. / 2.) * z + thrust * iz / m;
    const auto mdot = -c * thrust * exp(-1. / m / 1e16);

    return {prime(x) = xdot,  prime(y) = ydot,  prime(z) = zdot,  prime(vx) = vxdot,
            prime(vy) = vydot, prime(vz) = vzdot, prime(m) = mdot};
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_kep_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_kep_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_kep(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zoh_kep_mutex);

    // Lookup.
    if (auto it = ta_zoh_kep_cache.find(tol); it == ta_zoh_kep_cache.end()) {
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{zoh_kep_dyn(), init_state, heyoka::kw::tol = tol,
                                              heyoka::kw::pars = {1., 1., 0., 0., 0.}};
        return ta_zoh_kep_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_kep_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_kep_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_kep_var(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zoh_kep_var_mutex);

    // Lookup.
    if (auto it = ta_zoh_kep_var_cache.find(tol); it == ta_zoh_kep_var_cache.end()) {
        auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");
        auto vsys = var_ode_sys(zoh_kep_dyn(), {x, y, z, vx, vy, vz, m, par[0], par[1], par[2], par[3]}, 1);
        // Cache miss, create new one.
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_zoh_kep_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_zoh_kep_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zoh_kep_mutex);
    return ta_zoh_kep_cache.size();
}

size_t get_ta_zoh_kep_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zoh_kep_var_mutex);
    return ta_zoh_kep_var_cache.size();
}

} // namespace kep3::ta
