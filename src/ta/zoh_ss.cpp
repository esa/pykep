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
#include <heyoka/math/cos.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sin.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/ta/zoh_ss.hpp>

using heyoka::cos;
using heyoka::expression;
using heyoka::make_vars;
using heyoka::par;
using heyoka::pow;
using heyoka::prime;
using heyoka::sin;
using heyoka::sqrt;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
std::vector<std::pair<expression, expression>> zoh_ss_dyn()
{
    // State.
    auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");

    // Controls/parameters.
    const auto &alpha = par[0];
    const auto &beta = par[1];
    const auto &c = par[2];

    // Auxiliary expressions.
    const auto r2 = x * x + y * y + z * z;
    const auto r = sqrt(r2);

    // Angular momentum h = r x v.
    const auto h0 = y * vz - z * vy;
    const auto h1 = z * vx - x * vz;
    const auto h2 = x * vy - y * vx;
    const auto hnorm = sqrt(h0 * h0 + h1 * h1 + h2 * h2);

    // Unit vectors of RTN frame.
    const auto ir0 = x / r;
    const auto ir1 = y / r;
    const auto ir2 = z / r;

    const auto ih0 = h0 / hnorm;
    const auto ih1 = h1 / hnorm;
    const auto ih2 = h2 / hnorm;

    const auto it0 = ih1 * ir2 - ih2 * ir1;
    const auto it1 = ih2 * ir0 - ih0 * ir2;
    const auto it2 = ih0 * ir1 - ih1 * ir0;

    const auto thrust = c * (1. / r2) * pow(cos(alpha), 2.);
    const auto ar = cos(alpha) * thrust;
    const auto at = sin(alpha) * sin(beta) * thrust;
    const auto ah = sin(alpha) * cos(beta) * thrust;

    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot = -x * pow(r2, -3. / 2.) + ar * ir0 + at * it0 + ah * ih0;
    const auto vydot = -y * pow(r2, -3. / 2.) + ar * ir1 + at * it1 + ah * ih1;
    const auto vzdot = -z * pow(r2, -3. / 2.) + ar * ir2 + at * it2 + ah * ih2;

    return {
        prime(x) = xdot,
        prime(y) = ydot,
        prime(z) = zdot,
        prime(vx) = vxdot,
        prime(vy) = vydot,
        prime(vz) = vzdot,
    };
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_ss_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_ss_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_ss(double tol)
{
    std::lock_guard const lock(ta_zoh_ss_mutex);

    if (auto it = ta_zoh_ss_cache.find(tol); it == ta_zoh_ss_cache.end()) {
        const std::vector init_state = {1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{zoh_ss_dyn(), init_state, heyoka::kw::tol = tol, heyoka::kw::pars = {0., 0., 0.}};
        return ta_zoh_ss_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_ss_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_ss_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_ss_var(double tol)
{
    std::lock_guard const lock(ta_zoh_ss_var_mutex);

    if (auto it = ta_zoh_ss_var_cache.find(tol); it == ta_zoh_ss_var_cache.end()) {
        auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");
        auto vsys = var_ode_sys(zoh_ss_dyn(), {x, y, z, vx, vy, vz, par[0], par[1]}, 1);
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_zoh_ss_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

size_t get_ta_zoh_ss_cache_dim()
{
    std::lock_guard const lock(ta_zoh_ss_mutex);
    return ta_zoh_ss_cache.size();
}

size_t get_ta_zoh_ss_var_cache_dim()
{
    std::lock_guard const lock(ta_zoh_ss_var_mutex);
    return ta_zoh_ss_var_cache.size();
}

} // namespace kep3::ta
