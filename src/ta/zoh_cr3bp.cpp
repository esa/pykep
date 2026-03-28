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
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/ta/zoh_cr3bp.hpp>

using heyoka::exp;
using heyoka::expression;
using heyoka::make_vars;
using heyoka::par;
using heyoka::pow;
using heyoka::prime;
using heyoka::sqrt;
using heyoka::sum;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
std::vector<std::pair<expression, expression>> zoh_cr3bp_dyn()
{
    // The symbolic variables.
    auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");

    // Controls.
    const auto &thrust = par[0];
    const auto &ix = par[1];
    const auto &iy = par[2];
    const auto &iz = par[3];

    // Parameters.
    const auto &c = par[4];
    const auto &mu = par[5];

    // Distances to the bodies.
    const auto r_1 = sqrt(sum({pow(x + mu, 2.), pow(y, 2.), pow(z, 2.)}));
    const auto r_2 = sqrt(sum({pow(x - (1. - mu), 2.), pow(y, 2.), pow(z, 2.)}));

    // Equations of motion.
    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot
        = 2. * vy + x - (1. - mu) * (x + mu) / pow(r_1, 3.) - mu * (x + mu - 1.) / pow(r_2, 3.) + thrust * ix / m;
    const auto vydot = -2. * vx + y - (1. - mu) * y / pow(r_1, 3.) - mu * y / pow(r_2, 3.) + thrust * iy / m;
    const auto vzdot = -(1. - mu) * z / pow(r_1, 3.) - mu * z / pow(r_2, 3.) + thrust * iz / m;
    const auto mdot = -c * thrust * exp(-1. / m / 1e16);

    return {prime(x) = xdot, prime(y) = ydot, prime(z) = zdot, prime(vx) = vxdot,
            prime(vy) = vydot, prime(vz) = vzdot, prime(m) = mdot};
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_cr3bp_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_cr3bp_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_cr3bp(double tol)
{
    std::lock_guard const lock(ta_zoh_cr3bp_mutex);

    if (auto it = ta_zoh_cr3bp_cache.find(tol); it == ta_zoh_cr3bp_cache.end()) {
        const std::vector init_state = {1., 1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{zoh_cr3bp_dyn(), init_state, heyoka::kw::tol = tol,
                                              heyoka::kw::pars = {1., 1., 0., 0., 0., .01}};
        return ta_zoh_cr3bp_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_cr3bp_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_cr3bp_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_cr3bp_var(double tol)
{
    std::lock_guard const lock(ta_zoh_cr3bp_var_mutex);

    if (auto it = ta_zoh_cr3bp_var_cache.find(tol); it == ta_zoh_cr3bp_var_cache.end()) {
        auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");
        auto vsys = var_ode_sys(zoh_cr3bp_dyn(), {x, y, z, vx, vy, vz, m, par[0], par[1], par[2], par[3]}, 1);
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_zoh_cr3bp_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

size_t get_ta_zoh_cr3bp_cache_dim()
{
    std::lock_guard const lock(ta_zoh_cr3bp_mutex);
    return ta_zoh_cr3bp_cache.size();
}

size_t get_ta_zoh_cr3bp_var_cache_dim()
{
    std::lock_guard const lock(ta_zoh_cr3bp_var_mutex);
    return ta_zoh_cr3bp_var_cache.size();
}

} // namespace kep3::ta
