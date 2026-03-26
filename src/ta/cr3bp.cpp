// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <mutex>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/cr3bp.hpp>

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
std::tuple<std::vector<std::pair<expression, expression>>, expression, expression> expression_factory()
{
    // The symbolic variables.
    auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");

    // Renaming parameters.
    const auto &mu = par[0];

    // Distances to the bodies.
    auto r_1 = sqrt(sum({pow(x + par[0], 2.), pow(y, 2.), pow(z, 2.)}));
    auto r_2 = sqrt(sum({pow(x - (1. - par[0]), 2.), pow(y, 2.), pow(z, 2.)}));

    // The Equations of Motion.
    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot
        = 2. * vy + x - (1. - par[0]) * (x + par[0]) / (pow(r_1, 3.)) - par[0] * (x + par[0] - 1.) / pow(r_2, 3.);
    const auto vydot = -2. * vx + y - (1. - par[0]) * y / pow(r_1, 3.) - par[0] * y / pow(r_2, 3.);
    const auto vzdot = -(1. - par[0]) * z / pow(r_1, 3.) - par[0] * z / pow(r_2, 3.);

    // The effective potential. (note the sign convention here)
    const auto U = 1. / 2. * (pow(x, 2.) + pow(y, 2.)) + (1. - par[0]) / r_1 + par[0] / r_2;
    // The velocity squared (in rotating).
    const auto v2 = (pow(vx, 2.) + pow(vy, 2.) + pow(vz, 2.));
    // The Jacobi constant.
    const auto C = 2. * U - v2;

    return {
        {prime(x) = xdot, prime(y) = ydot, prime(z) = zdot, prime(vx) = vxdot, prime(vy) = vydot, prime(vz) = vzdot},
        U,
        C};
};

std::vector<std::pair<expression, expression>> cr3bp_dyn()
{
    return std::get<0>(expression_factory());
}

heyoka::expression cr3bp_effective_potential_U()
{
    return std::get<1>(expression_factory());
}

heyoka::expression cr3bp_jacobi_C()
{
    return std::get<2>(expression_factory());
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_cr3bp_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_cr3bp_cache;

const heyoka::taylor_adaptive<double> &get_ta_cr3bp(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_cr3bp_mutex);

    // Lookup.
    if (auto it = ta_cr3bp_cache.find(tol); it == ta_cr3bp_cache.end()) {
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{cr3bp_dyn(), init_state, heyoka::kw::tol = tol};
        return ta_cr3bp_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_cr3bp_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_cr3bp_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_cr3bp_var(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_cr3bp_var_mutex);

    // Lookup.
    if (auto it = ta_cr3bp_var_cache.find(tol); it == ta_cr3bp_var_cache.end()) {
        auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");
        auto vsys = var_ode_sys(cr3bp_dyn(), {x, y, z, vx, vy, vz}, 1);
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{vsys, init_state, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_cr3bp_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_cr3bp_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_cr3bp_mutex);
    return ta_cr3bp_cache.size();
}

size_t get_ta_cr3bp_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_cr3bp_mutex);
    return ta_cr3bp_var_cache.size();
}

} // namespace kep3::ta