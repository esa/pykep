// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <mutex>
#include <unordered_map>
#include <vector>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/cos.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sin.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/math/time.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/bcp.hpp>

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
std::vector<std::pair<expression, expression>> bcp_dyn()
{
    // The symbolic variables.
    auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");

    // Renaming parameters.
    const auto &mu = par[0];
    const auto &mu_sun = par[1];
    const auto &rho_sun = par[2];
    const auto &omega_sun = par[3];

    // Distances to the bodies.
    auto r_1 = sqrt(sum({pow(x + par[0], 2.), pow(y, 2.), pow(z, 2.)}));
    auto r_2 = sqrt(sum({pow(x - (1. - par[0]), 2.), pow(y, 2.), pow(z, 2.)}));
    auto x_sun = par[2] * heyoka::cos(omega_sun * heyoka::time);
    auto y_sun = par[2] * heyoka::sin(omega_sun * heyoka::time);
    auto r_sun = sqrt(sum({pow(x - x_sun, 2.), pow(y - y_sun, 2.), pow(z, 2.)}));

    // The Equations of Motion.
    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot = 2. * vy + x - (1. - par[0]) * (x + par[0]) / (pow(r_1, 3.))
                       - par[0] * (x + par[0] - 1.) / pow(r_2, 3.) - mu_sun / pow(r_sun, 3.) * (x - x_sun)
                       - mu_sun / pow(rho_sun, 2.) * heyoka::cos(omega_sun * heyoka::time);
    const auto vydot = -2. * vx + y - (1. - par[0]) * y / pow(r_1, 3.) - par[0] * y / pow(r_2, 3.)
                       - mu_sun / pow(r_sun, 3.) * (y - y_sun)
                       - mu_sun / pow(rho_sun, 2.) * heyoka::sin(omega_sun * heyoka::time);
    const auto vzdot = -(1. - par[0]) * z / pow(r_1, 3.) - par[0] * z / pow(r_2, 3.) - mu_sun / pow(r_sun, 3.) * z;

    return {prime(x) = xdot, prime(y) = ydot, prime(z) = zdot, prime(vx) = vxdot, prime(vy) = vydot, prime(vz) = vzdot};
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_bcp_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_bcp_cache;

const heyoka::taylor_adaptive<double> &get_ta_bcp(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_bcp_mutex);

    // Lookup.
    if (auto it = ta_bcp_cache.find(tol); it == ta_bcp_cache.end()) {
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{bcp_dyn(), init_state, heyoka::kw::tol = tol};
        return ta_bcp_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_bcp_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_bcp_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_bcp_var(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_bcp_var_mutex);

    // Lookup.
    if (auto it = ta_bcp_var_cache.find(tol); it == ta_bcp_var_cache.end()) {
        auto [x, y, z, vx, vy, vz] = make_vars("x", "y", "z", "vx", "vy", "vz");
        auto vsys = var_ode_sys(bcp_dyn(), {x, y, z, vx, vy, vz}, 1);
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{vsys, init_state, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_bcp_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_bcp_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_bcp_mutex);
    return ta_bcp_cache.size();
}

size_t get_ta_bcp_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_bcp_mutex);
    return ta_bcp_var_cache.size();
}

} // namespace kep3::ta