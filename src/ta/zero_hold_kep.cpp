// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <heyoka/kw.hpp>
#include <mutex>
#include <unordered_map>
#include <vector>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/relational.hpp>
#include <heyoka/math/select.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/zero_hold_kep.hpp>

using heyoka::eq;
using heyoka::expression;
using heyoka::make_vars;
using heyoka::par;
using heyoka::pow;
using heyoka::prime;
using heyoka::select;
using heyoka::sqrt;
using heyoka::sum;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
std::vector<std::pair<expression, expression>> zero_hold_kep_dyn()
{
    // The symbolic variables
    auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");

    // Renaming parameters
    const auto &mu = par[0];
    const auto &veff = par[1];
    const auto &[Tx, Ty, Tz] = std::array{par[2], par[3], par[4]};

    // The square of the radius
    const auto r2 = sum({pow(x, 2.), pow(y, 2.), pow(z, 2.)});

    // The thrust magnitude
    const auto u_norm = sqrt(sum({pow(Tx, 2.), pow(Ty, 2.), pow(Tz, 2.)}));

    // The Equations of Motion
    const auto xdot = vx;
    const auto ydot = vy;
    const auto zdot = vz;
    const auto vxdot = -mu * pow(r2, -3. / 2) * x + Tx / m;
    const auto vydot = -mu * pow(r2, -3. / 2) * y + Ty / m;
    const auto vzdot = -mu * pow(r2, -3. / 2) * z + Tz / m;
    // To avoid singularities in the corner case u_norm=0. we use a select here. Implications on performances should be
    // studied.
    const auto mdot = select(eq(u_norm, 0.), 0., -u_norm / veff);
    return {prime(x) = xdot,   prime(y) = ydot,   prime(z) = zdot, prime(vx) = vxdot,
            prime(vy) = vydot, prime(vz) = vzdot, prime(m) = mdot};
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zero_hold_kep_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zero_hold_kep_cache;

const heyoka::taylor_adaptive<double> &get_ta_zero_hold_kep(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zero_hold_kep_mutex);

    // Lookup.
    if (auto it = ta_zero_hold_kep_cache.find(tol); it == ta_zero_hold_kep_cache.end()) {
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{zero_hold_kep_dyn(), init_state, heyoka::kw::tol = tol,
                                              heyoka::kw::pars = {1., 1., 0., 0., 0.}};
        return ta_zero_hold_kep_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zero_hold_kep_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zero_hold_kep_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_zero_hold_kep_var(double tol)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zero_hold_kep_var_mutex);

    // Lookup.
    if (auto it = ta_zero_hold_kep_var_cache.find(tol); it == ta_zero_hold_kep_var_cache.end()) {
        auto [x, y, z, vx, vy, vz, m] = make_vars("x", "y", "z", "vx", "vy", "vz", "m");
        auto vsys = var_ode_sys(zero_hold_kep_dyn(), {x, y, z, vx, vy, vz, m, par[2], par[3], par[4]}, 1);
        // Cache miss, create new one.
        const std::vector init_state = {1., 1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{vsys, init_state, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true,
                                              heyoka::kw::pars = {1., 1., 0., 0., 1e-32}};
        return ta_zero_hold_kep_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_zero_hold_kep_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zero_hold_kep_mutex);
    return ta_zero_hold_kep_cache.size();
}

size_t get_ta_zero_hold_kep_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_zero_hold_kep_mutex);
    return ta_zero_hold_kep_var_cache.size();
}

} // namespace kep3::ta