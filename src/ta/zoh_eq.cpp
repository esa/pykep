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
#include <heyoka/math/exp.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sin.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/ta/zoh_eq.hpp>

using heyoka::cos;
using heyoka::exp;
using heyoka::expression;
using heyoka::make_vars;
using heyoka::par;
using heyoka::prime;
using heyoka::sin;
using heyoka::sqrt;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta
{
std::vector<std::pair<expression, expression>> zoh_eq_dyn()
{
    // The state.
    auto [p, f, g, h, k, L, m] = make_vars("p", "f", "g", "h", "k", "L", "m");

    // Parameters.
    const auto &thrust = par[0];
    const auto &i_r = par[1];
    const auto &i_t = par[2];
    const auto &i_n = par[3];
    const auto &c = par[4];

    const auto w = 1. + f * cos(L) + g * sin(L);
    const auto s2 = 1. + h * h + k * k;
    const auto hsk = h * sin(L) - k * cos(L);
    const auto sqrt_p = sqrt(p);

    // B(x) * T_vector / m + D(x) in expanded form.
    const auto dp = sqrt_p * (2. * p / w) * (i_t * thrust) / m;
    const auto df = sqrt_p
                    * (i_r * thrust * sin(L) + ((1. + w) * cos(L) + f) / w * (i_t * thrust)
                       - g / w * hsk * (i_n * thrust))
                    / m;
    const auto dg = sqrt_p
                    * (-i_r * thrust * cos(L) + ((1. + w) * sin(L) + g) / w * (i_t * thrust)
                       + f / w * hsk * (i_n * thrust))
                    / m;
    const auto dh = sqrt_p * (s2 / w / 2.) * cos(L) * (i_n * thrust) / m;
    const auto dk = sqrt_p * (s2 / w / 2.) * sin(L) * (i_n * thrust) / m;
    const auto dL = sqrt_p * (1. / w * hsk) * (i_n * thrust) / m + sqrt(1. / p / p / p) * w * w;
    const auto dm = -c * thrust * exp(-1. / m / 1e16);

    return {prime(p) = dp, prime(f) = df, prime(g) = dg, prime(h) = dh,
            prime(k) = dk, prime(L) = dL, prime(m) = dm};
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_eq_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_eq_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_eq(double tol)
{
    std::lock_guard const lock(ta_zoh_eq_mutex);

    if (auto it = ta_zoh_eq_cache.find(tol); it == ta_zoh_eq_cache.end()) {
        const std::vector init_state = {1., 1., 1., 1., 1., 1., 1.};
        auto new_ta = taylor_adaptive<double>{zoh_eq_dyn(), init_state, heyoka::kw::tol = tol,
                                              heyoka::kw::pars = {1., 1., 0., 0., 0.}};
        return ta_zoh_eq_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_zoh_eq_var_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unordered_map<double, taylor_adaptive<double>> ta_zoh_eq_var_cache;

const heyoka::taylor_adaptive<double> &get_ta_zoh_eq_var(double tol)
{
    std::lock_guard const lock(ta_zoh_eq_var_mutex);

    if (auto it = ta_zoh_eq_var_cache.find(tol); it == ta_zoh_eq_var_cache.end()) {
        auto [p, f, g, h, k, L, m] = make_vars("p", "f", "g", "h", "k", "L", "m");
        auto vsys = var_ode_sys(zoh_eq_dyn(), {p, f, g, h, k, L, m, par[0], par[1], par[2], par[3]}, 1);
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return ta_zoh_eq_var_cache.insert(std::make_pair(tol, std::move(new_ta))).first->second;
    } else {
        return it->second;
    }
}

size_t get_ta_zoh_eq_cache_dim()
{
    std::lock_guard const lock(ta_zoh_eq_mutex);
    return ta_zoh_eq_cache.size();
}

size_t get_ta_zoh_eq_var_cache_dim()
{
    std::lock_guard const lock(ta_zoh_eq_var_mutex);
    return ta_zoh_eq_var_cache.size();
}

} // namespace kep3::ta
