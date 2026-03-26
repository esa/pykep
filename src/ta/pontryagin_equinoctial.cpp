// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <heyoka/kw.hpp>
#include <map>
#include <mutex>
#include <tuple>
#include <vector>

#include <heyoka/config.hpp>
#include <heyoka/expression.hpp>
#include <heyoka/math/cos.hpp>
#include <heyoka/math/log.hpp>
#include <heyoka/math/pow.hpp>
#include <heyoka/math/sin.hpp>
#include <heyoka/math/sqrt.hpp>
#include <heyoka/math/sum.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/pontryagin_equinoctial.hpp>

using heyoka::cos;
using heyoka::diff;
using heyoka::expression;
using heyoka::log;
using heyoka::make_vars;
using heyoka::par;
using heyoka::sin;
using heyoka::sqrt;
using heyoka::taylor_adaptive;
using heyoka::var_ode_sys;

namespace kep3::ta::detail
{
// Custom hash function for std::pair<double, kep3::optimality_type>
struct pair_hash {
    std::size_t operator()(const std::pair<double, kep3::optimality_type> &p) const
    {
        std::size_t h1 = std::hash<double>{}(p.first);                 // Hash the double
        std::size_t h2 = std::hash<int>{}(static_cast<int>(p.second)); // Hash the enum
        return h1 ^ (h2 << 1);                                         // Combine hashes (XOR and shift)
    }
};
} // namespace kep3::ta::detail

namespace kep3::ta
{
// All the relevant expression in the TPBVP of the Equinoctial Low-Thrust indirect OCP are created here.
std::tuple<std::vector<std::pair<expression, expression>>, expression, expression, expression, std::vector<expression>,
           std::vector<expression>>
peq_expression_factory(kep3::optimality_type optimality)
{
    // The state
    auto [p, f, g, h, k, L, m] = make_vars("p", "f", "g", "h", "k", "L", "m");

    // The costate
    auto [lp, lf, lg, lh, lk, lL, lm] = make_vars("lp", "lf", "lg", "lh", "lk", "lL", "lm");
    // The controls
    auto [u, i_r, i_t, i_n] = make_vars("u", "ir", "it", "in");

    // Useful expressions

    auto w = 1. + f * cos(L) + g * sin(L);
    auto s2 = 1. + h * h + k * k;

    // B, matrix vector 6, 3
    std::array<std::array<expression, 3>, 6> B
        = {{{expression(0.), 2. * p / w, expression(0.)},

            {sin(L), ((1. + w) * cos(L) + f) / w, -g / w * (h * sin(L) - k * cos(L))},

            {-cos(L), ((1. + w) * sin(L) + g) / w, f / w * (h * sin(L) - k * cos(L))},

            {expression(0.), expression(0.), s2 / w / 2.0 * cos(L)},

            {expression(0.), expression(0.), s2 / w / 2.0 * sin(L)},

            {expression(0.), expression(0.), 1.0 / w * (h * sin(L) - k * cos(L))}}};
    // B = B * sqrt(p / par[0])
    for (auto &inner_arr : B) {
        for (auto &element : inner_arr) {
            element *= sqrt(p / par[0]);
        }
    }

    // i_vers, column vector 3, 1
    std::array<std::array<expression, 1>, 3> i_vers = {{{i_r}, {i_t}, {i_n}}};

    // B * i_vers * u * hy.par[1] / m
    std::array<expression, 6> fx;
    for (size_t i = 0; i < 6u; ++i) {
        fx[i] = (B[i][0] * i_vers[0][0] + B[i][1] * i_vers[1][0] + B[i][2] * i_vers[2][0]) * u * par[1] / m;
    }

    // fx + D, fm
    fx[5] = fx[5] + sqrt(par[0] / p / p / p) * w * w;
    auto fm = -par[1] / par[2] * u;

    // BTlam = B.T@lx
    std::array<expression, 6> lx = {lp, lf, lg, lh, lk, lL};
    std::array<expression, 3> BTlam;
    for (auto i = 0u; i < 3u; ++i) {
        BTlam[i] = expression(0.);
        for (auto j = 0u; j < 6u; ++j) {
            BTlam[i] += B[j][i] * lx[j];
        }
    }
    auto BTlam_norm = sqrt(BTlam[0] * BTlam[0] + BTlam[1] * BTlam[1] + BTlam[2] * BTlam[2]);

    expression rho{};
    expression H_full{};
    std::map<expression, expression> argmin_H_full{};

    if (optimality == kep3::optimality_type::MASS) {
        // Hamiltonian (mass optimal with log barrier)
        H_full = lx[0] * fx[0] + lx[1] * fx[1] + lx[2] * fx[2] + lx[3] * fx[3] + lx[4] * fx[4] + lx[5] * fx[5] + lm * fm
                 + par[4] * par[1] / par[2] * (u - par[3] * log(u * (1. - u)));
        // Switching function (mass optimal with log barrier)
        rho = 1. - par[2] * BTlam_norm / m / par[4] - lm / par[4];

        // We apply Pontryagin minimum principle (primer vector and u ^ * = 2eps / (rho + 2eps + sqrt(rho ^ 2 + 4 * eps
        // ^ 2)))
        argmin_H_full = {
            {i_r, -BTlam[0] / BTlam_norm},
            {i_t, -BTlam[1] / BTlam_norm},
            {i_n, -BTlam[2] / BTlam_norm},
            {u, 2. * par[3] / (rho + 2. * par[3] + sqrt(rho * rho + 4. * par[3] * par[3]))},
        };

    } else if (optimality == kep3::optimality_type::TIME) {
        // Hamiltonian (time optimal)
        // Hamiltonian (mass optimal with log barrier)
        H_full = lx[0] * fx[0] + lx[1] * fx[1] + lx[2] * fx[2] + lx[3] * fx[3] + lx[4] * fx[4] + lx[5] * fx[5] + lm * fm
                 + par[4] * par[1] / par[2];
        // Switching function (mass optimal with log barrier)
        rho = -par[2] * BTlam_norm / m / par[4] - lm / par[4];

        // We apply Pontryagin minimum principle (primer vector and u ^ * = 2eps / (rho + 2eps + sqrt(rho ^ 2 + 4 * eps
        // ^ 2)))
        argmin_H_full = {
            {i_r, -BTlam[0] / BTlam_norm},
            {i_t, -BTlam[1] / BTlam_norm},
            {i_n, -BTlam[2] / BTlam_norm},
            {u, expression(1.)},
        };
    }

    // #Augmented equations of motion
    auto rhs = std::vector<expression>{};
    rhs.reserve(14);
    for (const auto &var : {lp, lf, lg, lh, lk, lL, lm, p, f, g, h, k, L, m}) {
        rhs.push_back(diff(H_full, var));
    }
    // The costate equations have the minus sign (Hamiltonian formalism)
    for (auto j = 7u; j < 14; ++j) {
        rhs[j] = -rhs[j];
    }

    rhs = subs(rhs, argmin_H_full);

    // We build the Hamiltonian as a function of the state / co - state only
    //(i.e.no longer of controls)
    auto H = subs(H_full, argmin_H_full);

    // We assemble the system dynamics
    std::vector<std::pair<expression, expression>> sys{};
    sys.reserve(14);
    std::vector<expression> vars{p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm};
    for (auto i = 0u; i < 14; ++i) {
        sys.emplace_back(vars[i], rhs[i]);
    }
    return {sys, H, rho, argmin_H_full[u], {argmin_H_full[i_r], argmin_H_full[i_t], argmin_H_full[i_n]}, rhs};
}

std::vector<std::pair<expression, expression>> peq_dyn(kep3::optimality_type optimality)
{
    return std::get<0>(peq_expression_factory(optimality));
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_peq_mutex;
// Use a Function-Local static Variable (Lazy Initialization)
std::unordered_map<std::pair<double, kep3::optimality_type>, taylor_adaptive<double>, kep3::ta::detail::pair_hash> &
get_ta_peq_cache()
{
    static std::unordered_map<std::pair<double, kep3::optimality_type>, taylor_adaptive<double>,
                              kep3::ta::detail::pair_hash>
        ta_peq_cache;
    return ta_peq_cache;
}

const taylor_adaptive<double> &get_ta_peq(double tol, kep3::optimality_type optimality)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_peq_mutex);
    // Lookup.
    if (auto it = get_ta_peq_cache().find({tol, optimality}); it == get_ta_peq_cache().end()) {
        // Cache miss, create new one.
        auto new_ta = taylor_adaptive<double>{std::get<0>(peq_expression_factory(optimality)), heyoka::kw::tol = tol,
                                              heyoka::kw::compact_mode = true};
        return get_ta_peq_cache()
            .insert(std::make_pair(std::make_pair(tol, optimality), std::move(new_ta)))
            .first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ta_peq_var_mutex;
// Use a Function-Local static Variable (Lazy Initialization)
std::unordered_map<std::pair<double, kep3::optimality_type>, taylor_adaptive<double>, kep3::ta::detail::pair_hash> &
get_ta_peq_var_cache()
{
    static std::unordered_map<std::pair<double, kep3::optimality_type>, taylor_adaptive<double>,
                              kep3::ta::detail::pair_hash>
        ta_peq_cache;
    return ta_peq_cache;
}

const taylor_adaptive<double> &get_ta_peq_var(double tol, kep3::optimality_type optimality)
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_peq_var_mutex);

    // Lookup.
    if (auto it = get_ta_peq_var_cache().find({tol, optimality}); it == get_ta_peq_var_cache().end()) {
        auto [lp, lf, lg, lh, lk, lL, lm] = make_vars("lp", "lf", "lg", "lh", "lk", "lL", "lm");
        auto vsys
            = var_ode_sys(std::get<0>(peq_expression_factory(optimality)), {lp, lf, lg, lh, lk, lL, lm, par[4]}, 1);
        // Cache miss, create new one.
        auto new_ta = taylor_adaptive<double>{vsys, heyoka::kw::tol = tol, heyoka::kw::compact_mode = true};
        return get_ta_peq_var_cache()
            .insert(std::make_pair(std::make_pair(tol, optimality), std::move(new_ta)))
            .first->second;
    } else {
        // Cache hit, return existing.
        return it->second;
    }
}

size_t get_ta_peq_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_peq_mutex);
    return get_ta_peq_cache().size();
}

size_t get_ta_peq_var_cache_dim()
{
    // Lock down for access to cache.
    std::lock_guard const lock(ta_peq_mutex);
    return get_ta_peq_var_cache().size();
}

// Factory functions to help the static variable initialization later
auto peq_H_cfunc_factory(kep3::optimality_type optimality)
{
    auto [p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm]
        = make_vars("p", "f", "g", "h", "k", "L", "m", "lp", "lf", "lg", "lh", "lk", "lL", "lm");
    return heyoka::cfunc<double>({std::get<1>(peq_expression_factory(optimality))},
                                 {p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm});
}
auto peq_SF_cfunc_factory(kep3::optimality_type optimality)
{
    auto [p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm]
        = make_vars("p", "f", "g", "h", "k", "L", "m", "lp", "lf", "lg", "lh", "lk", "lL", "lm");
    return heyoka::cfunc<double>({std::get<2>(peq_expression_factory(optimality))},
                                 {p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm});
}
auto peq_u_cfunc_factory(kep3::optimality_type optimality)
{
    auto [p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm]
        = make_vars("p", "f", "g", "h", "k", "L", "m", "lp", "lf", "lg", "lh", "lk", "lL", "lm");
    return heyoka::cfunc<double>({std::get<3>(peq_expression_factory(optimality))},
                                 {p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm});
}
auto peq_i_vers_cfunc_factory(kep3::optimality_type optimality)
{
    auto [p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm]
        = make_vars("p", "f", "g", "h", "k", "L", "m", "lp", "lf", "lg", "lh", "lk", "lL", "lm");
    return heyoka::cfunc<double>({std::get<4>(peq_expression_factory(optimality))},
                                 {p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm});
}
auto peq_dyn_cfunc_factory(kep3::optimality_type optimality)
{
    auto [p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm]
        = make_vars("p", "f", "g", "h", "k", "L", "m", "lp", "lf", "lg", "lh", "lk", "lL", "lm");
    auto rhs = std::get<5>(peq_expression_factory(optimality));
    return heyoka::cfunc<double>({rhs[0], rhs[1], rhs[2], rhs[3], rhs[4], rhs[5], rhs[13]},
                                 {p, f, g, h, k, L, m, lp, lf, lg, lh, lk, lL, lm});
}

// Function-level static variable: it is initialised the first time the function is invoked
// and the initialisation is guaranteed to be thread-safe (the compiler will put the necessary
// mutex and locking).
const heyoka::cfunc<double> &get_peq_H_cfunc(kep3::optimality_type optimality)
{
    switch (optimality) {
        case kep3::optimality_type::MASS: {
            static const auto peq_H_cfunc_MASS = peq_H_cfunc_factory(kep3::optimality_type::MASS);
            return peq_H_cfunc_MASS;
        }
        case kep3::optimality_type::TIME: {
            static const auto peq_H_cfunc_TIME = peq_H_cfunc_factory(kep3::optimality_type::TIME);
            return peq_H_cfunc_TIME;
        }
        // LCOV_EXCL_START
        default:
            throw std::invalid_argument("Unhandled optimality type in get_peq_H_cfunc()");
            // LCOV_EXCL_STOP
    }
}

const heyoka::cfunc<double> &get_peq_SF_cfunc(kep3::optimality_type optimality)
{
    switch (optimality) {
        case kep3::optimality_type::MASS: {
            static const auto peq_SF_cfunc_MASS = peq_SF_cfunc_factory(kep3::optimality_type::MASS);
            return peq_SF_cfunc_MASS;
        }
        case kep3::optimality_type::TIME: {
            static const auto peq_SF_cfunc_TIME = peq_SF_cfunc_factory(kep3::optimality_type::TIME);
            return peq_SF_cfunc_TIME;
        }
        // LCOV_EXCL_START
        default:
            throw std::invalid_argument("Unhandled optimality type in get_peq_SF_cfunc()");
            // LCOV_EXCL_STOP
    }
}

const heyoka::cfunc<double> &get_peq_u_cfunc(kep3::optimality_type optimality)
{
    switch (optimality) {
        case kep3::optimality_type::MASS: {
            static const auto peq_u_cfunc_MASS = peq_u_cfunc_factory(kep3::optimality_type::MASS);
            return peq_u_cfunc_MASS;
        }
        case kep3::optimality_type::TIME: {
            static const auto peq_u_cfunc_TIME = peq_u_cfunc_factory(kep3::optimality_type::TIME);
            return peq_u_cfunc_TIME;
        }
        // LCOV_EXCL_START
        default:
            throw std::invalid_argument("Unhandled optimality type in get_peq_u_cfunc()");
            // LCOV_EXCL_STOP
    }
}
const heyoka::cfunc<double> &get_peq_i_vers_cfunc(kep3::optimality_type optimality)
{
    switch (optimality) {
        case kep3::optimality_type::MASS: {
            static const auto peq_i_vers_cfunc_MASS = peq_i_vers_cfunc_factory(kep3::optimality_type::MASS);
            return peq_i_vers_cfunc_MASS;
        }
        case kep3::optimality_type::TIME: {
            static const auto peq_i_vers_cfunc_TIME = peq_i_vers_cfunc_factory(kep3::optimality_type::TIME);
            return peq_i_vers_cfunc_TIME;
        }
        // LCOV_EXCL_START
        default:
            throw std::invalid_argument("Unhandled optimality type in get_peq_i_vers_cfunc()");
            // LCOV_EXCL_STOP
    }
}
const heyoka::cfunc<double> &get_peq_dyn_cfunc(kep3::optimality_type optimality)
{
    switch (optimality) {
        case kep3::optimality_type::MASS: {
            static const auto peq_dyn_cfunc_MASS = peq_dyn_cfunc_factory(kep3::optimality_type::MASS);
            return peq_dyn_cfunc_MASS;
        }
        case kep3::optimality_type::TIME: {
            static const auto peq_dyn_cfunc_TIME = peq_dyn_cfunc_factory(kep3::optimality_type::TIME);
            return peq_dyn_cfunc_TIME;
        }
        // LCOV_EXCL_START
        default:
            throw std::invalid_argument("Unhandled optimality type in get_peq_dyn_cfunc()");
            // LCOV_EXCL_STOP
    }
}

} // namespace kep3::ta