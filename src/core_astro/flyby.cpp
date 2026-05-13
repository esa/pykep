// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include <heyoka/expression.hpp>
#include <heyoka/math/cos.hpp>
#include <heyoka/math/sqrt.hpp>


#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>

#include <kep3/core_astro/flyby.hpp>
#include <kep3/linalg.hpp>

namespace kep3
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out,
                                 double mu, double safe_radius)
{
    const double Vin2 = v_rel_in[0] * v_rel_in[0] + v_rel_in[1] * v_rel_in[1] + v_rel_in[2] * v_rel_in[2];
    const double Vout2 = v_rel_out[0] * v_rel_out[0] + v_rel_out[1] * v_rel_out[1] + v_rel_out[2] * v_rel_out[2];
    const double eq_V2 = Vin2 - Vout2;

    const double e_min = 1 + safe_radius / mu * Vin2;
    // Using Giacomo Acciarini's (ACT) differentiable version as anew pykep baseline.
    const double cosalpha = (v_rel_in[0] * v_rel_out[0] + v_rel_in[1] * v_rel_out[1] + v_rel_in[2] * v_rel_out[2])
                            / std::sqrt(Vin2 * Vout2);
    const double ineq_delta = (1. - 2. / e_min / e_min) - cosalpha;
    return {eq_V2, ineq_delta};
}

std::pair<std::vector<heyoka::expression>, std::optional<std::vector<heyoka::expression>>> fb_con(bool jacobian)
{
    using namespace heyoka;
    // The symbolic variables.
    auto [v_rel_in_x, v_rel_in_y, v_rel_in_z, v_rel_out_x, v_rel_out_y, v_rel_out_z] = make_vars("vx_i", "vy_i", "vz_i", "vx_o", "vy_o", "vz_o");
    auto mu = par[0];
    auto safe_radius = heyoka::par[1];

    // The equality constraint
    auto Vin2 = v_rel_in_x * v_rel_in_x + v_rel_in_y * v_rel_in_y + v_rel_in_z * v_rel_in_z;
    auto Vout2 = v_rel_out_x * v_rel_out_x + v_rel_out_y * v_rel_out_y + v_rel_out_z * v_rel_out_z;
    auto eq_V2 = Vin2 - Vout2;  

    // The inequality constraint
    auto e_min = 1_dbl + safe_radius / mu * Vin2;
    auto cosalpha = (v_rel_in_x * v_rel_out_x + v_rel_in_y * v_rel_out_y + v_rel_in_z * v_rel_out_z) / sqrt(Vin2 * Vout2);
    auto ineq_delta = (1_dbl - 2_dbl / e_min / e_min) - cosalpha;
    std::vector<heyoka::expression> retval_1{eq_V2, ineq_delta};

    if (jacobian) {
        auto dt = diff_tensors(retval_1, {v_rel_in_x, v_rel_in_y, v_rel_in_z, v_rel_out_x, v_rel_out_y, v_rel_out_z}, kw::diff_order = 1);
        return {retval_1, dt.get_jacobian()};
    } else {
        return {retval_1, std::nullopt};
    }
}                                                       

std::pair<double, double> fb_con(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out,
                                 const kep3::planet &pl)
{
    return fb_con(v_rel_in, v_rel_out, pl.get_mu_self(), pl.get_safe_radius());
}

double fb_dv(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out, double mu,
             double safe_radius)
{
    const double Vin2 = v_rel_in[0] * v_rel_in[0] + v_rel_in[1] * v_rel_in[1] + v_rel_in[2] * v_rel_in[2];
    const double Vout2 = v_rel_out[0] * v_rel_out[0] + v_rel_out[1] * v_rel_out[1] + v_rel_out[2] * v_rel_out[2];
    // eq_V2 = Vin2 - Vout2;

    const double e_min = 1 + safe_radius / mu * Vin2;
    const double alpha
        = std::acos((v_rel_in[0] * v_rel_out[0] + v_rel_in[1] * v_rel_out[1] + v_rel_in[2] * v_rel_out[2])
                    / std::sqrt(Vin2 * Vout2));
    const double ineq_delta = alpha - 2 * std::asin(1 / e_min);

    double dv = 0.;
    if (ineq_delta > 0.0) {
        dv = std::sqrt(Vout2 + Vin2 - 2.0 * std::sqrt(Vout2 * Vin2) * std::cos(ineq_delta));
    } else {
        dv = std::abs(sqrt(Vout2) - std::sqrt(Vin2));
    }
    return dv;
}

double fb_dv(const std::array<double, 3> &v_rel_in, const std::array<double, 3> &v_rel_out, const kep3::planet &pl)
{
    return fb_dv(v_rel_in, v_rel_out, pl.get_mu_self(), pl.get_safe_radius());
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::array<double, 3> fb_vout(const std::array<double, 3> &v_in, const std::array<double, 3> &v_pla, double rp,
                              double beta, double mu)
{
    std::array<double, 3> v_rel_in = {v_in[0] - v_pla[0], v_in[1] - v_pla[1], v_in[2] - v_pla[2]};
    std::array<double, 3> v_out{};

    const double v_rel_in2 = v_rel_in[0] * v_rel_in[0] + v_rel_in[1] * v_rel_in[1] + v_rel_in[2] * v_rel_in[2];
    const double v_rel_in_norm = std::sqrt(v_rel_in2);
    const double ecc = 1 + rp / mu * v_rel_in2;
    const double delta = 2 * std::asin(1.0 / ecc);
    const std::array<double, 3> i_hat
        = {{v_rel_in[0] / v_rel_in_norm, v_rel_in[1] / v_rel_in_norm, v_rel_in[2] / v_rel_in_norm}};

    std::array<double, 3> j_hat = linalg::_cross(i_hat, v_pla);
    linalg::_normalize(j_hat);
    std::array<double, 3> k_hat = linalg::_cross(i_hat, j_hat);

    v_out[0] = v_pla[0] + v_rel_in_norm * std::cos(delta) * i_hat[0]
               + v_rel_in_norm * std::cos(beta) * std::sin(delta) * j_hat[0]
               + v_rel_in_norm * std::sin(beta) * std::sin(delta) * k_hat[0];
    v_out[1] = v_pla[1] + v_rel_in_norm * std::cos(delta) * i_hat[1]
               + v_rel_in_norm * std::cos(beta) * std::sin(delta) * j_hat[1]
               + v_rel_in_norm * std::sin(beta) * std::sin(delta) * k_hat[1];
    v_out[2] = v_pla[2] + v_rel_in_norm * std::cos(delta) * i_hat[2]
               + v_rel_in_norm * std::cos(beta) * std::sin(delta) * j_hat[2]
               + v_rel_in_norm * std::sin(beta) * std::sin(delta) * k_hat[2];
    return v_out;
}

} // namespace kep3
