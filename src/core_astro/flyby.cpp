// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <cmath>
#include <utility>

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
    const double alpha
        = std::acos((v_rel_in[0] * v_rel_out[0] + v_rel_in[1] * v_rel_out[1] + v_rel_in[2] * v_rel_out[2])
                    / std::sqrt(Vin2 * Vout2));
    const double ineq_delta = alpha - 2 * std::asin(1 / e_min);
    return {eq_V2, ineq_delta};
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
