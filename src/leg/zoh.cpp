// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <heyoka/kw.hpp>
#include <heyoka/taylor.hpp>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/leg/zoh.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg
{

namespace
{

using kep3::linalg::_dot;

bool propagate_until_safe_impl(heyoka::taylor_adaptive<double> &ta, double t, const std::optional<unsigned> &max_steps);

template <std::size_t D, std::size_t C, typename MatDD, typename MatDC, typename MatD1>
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> compute_mc_grad_fixed_impl(
    heyoka::taylor_adaptive<double> &ta_var, const std::vector<double> &state0, const std::vector<double> &controls,
    const std::vector<double> &state1, const std::vector<double> &tgrid, const std::vector<double> &ic_var,
    const std::vector<double> &pars_no_control, const std::optional<unsigned> &max_steps,
    const heyoka::cfunc<double> &dyn_cfunc, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck)
{
    // Prepare output containers.
    std::vector<double> dmc_dx0(D * D, 0.0);
    std::vector<double> dmc_dx1(D * D, 0.0);
    std::vector<double> dmc_dcontrols(D * C * nseg, 0.0);
    std::vector<double> dmc_dtgrid(D * (nseg + 1u), 0.0);

    // Views/shapes used to map flat state blocks and flattened outputs.
    const std::array<std::size_t, 2u> mat_dd_shape{D, D};
    const std::array<std::size_t, 2u> stm_shape{D, D + C};

    // Forward segments.
    std::vector<MatDD> M_seg_fwd;
    std::vector<MatDC> C_seg_fwd;
    std::vector<MatD1> dyn_fwd;
    std::vector<MatDD> M_fwd;

    ta_var.set_time(tgrid.front());
    std::copy(state0.begin(), state0.end(), ta_var.get_state_data());
    std::copy(pars_no_control.begin(), pars_no_control.end(), ta_var.get_pars_data() + C);

    unsigned successful_fwd = 0u;
    for (unsigned i = 0u; i < nseg_fwd; ++i) {
        std::copy(ic_var.begin(), ic_var.end(), ta_var.get_state_data() + D);

        const auto start = static_cast<std::size_t>(C * i);
        std::copy(controls.begin() + static_cast<std::ptrdiff_t>(start),
                  controls.begin() + static_cast<std::ptrdiff_t>(start + C), ta_var.get_pars_data());
        std::copy(pars_no_control.begin(), pars_no_control.end(), ta_var.get_pars_data() + C);

        if (!propagate_until_safe_impl(ta_var, tgrid[i + 1u], max_steps)) {
            break;
        }

        // Extract STM and control sensitivity.
        const auto &x_var_vec = ta_var.get_state();
        auto x_var
            = xt::adapt(x_var_vec.data() + static_cast<std::ptrdiff_t>(D), static_cast<std::size_t>(D * (D + C)),
                        xt::no_ownership(), stm_shape);

        MatDD M{};
        MatDC C_mat{};
        M = xt::eval(xt::view(x_var, xt::all(), xt::range(0u, D)));
        C_mat = xt::eval(xt::view(x_var, xt::all(), xt::range(D, D + C)));
        M_seg_fwd.push_back(std::move(M));
        C_seg_fwd.push_back(std::move(C_mat));

        std::vector<double> state_arr(D, 0.0);
        std::copy_n(ta_var.get_state().begin(), static_cast<std::ptrdiff_t>(D), state_arr.begin());

        std::vector<double> pars_vec;
        pars_vec.reserve(C + pars_no_control.size());
        pars_vec.insert(pars_vec.end(), controls.begin() + static_cast<std::ptrdiff_t>(start),
                        controls.begin() + static_cast<std::ptrdiff_t>(start + C));
        pars_vec.insert(pars_vec.end(), pars_no_control.begin(), pars_no_control.end());

        // Compute dynamics at this point using dyn_cfunc.
        std::vector<double> dyn_out(D, 0.0);
        dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);

        MatD1 dyn{};
        for (std::size_t j = 0u; j < D; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_fwd.push_back(std::move(dyn));

        ++successful_fwd;
    }

    // Compute M_fwd chain.
    MatDD cur = kep3::linalg::_eye<D>();
    for (auto it = M_seg_fwd.rbegin(); it != M_seg_fwd.rend(); ++it) {
        cur = _dot<D, D, D>(cur, *it);
        M_fwd.push_back(cur);
    }
    std::reverse(M_fwd.begin(), M_fwd.end());
    M_fwd.push_back(kep3::linalg::_eye<D>());

    // 1. dmc/dx0.
    if (successful_fwd == 0u) {
        const auto I = kep3::linalg::_eye<D>();
        std::copy(I.begin(), I.end(), dmc_dx0.begin());
    } else {
        std::copy(M_fwd[0].begin(), M_fwd[0].end(), dmc_dx0.begin());
    }

    // 2. dmc/dcontrols (forward).
    std::vector<double> C_fwd(D * C * nseg_fwd, 0.0);
    auto C_fwd_view = xt::adapt(C_fwd, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(C * nseg_fwd)});
    for (unsigned i = 0u; i < successful_fwd; ++i) {
        const auto prod = _dot<D, D, C>(M_fwd[i + 1u], C_seg_fwd[i]);
        xt::view(C_fwd_view, xt::all(), xt::range(C * i, C * (i + 1u))) = prod;
    }

    // 3. dmc/dtgrid (forward).
    std::vector<double> dmcdtgrid_fwd(D * (nseg + 1u), 0.0);
    auto dmcdtgrid_fwd_view
        = xt::adapt(dmcdtgrid_fwd, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(nseg + 1u)});
    if (successful_fwd > 0u) {
        auto tmp = _dot<D, D, 1u>(M_fwd[1u], dyn_fwd[0u]);
        xt::view(dmcdtgrid_fwd_view, xt::all(), 0u) = -xt::view(tmp, xt::all(), 0u);

        tmp = _dot<D, D, 1u>(M_fwd[successful_fwd], dyn_fwd.back());
        xt::view(dmcdtgrid_fwd_view, xt::all(), successful_fwd) = xt::view(tmp, xt::all(), 0u);

        for (unsigned i = 1u; i < successful_fwd; ++i) {
            auto a = _dot<D, D, 1u>(M_seg_fwd[i], dyn_fwd[i - 1u]);
            a -= dyn_fwd[i];
            tmp = _dot<D, D, 1u>(M_fwd[i + 1u], a);
            xt::view(dmcdtgrid_fwd_view, xt::all(), i) = xt::view(tmp, xt::all(), 0u);
        }
    }

    // Backward pass.
    std::vector<MatDD> M_seg_bck;
    std::vector<MatDC> C_seg_bck;
    std::vector<MatD1> dyn_bck;
    std::vector<MatDD> M_bck;

    ta_var.set_time(tgrid.back());
    std::copy(state1.begin(), state1.end(), ta_var.get_state_data());
    std::copy(pars_no_control.begin(), pars_no_control.end(), ta_var.get_pars_data() + C);

    unsigned successful_bck = 0u;
    for (unsigned i = 0u; i < nseg_bck; ++i) {
        std::copy(ic_var.begin(), ic_var.end(), ta_var.get_state_data() + D);

        const auto start = controls.size() - static_cast<std::size_t>(C * (i + 1u));
        std::copy(controls.begin() + static_cast<std::ptrdiff_t>(start),
                  controls.begin() + static_cast<std::ptrdiff_t>(start + C), ta_var.get_pars_data());
        std::copy(pars_no_control.begin(), pars_no_control.end(), ta_var.get_pars_data() + C);

        if (!propagate_until_safe_impl(ta_var, tgrid[tgrid.size() - static_cast<std::size_t>(2u + i)], max_steps)) {
            break;
        }

        const auto &x_var_vec = ta_var.get_state();
        auto x_var
            = xt::adapt(x_var_vec.data() + static_cast<std::ptrdiff_t>(D), static_cast<std::size_t>(D * (D + C)),
                        xt::no_ownership(), stm_shape);

        MatDD M{};
        MatDC C_mat{};
        M = xt::eval(xt::view(x_var, xt::all(), xt::range(0u, D)));
        C_mat = xt::eval(xt::view(x_var, xt::all(), xt::range(D, D + C)));
        M_seg_bck.push_back(std::move(M));
        C_seg_bck.push_back(std::move(C_mat));

        std::vector<double> state_arr(D, 0.0);
        std::copy_n(ta_var.get_state().begin(), static_cast<std::ptrdiff_t>(D), state_arr.begin());

        std::vector<double> pars_vec;
        pars_vec.reserve(C + pars_no_control.size());
        pars_vec.insert(pars_vec.end(), controls.begin() + static_cast<std::ptrdiff_t>(start),
                        controls.begin() + static_cast<std::ptrdiff_t>(start + C));
        pars_vec.insert(pars_vec.end(), pars_no_control.begin(), pars_no_control.end());

        std::vector<double> dyn_out(D, 0.0);
        dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);

        MatD1 dyn{};
        for (std::size_t j = 0u; j < D; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_bck.push_back(std::move(dyn));

        ++successful_bck;
    }

    cur = kep3::linalg::_eye<D>();
    for (auto it = M_seg_bck.rbegin(); it != M_seg_bck.rend(); ++it) {
        cur = _dot<D, D, D>(cur, *it);
        M_bck.push_back(cur);
    }
    std::reverse(M_bck.begin(), M_bck.end());
    M_bck.push_back(kep3::linalg::_eye<D>());

    // 4. dmc/dx1.
    if (successful_bck == 0u) {
        for (std::size_t r = 0u; r < D; ++r) {
            for (std::size_t cc = 0u; cc < D; ++cc) {
                dmc_dx1[r * D + cc] = (r == cc) ? -1.0 : 0.0;
            }
        }
    } else {
        std::copy(M_bck[0].begin(), M_bck[0].end(), dmc_dx1.begin());
        auto dmc_dx1_view = xt::adapt(dmc_dx1, mat_dd_shape);
        dmc_dx1_view *= -1.0;
    }

    // 5. dmc/dcontrols (backward).
    // Fill C_bck from right to left so that the first backward segment's gradient
    // appears in the last block (Python-style).
    std::vector<double> C_bck(D * C * nseg_bck, 0.0);
    auto C_bck_view = xt::adapt(C_bck, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(C * nseg_bck)});
    for (unsigned i = 0u; i < successful_bck; ++i) {
        const auto prod = _dot<D, D, C>(M_bck[i + 1u], C_seg_bck[i]);
        const auto block = nseg_bck - 1u - i;
        xt::view(C_bck_view, xt::all(), xt::range(C * block, C * (block + 1u))) = prod;
    }

    // 6. dmc/dtgrid (backward).
    std::vector<double> dmcdtgrid_bck(D * (nseg + 1u), 0.0);
    auto dmcdtgrid_bck_view
        = xt::adapt(dmcdtgrid_bck, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(nseg + 1u)});
    if (successful_bck > 0u) {
        auto tmp = _dot<D, D, 1u>(M_bck[1u], dyn_bck[0u]);
        xt::view(dmcdtgrid_bck_view, xt::all(), nseg) = xt::view(tmp, xt::all(), 0u);

        tmp = _dot<D, D, 1u>(M_bck[successful_bck], dyn_bck.back());
        xt::view(dmcdtgrid_bck_view, xt::all(), (nseg - successful_bck)) -= xt::view(tmp, xt::all(), 0u);

        for (unsigned i = 1u; i < successful_bck; ++i) {
            auto a = _dot<D, D, 1u>(M_seg_bck[i], dyn_bck[i - 1u]);
            a -= dyn_bck[i];
            tmp = _dot<D, D, 1u>(M_bck[i + 1u], a);
            xt::view(dmcdtgrid_bck_view, xt::all(), (nseg - i)) = -xt::view(tmp, xt::all(), 0u);
        }
    }

    // 7. Assemble dmc/dcontrols = [C_fwd, -C_bck].
    auto dmc_dcontrols_view
        = xt::adapt(dmc_dcontrols, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(C * nseg)});
    if (nseg_fwd > 0u) {
        xt::view(dmc_dcontrols_view, xt::all(), xt::range(0u, C * nseg_fwd)) = C_fwd_view;
    }
    if (nseg_bck > 0u) {
        xt::view(dmc_dcontrols_view, xt::all(), xt::range(C * nseg_fwd, C * nseg)) = -C_bck_view;
    }

    // 8. Assemble dmc/dtgrid.
    auto dmc_dtgrid_view = xt::adapt(dmc_dtgrid, std::array<std::size_t, 2u>{D, static_cast<std::size_t>(nseg + 1u)});
    dmc_dtgrid_view = dmcdtgrid_fwd_view + dmcdtgrid_bck_view;

    return {dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid};
}

bool propagate_until_safe_impl(heyoka::taylor_adaptive<double> &ta, double t, const std::optional<unsigned> &max_steps)
{
    const auto prev_time = ta.get_time();
    const auto prev_state = ta.get_state();

    try {
        if (max_steps) {
            ta.propagate_until(t, heyoka::kw::max_steps = *max_steps);
        } else {
            ta.propagate_until(t);
        }
    } catch (...) {
        ta.set_time(prev_time);
        std::copy(prev_state.begin(), prev_state.end(), ta.get_state_data());
        return false;
    }

    return true;
}

} // namespace

zoh::zoh(const std::vector<double> &state0, const std::vector<double> &controls, const std::vector<double> &state1,
         const std::vector<double> &tgrid, double cut,
         const std::pair<heyoka::taylor_adaptive<double>, std::optional<heyoka::taylor_adaptive<double>>> &tas,
         std::optional<unsigned> max_steps, unsigned dim_dynamics, unsigned dim_controls)
    : m_state0(state0), m_controls(controls), m_state1(state1), m_tgrid(tgrid), m_cut(cut), m_max_steps(max_steps),
      m_dim_dynamics(dim_dynamics), m_dim_controls(dim_controls), m_ta(tas.first), m_ta_var(tas.second)
{
    update_nseg();
    update_ic_var();
    update_pars_no_control();

    const auto &sys = m_ta.get_sys();
    std::vector<heyoka::expression> dyn, vars;
    for (const auto &pair : sys) {
        vars.push_back(pair.first);
        dyn.push_back(pair.second);
    }
    m_dyn_cfunc = heyoka::cfunc<double>(dyn, vars, heyoka::kw::compact_mode = true);

    sanity_checks();
}

void zoh::set_state0(const std::vector<double> &state0)
{
    m_state0 = state0;
    sanity_checks();
}

void zoh::set_state1(const std::vector<double> &state1)
{
    m_state1 = state1;
    sanity_checks();
}

void zoh::set_controls(const std::vector<double> &controls)
{
    m_controls = controls;
    update_nseg();
    sanity_checks();
}

void zoh::set_tgrid(const std::vector<double> &tgrid)
{
    m_tgrid = tgrid;
    sanity_checks();
}

void zoh::set_cut(double cut)
{
    m_cut = cut;
    update_nseg();
    sanity_checks();
}

void zoh::set_max_steps(std::optional<unsigned> max_steps)
{
    m_max_steps = max_steps;
}

void zoh::set(const std::vector<double> &state0, const std::vector<double> &controls, const std::vector<double> &state1,
              const std::vector<double> &tgrid, double cut, std::optional<unsigned> max_steps)
{
    m_state0 = state0;
    m_controls = controls;
    m_state1 = state1;
    m_tgrid = tgrid;
    m_cut = cut;
    m_max_steps = max_steps;
    update_nseg();
    sanity_checks();
}

const std::vector<double> &zoh::get_state0() const
{
    return m_state0;
}

const std::vector<double> &zoh::get_state1() const
{
    return m_state1;
}

const std::vector<double> &zoh::get_controls() const
{
    return m_controls;
}

const std::vector<double> &zoh::get_tgrid() const
{
    return m_tgrid;
}

double zoh::get_cut() const
{
    return m_cut;
}

unsigned zoh::get_dim_dynamics() const
{
    return m_dim_dynamics;
}

unsigned zoh::get_dim_controls() const
{
    return m_dim_controls;
}

std::optional<unsigned> zoh::get_max_steps() const
{
    return m_max_steps;
}

const heyoka::taylor_adaptive<double> &zoh::get_ta() const
{
    return m_ta;
}

const std::optional<heyoka::taylor_adaptive<double>> &zoh::get_ta_var() const
{
    return m_ta_var;
}

bool zoh::has_ta_var() const
{
    return m_ta_var.has_value();
}

unsigned zoh::get_nseg() const
{
    return m_nseg;
}

unsigned zoh::get_nseg_fwd() const
{
    return m_nseg_fwd;
}

unsigned zoh::get_nseg_bck() const
{
    return m_nseg_bck;
}

std::vector<double> zoh::compute_mismatch_constraints() const
{
    auto &ta = m_ta;

    // Forward propagation.
    ta.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        const auto start = static_cast<std::size_t>(m_dim_controls * i);
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + m_dim_controls), ta.get_pars_data());
        if (!propagate_until_safe_impl(ta, m_tgrid[i + 1u], m_max_steps)) {
            break;
        }
    }

    std::vector<double> state_fwd(m_dim_dynamics, 0.0);
    std::copy(ta.get_state().begin(), ta.get_state().begin() + static_cast<std::ptrdiff_t>(m_dim_dynamics),
              state_fwd.begin());

    // Backward propagation.
    ta.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        const auto start = m_controls.size() - static_cast<std::size_t>(m_dim_controls * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + m_dim_controls), ta.get_pars_data());
        if (!propagate_until_safe_impl(ta, m_tgrid[m_tgrid.size() - static_cast<std::size_t>(2u + i)], m_max_steps)) {
            break;
        }
    }

    std::vector<double> state_bck(m_dim_dynamics, 0.0);
    std::copy(ta.get_state().begin(), ta.get_state().begin() + static_cast<std::ptrdiff_t>(m_dim_dynamics),
              state_bck.begin());

    std::vector<double> retval(m_dim_dynamics, 0.0);
    for (unsigned i = 0u; i < m_dim_dynamics; ++i) {
        retval[i] = state_fwd[i] - state_bck[i];
    }
    return retval;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
zoh::compute_mc_grad() const
{
    using kep3::linalg::mat61;
    using kep3::linalg::mat62;
    using kep3::linalg::mat66;
    using kep3::linalg::mat71;
    using kep3::linalg::mat74;
    using kep3::linalg::mat77;

    if (!m_ta_var) {
        throw std::logic_error("zoh::compute_mc_grad() requires a variational integrator (ta_var)");
    }

    const auto d = m_dim_dynamics;
    const auto c = m_dim_controls;

    if (d == 7u && c == 4u) {
        if (m_ta_var->get_dim() != 7u + 7u * 7u + 7u * 4u) {
            throw std::logic_error("zoh::compute_mc_grad() requires ta_var with compatible variational state dimension");
        }
        return compute_mc_grad_fixed_impl<7u, 4u, mat77, mat74, mat71>(
            *m_ta_var, m_state0, m_controls, m_state1, m_tgrid, m_ic_var, m_pars_no_control, m_max_steps,
            m_dyn_cfunc, m_nseg, m_nseg_fwd, m_nseg_bck);
    }

    if (d == 6u && c == 2u) {
        if (m_ta_var->get_dim() != 6u + 6u * 6u + 6u * 2u) {
            throw std::logic_error("zoh::compute_mc_grad() requires ta_var with compatible variational state dimension");
        }
        return compute_mc_grad_fixed_impl<6u, 2u, mat66, mat62, mat61>(
            *m_ta_var, m_state0, m_controls, m_state1, m_tgrid, m_ic_var, m_pars_no_control, m_max_steps,
            m_dyn_cfunc, m_nseg, m_nseg_fwd, m_nseg_bck);
    }

    throw std::logic_error(fmt::format(
        "zoh::compute_mc_grad() not implemented for dim_dynamics={}, dim_controls={}", d, c));
}

std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, bool>
zoh::get_state_info(unsigned N) const
{
    if (N == 0u) {
        throw std::logic_error("zoh::get_state_info() requires N >= 1");
    }

    bool success = true;
    auto &ta = m_ta;

    std::vector<std::vector<std::vector<double>>> state_fwd;
    ta.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        const auto start = static_cast<std::size_t>(m_dim_controls * i);
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

        std::vector<std::vector<double>> segment_states;
        const double t0 = m_tgrid[i];
        const double t1 = m_tgrid[i + 1u];

        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? t0 : (t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(N - 1u));
            if (!propagate_until_safe_impl(ta, t, m_max_steps)) {
                success = false;
                break;
            }
            std::vector<double> state(m_dim_dynamics, 0.0);
            std::copy_n(ta.get_state().begin(), static_cast<std::ptrdiff_t>(m_dim_dynamics), state.begin());
            segment_states.push_back(std::move(state));
        }

        if (!segment_states.empty()) {
            state_fwd.push_back(std::move(segment_states));
        }
        if (!success) {
            break;
        }
    }

    std::vector<std::vector<std::vector<double>>> state_bck;
    ta.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        const auto start = m_controls.size() - static_cast<std::size_t>(m_dim_controls * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

        std::vector<std::vector<double>> segment_states;
        const double t0 = m_tgrid[m_tgrid.size() - 1u - i];
        const double t1 = m_tgrid[m_tgrid.size() - 2u - i];

        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? t0 : (t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(N - 1u));
            if (!propagate_until_safe_impl(ta, t, m_max_steps)) {
                success = false;
                break;
            }
            std::vector<double> state(m_dim_dynamics, 0.0);
            std::copy_n(ta.get_state().begin(), static_cast<std::ptrdiff_t>(m_dim_dynamics), state.begin());
            segment_states.push_back(std::move(state));
        }

        if (!segment_states.empty()) {
            state_bck.push_back(std::move(segment_states));
        }
        if (!success) {
            break;
        }
    }

    return {state_fwd, state_bck, success};
}

void zoh::update_nseg()
{
    m_nseg = static_cast<unsigned>(m_controls.size() / m_dim_controls);
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}

void zoh::update_ic_var()
{
    m_ic_var.assign(m_dim_dynamics * (m_dim_dynamics + m_dim_controls), 0.0);
    for (unsigned i = 0u; i < m_dim_dynamics; ++i) {
        m_ic_var[i * (m_dim_dynamics + m_dim_controls) + i] = 1.0;
    }
}

void zoh::update_pars_no_control()
{
    const auto &pars = m_ta.get_pars();
    m_pars_no_control.assign(pars.begin() + static_cast<std::ptrdiff_t>(m_dim_controls), pars.end());
}

void zoh::sanity_checks() const
{
    if (m_dim_dynamics == 0u || m_dim_controls == 0u) {
        throw std::logic_error("dim_dynamics and dim_controls must be positive.");
    }

    if (m_cut < 0. || m_cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh leg must be in the [0, 1] interval.");
    }

    if (m_state0.size() != m_dim_dynamics || m_state1.size() != m_dim_dynamics) {
        throw std::logic_error("state0/state1 sizes must match dim_dynamics.");
    }

    if (m_ta.get_dim() != m_dim_dynamics) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a Taylor adaptive integrator state "
                                           "dimension of {}, while {} is required.",
                                           m_ta.get_dim(), m_dim_dynamics));
    }

    if (m_ta.get_pars().size() < m_dim_controls) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a Taylor adaptive integrator "
                                           "parameters dimension of {}, while >= {} is required.",
                                           m_ta.get_pars().size(), m_dim_controls));
    }

    if ((m_controls.size() % m_dim_controls) != 0u) {
        throw std::logic_error("In a zoh leg, controls size must be a multiple of dim_controls.");
    }

    if (m_tgrid.size() != m_nseg + 1u) {
        throw std::logic_error(
            "The tgrid and controls have incompatible sizes. They must be nseg + 1 and dim_controls * nseg.");
    }

    if (m_ta_var) {
        const auto expected = m_dim_dynamics + m_dim_dynamics * m_dim_dynamics + m_dim_dynamics * m_dim_controls;
        if (m_ta_var->get_dim() != expected) {
            throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a variational Taylor adaptive "
                                               "integrator state dimension of {}, while {} is required.",
                                               m_ta_var->get_dim(), expected));
        }
        if (m_ta_var->get_pars().size() != m_ta.get_pars().size()) {
            throw std::logic_error("The variational and nominal Taylor adaptive integrators for a zoh leg must expose "
                                   "the same number of parameters.");
        }
    }
}

std::ostream &operator<<(std::ostream &s, const zoh &leg)
{
    s << fmt::format("Dynamics dimension: {}\n", leg.get_dim_dynamics());
    s << fmt::format("Control dimension: {}\n", leg.get_dim_controls());
    s << fmt::format("Number of segments: {}\n", leg.get_nseg());
    s << fmt::format("Number of fwd segments: {}\n", leg.get_nseg_fwd());
    s << fmt::format("Number of bck segments: {}\n", leg.get_nseg_bck());
    s << fmt::format("Cut parameter: {}\n", leg.get_cut());
    s << fmt::format("Variational integrator available: {}\n", leg.has_ta_var());
    if (leg.get_max_steps()) {
        s << fmt::format("Maximum propagation steps: {}\n\n", *leg.get_max_steps());
    } else {
        s << "Maximum propagation steps: none\n\n";
    }
    s << fmt::format("Initial state: {}\n", leg.get_state0());
    s << fmt::format("Final state: {}\n", leg.get_state1());
    s << fmt::format("Time grid: {}\n", leg.get_tgrid());
    s << fmt::format("Control values: {}\n\n", leg.get_controls());
    s << fmt::format("Mismatch constraints: {}\n", leg.compute_mismatch_constraints());
    return s;
}

} // namespace kep3::leg
