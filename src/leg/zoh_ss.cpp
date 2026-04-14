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
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/leg/zoh_ss.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg
{

namespace
{

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

zoh_ss::zoh_ss(const std::array<double, 6> &state0, const std::vector<double> &controls,
               const std::array<double, 6> &state1, const std::vector<double> &tgrid, double cut,
               const std::pair<heyoka::taylor_adaptive<double>, std::optional<heyoka::taylor_adaptive<double>>> &tas,
               std::optional<unsigned> max_steps)
    : m_state0(state0), m_controls(controls), m_state1(state1), m_tgrid(tgrid), m_cut(cut), m_max_steps(max_steps),
      m_ta(tas.first), m_ta_var(tas.second)
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

void zoh_ss::set_state0(const std::array<double, 6> &state0)
{
    m_state0 = state0;
}

void zoh_ss::set_state1(const std::array<double, 6> &state1)
{
    m_state1 = state1;
}

void zoh_ss::set_controls(const std::vector<double> &controls)
{
    m_controls = controls;
    update_nseg();
    sanity_checks();
}

void zoh_ss::set_tgrid(const std::vector<double> &tgrid)
{
    m_tgrid = tgrid;
    sanity_checks();
}

void zoh_ss::set_cut(double cut)
{
    m_cut = cut;
    update_nseg();
    sanity_checks();
}

void zoh_ss::set_max_steps(std::optional<unsigned> max_steps)
{
    m_max_steps = max_steps;
}

void zoh_ss::set(const std::array<double, 6> &state0, const std::vector<double> &controls,
                 const std::array<double, 6> &state1, const std::vector<double> &tgrid, double cut,
                 std::optional<unsigned> max_steps)
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

const std::array<double, 6> &zoh_ss::get_state0() const
{
    return m_state0;
}

const std::array<double, 6> &zoh_ss::get_state1() const
{
    return m_state1;
}

const std::vector<double> &zoh_ss::get_controls() const
{
    return m_controls;
}

const std::vector<double> &zoh_ss::get_tgrid() const
{
    return m_tgrid;
}

double zoh_ss::get_cut() const
{
    return m_cut;
}

std::optional<unsigned> zoh_ss::get_max_steps() const
{
    return m_max_steps;
}

const heyoka::taylor_adaptive<double> &zoh_ss::get_ta() const
{
    return m_ta;
}

const std::optional<heyoka::taylor_adaptive<double>> &zoh_ss::get_ta_var() const
{
    return m_ta_var;
}

bool zoh_ss::has_ta_var() const
{
    return m_ta_var.has_value();
}

unsigned zoh_ss::get_nseg() const
{
    return m_nseg;
}

unsigned zoh_ss::get_nseg_fwd() const
{
    return m_nseg_fwd;
}

unsigned zoh_ss::get_nseg_bck() const
{
    return m_nseg_bck;
}

std::array<double, 6> zoh_ss::compute_mismatch_constraints() const
{
    auto &ta = m_ta;

    // Forward segments.
    ta.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

    for (decltype(m_nseg_fwd) i = 0u; i < m_nseg_fwd; ++i) {
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i + 2u), ta.get_pars_data());
        if (!propagate_until_safe_impl(ta, m_tgrid[static_cast<std::size_t>(i) + 1u], m_max_steps)) {
            break;
        }
    }

    std::array<double, 6> state_fwd{};
    std::copy(ta.get_state().begin(), ta.get_state().begin() + 6, state_fwd.begin());

    // Backward segments.
    ta.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

    for (decltype(m_nseg_bck) i = 0u; i < m_nseg_bck; ++i) {
        const auto control_idx = m_controls.size() - static_cast<std::size_t>(2u * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx + 2u), ta.get_pars_data());
        if (!propagate_until_safe_impl(ta, m_tgrid[m_tgrid.size() - static_cast<std::size_t>(2u + i)], m_max_steps)) {
            break;
        }
    }

    std::array<double, 6> state_bck{};
    std::copy(ta.get_state().begin(), ta.get_state().begin() + 6, state_bck.begin());

    std::array<double, 6> retval{};
    for (std::size_t i = 0u; i < retval.size(); ++i) {
        retval[i] = state_fwd[i] - state_bck[i];
    }
    return retval;
}

std::tuple<std::array<double, 36>, std::array<double, 36>, std::vector<double>, std::vector<double>>
zoh_ss::compute_mc_grad() const
{
    using kep3::linalg::_dot;
    using kep3::linalg::mat61;
    using kep3::linalg::mat62;
    using kep3::linalg::mat66;
    using kep3::linalg::_eye;
    using xt::adapt;

    if (!m_ta_var) {
        throw std::logic_error("zoh_ss::compute_mc_grad() requires a variational integrator (ta_var)");
    }
    if (m_ta_var->get_dim() != 54u) {
        throw std::logic_error("zoh_ss::compute_mc_grad() requires ta_var with state dim 54");
    }

    std::array<double, 36> dmc_dx0{};
    std::array<double, 36> dmc_dx1{};
    std::vector<double> dmc_dcontrols(6 * 2 * m_nseg, 0.0);
    std::vector<double> dmc_dtgrid(6 * (m_nseg + 1), 0.0);

    std::vector<mat66> M_seg_fwd;
    std::vector<mat62> C_seg_fwd;
    std::vector<mat61> dyn_fwd;
    std::vector<mat66> M_fwd;

    auto &ta_var = *m_ta_var;
    ta_var.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta_var.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 2);

    unsigned successful_fwd = 0u;
    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + 6);
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i + 2u), ta_var.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 2);

        if (!propagate_until_safe_impl(ta_var, m_tgrid[i + 1], m_max_steps)) {
            break;
        }

        auto x_var = adapt(ta_var.get_state(), {54u});
        mat66 M{};
        mat62 C{};
        for (unsigned r = 0; r < 6; ++r) {
            for (unsigned c = 0; c < 6; ++c) {
                M(r, c) = x_var(6 + r * 8 + c);
            }
            for (unsigned c = 0; c < 2; ++c) {
                C(r, c) = x_var(6 + r * 8 + 6 + c);
            }
        }
        M_seg_fwd.push_back(M);
        C_seg_fwd.push_back(C);

        mat61 dyn{};
        std::array<double, 6> state_arr{};
        std::copy_n(ta_var.get_state().begin(), 6, state_arr.begin());
        std::vector<double> pars_vec;
        pars_vec.reserve(2 + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i),
                        m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i + 2u));
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());
        std::array<double, 6> dyn_out{};
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        for (unsigned j = 0; j < 6; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_fwd.push_back(dyn);

        ++successful_fwd;
    }

    mat66 cur = _eye<6>();
    for (auto it = M_seg_fwd.rbegin(); it != M_seg_fwd.rend(); ++it) {
        cur = _dot<6, 6, 6>(cur, *it);
        M_fwd.push_back(cur);
    }
    std::reverse(M_fwd.begin(), M_fwd.end());
    M_fwd.push_back(_eye<6>());

    std::copy(M_fwd[0].data(), M_fwd[0].data() + 36, dmc_dx0.begin());
    if (successful_fwd == 0u) {
        dmc_dx0.fill(0.);
        for (unsigned i = 0u; i < 6u; ++i) {
            dmc_dx0[6u * i + i] = 1.;
        }
    }

    xt::xarray<double> C_fwd = xt::zeros<double>({6u, 2u * m_nseg_fwd});
    for (unsigned i = 0u; i < successful_fwd; ++i) {
        auto prod = _dot<6, 6, 2>(M_fwd[i + 1], C_seg_fwd[i]);
        xt::view(C_fwd, xt::all(), xt::range(2 * i, 2 * i + 2)) = prod;
    }

    xt::xarray<double> dmcdtgrid_fwd = xt::zeros<double>({6u, m_nseg + 1, 1u});
    if (successful_fwd > 0u) {
        xt::view(dmcdtgrid_fwd, xt::all(), 0, xt::all()) = xt::eval(-_dot<6, 6, 1>(M_fwd[1], dyn_fwd[0]));
        xt::view(dmcdtgrid_fwd, xt::all(), successful_fwd, xt::all())
            = xt::eval(_dot<6, 6, 1>(M_fwd[successful_fwd], dyn_fwd.back()));

        for (unsigned i = 1u; i < successful_fwd; ++i) {
            auto tmp = _dot<6, 6, 1>(M_seg_fwd[i], dyn_fwd[i - 1]) - dyn_fwd[i];
            xt::view(dmcdtgrid_fwd, xt::all(), i, xt::all()) = xt::eval(_dot<6, 6, 1>(M_fwd[i + 1], tmp));
        }
    }

    std::vector<mat66> M_seg_bck;
    std::vector<mat62> C_seg_bck;
    std::vector<mat61> dyn_bck;
    std::vector<mat66> M_bck;

    ta_var.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta_var.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 2);

    unsigned successful_bck = 0u;
    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + 6);
        const auto control_idx = m_controls.size() - static_cast<std::size_t>(2u * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx + 2u), ta_var.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 2);

        if (!propagate_until_safe_impl(ta_var, m_tgrid[m_tgrid.size() - 2u - i], m_max_steps)) {
            break;
        }

        auto x_var = adapt(ta_var.get_state(), {54u});
        mat66 M{};
        mat62 C{};
        for (unsigned r = 0; r < 6; ++r) {
            for (unsigned c = 0; c < 6; ++c) {
                M(r, c) = x_var(6 + r * 8 + c);
            }
            for (unsigned c = 0; c < 2; ++c) {
                C(r, c) = x_var(6 + r * 8 + 6 + c);
            }
        }
        M_seg_bck.push_back(M);
        C_seg_bck.push_back(C);

        mat61 dyn{};
        std::array<double, 6> state_arr{};
        std::copy_n(ta_var.get_state().begin(), 6, state_arr.begin());
        std::vector<double> pars_vec;
        pars_vec.reserve(2 + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx),
                        m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx + 2u));
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());
        std::array<double, 6> dyn_out{};
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        for (unsigned j = 0; j < 6; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_bck.push_back(dyn);

        ++successful_bck;
    }

    cur = _eye<6>();
    for (auto it = M_seg_bck.rbegin(); it != M_seg_bck.rend(); ++it) {
        cur = _dot<6, 6, 6>(cur, *it);
        M_bck.push_back(cur);
    }
    std::reverse(M_bck.begin(), M_bck.end());
    M_bck.push_back(_eye<6>());

    for (unsigned i = 0; i < 36; ++i) {
        dmc_dx1[i] = -M_bck[0].data()[i];
    }
    if (successful_bck == 0u) {
        dmc_dx1.fill(0.);
        for (unsigned i = 0u; i < 6u; ++i) {
            dmc_dx1[6u * i + i] = -1.;
        }
    }

    xt::xarray<double> C_bck = xt::zeros<double>({6u, 2u * m_nseg_bck});
    for (unsigned i = 0u; i < successful_bck; ++i) {
        auto prod = _dot<6, 6, 2>(M_bck[i + 1], C_seg_bck[i]);
        const unsigned block = m_nseg_bck - 1u - i;
        xt::view(C_bck, xt::all(), xt::range(2 * block, 2 * block + 2)) = prod;
    }

    xt::xarray<double> dmcdtgrid_bck = xt::zeros<double>({6u, m_nseg + 1, 1u});
    if (successful_bck > 0u) {
        xt::view(dmcdtgrid_bck, xt::all(), m_nseg, xt::all()) = xt::eval(_dot<6, 6, 1>(M_bck[1], dyn_bck[0]));
        xt::view(dmcdtgrid_bck, xt::all(), m_nseg - successful_bck, xt::all())
            -= xt::eval(_dot<6, 6, 1>(M_bck[successful_bck], dyn_bck.back()));

        for (unsigned i = 1u; i < successful_bck; ++i) {
            auto tmp = _dot<6, 6, 1>(M_seg_bck[i], dyn_bck[i - 1]) - dyn_bck[i];
            xt::view(dmcdtgrid_bck, xt::all(), m_nseg - i, xt::all()) = xt::eval(-_dot<6, 6, 1>(M_bck[i + 1], tmp));
        }
    }

    auto x_dmc_dcontrols = adapt(dmc_dcontrols, {6u, 2u * m_nseg});
    if (m_nseg_fwd > 0u) {
        xt::view(x_dmc_dcontrols, xt::all(), xt::range(0, 2 * m_nseg_fwd)) = C_fwd;
    }
    if (m_nseg_bck > 0u) {
        xt::view(x_dmc_dcontrols, xt::all(), xt::range(2 * m_nseg_fwd, 2 * m_nseg)) = -C_bck;
    }

    auto x_dmc_dtgrid = adapt(dmc_dtgrid, {6u, m_nseg + 1, 1u});
    x_dmc_dtgrid = dmcdtgrid_fwd + dmcdtgrid_bck;

    return {dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid};
}

std::tuple<std::vector<std::vector<std::array<double, 6>>>, std::vector<std::vector<std::array<double, 6>>>, bool>
zoh_ss::get_state_info(unsigned N) const
{
    if (N == 0u) {
        throw std::logic_error("zoh_ss::get_state_info() requires N >= 1");
    }

    bool success = true;

    std::vector<std::vector<std::array<double, 6>>> state_fwd;
    auto &ta = m_ta;
    ta.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(2u * i + 2u), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

        std::vector<std::array<double, 6>> segment_states;
        const double t0 = m_tgrid[i];
        const double t1 = m_tgrid[i + 1u];

        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? t0 : (t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(N - 1u));
            if (!propagate_until_safe_impl(ta, t, m_max_steps)) {
                success = false;
                break;
            }
            std::array<double, 6> state{};
            std::copy(ta.get_state().begin(), ta.get_state().begin() + 6, state.begin());
            segment_states.push_back(state);
        }

        if (!success) {
            break;
        }
        state_fwd.push_back(segment_states);
    }

    std::vector<std::vector<std::array<double, 6>>> state_bck;
    ta.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        const auto control_idx = m_controls.size() - static_cast<std::size_t>(2u * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx + 2u), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 2);

        std::vector<std::array<double, 6>> segment_states;
        const double t0 = m_tgrid[m_tgrid.size() - 1u - i];
        const double t1 = m_tgrid[m_tgrid.size() - 2u - i];

        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? t0 : (t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(N - 1u));
            if (!propagate_until_safe_impl(ta, t, m_max_steps)) {
                success = false;
                break;
            }
            std::array<double, 6> state{};
            std::copy(ta.get_state().begin(), ta.get_state().begin() + 6, state.begin());
            segment_states.push_back(state);
        }

        if (!success) {
            break;
        }
        state_bck.push_back(segment_states);
    }

    return {state_fwd, state_bck, success};
}

void zoh_ss::update_nseg()
{
    m_nseg = static_cast<unsigned>(m_controls.size() / 2u);
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}

void zoh_ss::update_ic_var()
{
    m_ic_var.fill(0.);
    for (std::size_t i = 0u; i < 6u; ++i) {
        m_ic_var[i * 8u + i] = 1.;
    }
}

void zoh_ss::update_pars_no_control()
{
    const auto &pars = m_ta.get_pars();
    m_pars_no_control.assign(pars.begin() + 2, pars.end());
}

void zoh_ss::sanity_checks() const
{
    if (m_cut < 0. || m_cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh_ss leg must be in the [0, 1] interval.");
    }

    if (m_ta.get_dim() != 6u) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh_ss leg with a Taylor adaptive integrator state "
                                           "dimension of {}, while 6 is required.",
                                           m_ta.get_dim()));
    }

    if (m_ta.get_pars().size() <= 2u) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh_ss leg with a Taylor adaptive integrator "
                                           "parameters dimension of {}, while >2 is required.",
                                           m_ta.get_pars().size()));
    }

    if ((m_controls.size() % 2u) != 0u) {
        throw std::logic_error("In a zoh_ss leg the controls size must be a multiple of 2: [alpha, beta] * nseg.");
    }

    if (m_tgrid.size() != m_nseg + 1u) {
        throw std::logic_error("The tgrid and controls have incompatible sizes. They must be nseg + 1 and 2 * nseg.");
    }

    if (m_ta_var) {
        if (m_ta_var->get_dim() != 54u) {
            throw std::logic_error(fmt::format("Attempting to construct a zoh_ss leg with a variational Taylor adaptive "
                                               "integrator state dimension of {}, while 54 is required.",
                                               m_ta_var->get_dim()));
        }
        if (m_ta_var->get_pars().size() != m_ta.get_pars().size()) {
            throw std::logic_error("The variational and nominal Taylor adaptive integrators for a zoh_ss leg must expose "
                                   "the same number of parameters.");
        }
    }
}

std::ostream &operator<<(std::ostream &s, const zoh_ss &leg)
{
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
