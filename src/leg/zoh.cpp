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
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <heyoka/kw.hpp>
#include <heyoka/taylor.hpp>

#include <kep3/leg/zoh.hpp>

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

std::vector<double> make_identity(unsigned n)
{
    std::vector<double> out(n * n, 0.0);
    for (unsigned i = 0u; i < n; ++i) {
        out[i * n + i] = 1.0;
    }
    return out;
}

std::vector<double> matmul(const std::vector<double> &A, const std::vector<double> &B, unsigned m, unsigned n, unsigned p)
{
    std::vector<double> out(m * p, 0.0);
    for (unsigned i = 0u; i < m; ++i) {
        for (unsigned k = 0u; k < n; ++k) {
            const auto aik = A[i * n + k];
            if (aik == 0.0) {
                continue;
            }
            for (unsigned j = 0u; j < p; ++j) {
                out[i * p + j] += aik * B[k * p + j];
            }
        }
    }
    return out;
}

std::vector<double> matvec(const std::vector<double> &A, const std::vector<double> &x, unsigned m, unsigned n)
{
    std::vector<double> out(m, 0.0);
    for (unsigned i = 0u; i < m; ++i) {
        for (unsigned k = 0u; k < n; ++k) {
            out[i] += A[i * n + k] * x[k];
        }
    }
    return out;
}

void axpy_inplace(std::vector<double> &dst, const std::vector<double> &src, double alpha)
{
    for (std::size_t i = 0u; i < dst.size(); ++i) {
        dst[i] += alpha * src[i];
    }
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

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> zoh::compute_mc_grad() const
{
    if (!m_ta_var) {
        throw std::logic_error("zoh::compute_mc_grad() requires a variational integrator (ta_var)");
    }

    const auto d = m_dim_dynamics;
    const auto c = m_dim_controls;
    const auto stm_cols = d + c;

    if (m_ta_var->get_dim() != d + d * d + d * c) {
        throw std::logic_error("zoh::compute_mc_grad() requires ta_var with compatible variational state dimension");
    }

    std::vector<double> dmc_dx0(d * d, 0.0);
    std::vector<double> dmc_dx1(d * d, 0.0);
    std::vector<double> dmc_dcontrols(d * c * m_nseg, 0.0);
    std::vector<double> dmc_dtgrid(d * (m_nseg + 1u), 0.0);

    auto &ta_var = *m_ta_var;

    // Forward pass.
    std::vector<std::vector<double>> M_seg_fwd;
    std::vector<std::vector<double>> C_seg_fwd;
    std::vector<std::vector<double>> dyn_fwd;
    std::vector<std::vector<double>> M_fwd;

    ta_var.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta_var.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + c);

    unsigned successful_fwd = 0u;
    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + d);

        const auto start = static_cast<std::size_t>(c * i);
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + c), ta_var.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + c);

        if (!propagate_until_safe_impl(ta_var, m_tgrid[i + 1u], m_max_steps)) {
            break;
        }

        std::vector<double> M(d * d, 0.0);
        std::vector<double> C(d * c, 0.0);
        const auto &x_var = ta_var.get_state();
        for (unsigned r = 0u; r < d; ++r) {
            for (unsigned cc = 0u; cc < d; ++cc) {
                M[r * d + cc] = x_var[d + r * stm_cols + cc];
            }
            for (unsigned cc = 0u; cc < c; ++cc) {
                C[r * c + cc] = x_var[d + r * stm_cols + d + cc];
            }
        }
        M_seg_fwd.push_back(std::move(M));
        C_seg_fwd.push_back(std::move(C));

        std::vector<double> state_arr(d, 0.0);
        std::copy_n(ta_var.get_state().begin(), static_cast<std::ptrdiff_t>(d), state_arr.begin());

        std::vector<double> pars_vec;
        pars_vec.reserve(c + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                        m_controls.begin() + static_cast<std::ptrdiff_t>(start + c));
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());

        std::vector<double> dyn_out(d, 0.0);
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        dyn_fwd.push_back(std::move(dyn_out));

        ++successful_fwd;
    }

    auto cur = make_identity(d);
    for (auto it = M_seg_fwd.rbegin(); it != M_seg_fwd.rend(); ++it) {
        cur = matmul(cur, *it, d, d, d);
        M_fwd.push_back(cur);
    }
    std::reverse(M_fwd.begin(), M_fwd.end());
    M_fwd.push_back(make_identity(d));

    if (successful_fwd == 0u) {
        dmc_dx0 = make_identity(d);
    } else {
        dmc_dx0 = M_fwd[0];
    }

    std::vector<double> C_fwd(d * c * m_nseg_fwd, 0.0);
    for (unsigned i = 0u; i < successful_fwd; ++i) {
        const auto prod = matmul(M_fwd[i + 1u], C_seg_fwd[i], d, d, c);
        for (unsigned r = 0u; r < d; ++r) {
            for (unsigned cc = 0u; cc < c; ++cc) {
                C_fwd[r * (c * m_nseg_fwd) + c * i + cc] = prod[r * c + cc];
            }
        }
    }

    std::vector<double> dmcdtgrid_fwd(d * (m_nseg + 1u), 0.0);
    if (successful_fwd > 0u) {
        auto tmp = matvec(M_fwd[1u], dyn_fwd[0u], d, d);
        for (unsigned r = 0u; r < d; ++r) {
            dmcdtgrid_fwd[r * (m_nseg + 1u)] = -tmp[r];
        }

        tmp = matvec(M_fwd[successful_fwd], dyn_fwd.back(), d, d);
        for (unsigned r = 0u; r < d; ++r) {
            dmcdtgrid_fwd[r * (m_nseg + 1u) + successful_fwd] = tmp[r];
        }

        for (unsigned i = 1u; i < successful_fwd; ++i) {
            auto a = matvec(M_seg_fwd[i], dyn_fwd[i - 1u], d, d);
            for (unsigned r = 0u; r < d; ++r) {
                a[r] -= dyn_fwd[i][r];
            }
            tmp = matvec(M_fwd[i + 1u], a, d, d);
            for (unsigned r = 0u; r < d; ++r) {
                dmcdtgrid_fwd[r * (m_nseg + 1u) + i] = tmp[r];
            }
        }
    }

    // Backward pass.
    std::vector<std::vector<double>> M_seg_bck;
    std::vector<std::vector<double>> C_seg_bck;
    std::vector<std::vector<double>> dyn_bck;
    std::vector<std::vector<double>> M_bck;

    ta_var.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta_var.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + c);

    unsigned successful_bck = 0u;
    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + d);

        const auto start = m_controls.size() - static_cast<std::size_t>(c * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(start + c), ta_var.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + c);

        if (!propagate_until_safe_impl(ta_var, m_tgrid[m_tgrid.size() - static_cast<std::size_t>(2u + i)], m_max_steps)) {
            break;
        }

        std::vector<double> M(d * d, 0.0);
        std::vector<double> C(d * c, 0.0);
        const auto &x_var = ta_var.get_state();
        for (unsigned r = 0u; r < d; ++r) {
            for (unsigned cc = 0u; cc < d; ++cc) {
                M[r * d + cc] = x_var[d + r * stm_cols + cc];
            }
            for (unsigned cc = 0u; cc < c; ++cc) {
                C[r * c + cc] = x_var[d + r * stm_cols + d + cc];
            }
        }
        M_seg_bck.push_back(std::move(M));
        C_seg_bck.push_back(std::move(C));

        std::vector<double> state_arr(d, 0.0);
        std::copy_n(ta_var.get_state().begin(), static_cast<std::ptrdiff_t>(d), state_arr.begin());

        std::vector<double> pars_vec;
        pars_vec.reserve(c + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.begin() + static_cast<std::ptrdiff_t>(start),
                        m_controls.begin() + static_cast<std::ptrdiff_t>(start + c));
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());

        std::vector<double> dyn_out(d, 0.0);
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        dyn_bck.push_back(std::move(dyn_out));

        ++successful_bck;
    }

    cur = make_identity(d);
    for (auto it = M_seg_bck.rbegin(); it != M_seg_bck.rend(); ++it) {
        cur = matmul(cur, *it, d, d, d);
        M_bck.push_back(cur);
    }
    std::reverse(M_bck.begin(), M_bck.end());
    M_bck.push_back(make_identity(d));

    if (successful_bck == 0u) {
        dmc_dx1.assign(d * d, 0.0);
        for (unsigned i = 0u; i < d; ++i) {
            dmc_dx1[i * d + i] = -1.0;
        }
    } else {
        dmc_dx1 = M_bck[0];
        for (auto &v : dmc_dx1) {
            v = -v;
        }
    }

    std::vector<double> C_bck(d * c * m_nseg_bck, 0.0);
    for (unsigned i = 0u; i < successful_bck; ++i) {
        const auto prod = matmul(M_bck[i + 1u], C_seg_bck[i], d, d, c);
        const auto block = m_nseg_bck - 1u - i;
        for (unsigned r = 0u; r < d; ++r) {
            for (unsigned cc = 0u; cc < c; ++cc) {
                C_bck[r * (c * m_nseg_bck) + c * block + cc] = prod[r * c + cc];
            }
        }
    }

    std::vector<double> dmcdtgrid_bck(d * (m_nseg + 1u), 0.0);
    if (successful_bck > 0u) {
        auto tmp = matvec(M_bck[1u], dyn_bck[0u], d, d);
        for (unsigned r = 0u; r < d; ++r) {
            dmcdtgrid_bck[r * (m_nseg + 1u) + m_nseg] = tmp[r];
        }

        tmp = matvec(M_bck[successful_bck], dyn_bck.back(), d, d);
        for (unsigned r = 0u; r < d; ++r) {
            dmcdtgrid_bck[r * (m_nseg + 1u) + (m_nseg - successful_bck)] -= tmp[r];
        }

        for (unsigned i = 1u; i < successful_bck; ++i) {
            auto a = matvec(M_seg_bck[i], dyn_bck[i - 1u], d, d);
            for (unsigned r = 0u; r < d; ++r) {
                a[r] -= dyn_bck[i][r];
            }
            tmp = matvec(M_bck[i + 1u], a, d, d);
            for (unsigned r = 0u; r < d; ++r) {
                dmcdtgrid_bck[r * (m_nseg + 1u) + (m_nseg - i)] = -tmp[r];
            }
        }
    }

    // Assemble dmc/dcontrols = [C_fwd, -C_bck].
    for (unsigned r = 0u; r < d; ++r) {
        for (unsigned j = 0u; j < c * m_nseg_fwd; ++j) {
            dmc_dcontrols[r * (c * m_nseg) + j] = C_fwd[r * (c * m_nseg_fwd) + j];
        }
        for (unsigned j = 0u; j < c * m_nseg_bck; ++j) {
            dmc_dcontrols[r * (c * m_nseg) + (c * m_nseg_fwd + j)] = -C_bck[r * (c * m_nseg_bck) + j];
        }
    }

    // Assemble dmc/dtgrid.
    for (unsigned r = 0u; r < d; ++r) {
        for (unsigned j = 0u; j < m_nseg + 1u; ++j) {
            dmc_dtgrid[r * (m_nseg + 1u) + j] = dmcdtgrid_fwd[r * (m_nseg + 1u) + j] + dmcdtgrid_bck[r * (m_nseg + 1u) + j];
        }
    }

    return {dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid};
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
        throw std::logic_error("The tgrid and controls have incompatible sizes. They must be nseg + 1 and dim_controls * nseg.");
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
