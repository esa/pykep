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

#include <kep3/leg/zoh_ms.hpp>

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

zoh_ms::zoh_ms(const std::vector<double> &states, const std::vector<double> &controls,
               const std::vector<double> &tgrid, double cut,
               const std::pair<heyoka::taylor_adaptive<double>, std::optional<heyoka::taylor_adaptive<double>>> &tas,
               std::optional<unsigned> max_steps, unsigned dim_dynamics, unsigned dim_controls)
    : m_states(states), m_controls(controls), m_tgrid(tgrid), m_cut(cut), m_max_steps(max_steps),
      m_dim_dynamics(dim_dynamics), m_dim_controls(dim_controls), m_ta(tas.first), m_ta_var(tas.second)
{
    update_nseg();
    update_pars_no_control();

    sanity_checks();
}

void zoh_ms::set_states(const std::vector<double> &states)
{
    if (states.size() != m_states.size()) {
        throw std::logic_error("states size must match the existing states size.");
    }
    m_states = states;
}

void zoh_ms::set_controls(const std::vector<double> &controls)
{
    if (controls.size() != m_controls.size()) {
        throw std::logic_error("controls size must match the existing controls size.");
    }
    m_controls = controls;
}

void zoh_ms::set_tgrid(const std::vector<double> &tgrid)
{
    if (tgrid.size() != m_tgrid.size()) {
        throw std::logic_error("tgrid size must match the existing tgrid size.");
    }
    m_tgrid = tgrid;
}

void zoh_ms::set_cut(double cut)
{
    if (cut < 0. || cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh_ms leg must be in the [0, 1] interval.");
    }
    m_cut = cut;
    update_nseg();
}

void zoh_ms::set_max_steps(std::optional<unsigned> max_steps)
{
    m_max_steps = max_steps;
}

void zoh_ms::set(const std::vector<double> &states, const std::vector<double> &controls,
                 const std::vector<double> &tgrid, std::optional<double> cut, std::optional<unsigned> max_steps)
{
    if (m_dim_dynamics == 0u || m_dim_controls == 0u) {
        throw std::logic_error("dim_dynamics and dim_controls must be positive.");
    }
    if (cut && (*cut < 0. || *cut > 1.)) {
        throw std::logic_error("The cut parameter of a zoh_ms leg must be in the [0, 1] interval.");
    }
    if ((controls.size() % m_dim_controls) != 0u || tgrid.size() != controls.size() / m_dim_controls + 1u
        || states.size() != (controls.size() / m_dim_controls + 1u) * m_dim_dynamics) {
        throw std::logic_error(
            "The states, tgrid and controls have incompatible sizes. They must be (nseg + 1) * dim_dynamics, "
            "nseg + 1 and dim_controls * nseg.");
    }

    m_states = states;
    m_controls = controls;
    m_tgrid = tgrid;
    if (cut) {
        m_cut = *cut;
    }
    m_max_steps = max_steps;
    update_nseg();
}

const std::vector<double> &zoh_ms::get_states() const
{
    return m_states;
}

const std::vector<double> &zoh_ms::get_controls() const
{
    return m_controls;
}

const std::vector<double> &zoh_ms::get_tgrid() const
{
    return m_tgrid;
}

double zoh_ms::get_cut() const
{
    return m_cut;
}

unsigned zoh_ms::get_dim_dynamics() const
{
    return m_dim_dynamics;
}

unsigned zoh_ms::get_dim_controls() const
{
    return m_dim_controls;
}

std::optional<unsigned> zoh_ms::get_max_steps() const
{
    return m_max_steps;
}

const heyoka::taylor_adaptive<double> &zoh_ms::get_ta() const
{
    return m_ta;
}

const std::optional<heyoka::taylor_adaptive<double>> &zoh_ms::get_ta_var() const
{
    return m_ta_var;
}

bool zoh_ms::has_ta_var() const
{
    return m_ta_var.has_value();
}

unsigned zoh_ms::get_nseg() const
{
    return m_nseg;
}

unsigned zoh_ms::get_nseg_fwd() const
{
    return m_nseg_fwd;
}

unsigned zoh_ms::get_nseg_bck() const
{
    return m_nseg_bck;
}

std::vector<double> zoh_ms::compute_defects() const
{
    auto &ta = m_ta;
    std::vector<double> defects(m_nseg * m_dim_dynamics, 0.0);

    for (unsigned i = 0u; i < m_nseg; ++i) {
        const bool forward = i < m_nseg_fwd;
        const auto start_node = forward ? i : i + 1u;
        const auto end_node = forward ? i + 1u : i;
        const auto state_start = static_cast<std::size_t>(m_dim_dynamics * start_node);
        const auto state_end = static_cast<std::size_t>(m_dim_dynamics * end_node);
        const auto control_start = static_cast<std::size_t>(m_dim_controls * i);

        ta.set_time(m_tgrid[start_node]);
        std::copy(m_states.begin() + static_cast<std::ptrdiff_t>(state_start),
                  m_states.begin() + static_cast<std::ptrdiff_t>(state_start + m_dim_dynamics), ta.get_state_data());
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

        propagate_until_safe_impl(ta, m_tgrid[end_node], m_max_steps);

        const auto defect_start = static_cast<std::size_t>(m_dim_dynamics * i);
        for (unsigned j = 0u; j < m_dim_dynamics; ++j) {
            const auto propagated = ta.get_state()[j];
            const auto target = m_states[state_end + j];
            defects[defect_start + j] = forward ? propagated - target : target - propagated;
        }
    }

    return defects;
}

void zoh_ms::set_initial_guess(bool ballistic)
{
    auto &ta = m_ta;
    const auto &controls = ballistic ? std::vector<double>(m_controls.size(), 0.0) : m_controls;

    ta.set_time(m_tgrid.front());
    std::copy(m_states.begin(), m_states.begin() + static_cast<std::ptrdiff_t>(m_dim_dynamics), ta.get_state_data());

    const auto nseg_fwd = std::min(m_nseg_fwd, m_nseg > 0u ? m_nseg - 1u : 0u);
    for (unsigned i = 0u; i < nseg_fwd; ++i) {
        const auto control_start = static_cast<std::size_t>(m_dim_controls * i);
        std::copy(controls.begin() + static_cast<std::ptrdiff_t>(control_start),
                  controls.begin() + static_cast<std::ptrdiff_t>(control_start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);
        propagate_until_safe_impl(ta, m_tgrid[i + 1u], m_max_steps);

        const auto state_start = static_cast<std::size_t>(m_dim_dynamics * (i + 1u));
        std::copy_n(ta.get_state().begin(), static_cast<std::ptrdiff_t>(m_dim_dynamics),
                    m_states.begin() + static_cast<std::ptrdiff_t>(state_start));
    }

    ta.set_time(m_tgrid.back());
    std::copy(m_states.end() - static_cast<std::ptrdiff_t>(m_dim_dynamics), m_states.end(), ta.get_state_data());

    for (unsigned i = m_nseg; i > m_nseg_fwd + 1u; --i) {
        const auto segment = i - 1u;
        const auto control_start = static_cast<std::size_t>(m_dim_controls * segment);
        std::copy(controls.begin() + static_cast<std::ptrdiff_t>(control_start),
                  controls.begin() + static_cast<std::ptrdiff_t>(control_start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);
        propagate_until_safe_impl(ta, m_tgrid[segment], m_max_steps);

        const auto state_start = static_cast<std::size_t>(m_dim_dynamics * segment);
        std::copy_n(ta.get_state().begin(), static_cast<std::ptrdiff_t>(m_dim_dynamics),
                    m_states.begin() + static_cast<std::ptrdiff_t>(state_start));
    }
}

std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, bool>
zoh_ms::get_state_info(unsigned N) const
{
    if (N == 0u) {
        throw std::logic_error("zoh_ms::get_state_info() requires N >= 1");
    }

    bool success = true;
    auto &ta = m_ta;
    std::vector<std::vector<std::vector<double>>> state_fwd;
    std::vector<std::vector<std::vector<double>>> state_bck;

    for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
        const auto state_start = static_cast<std::size_t>(m_dim_dynamics * i);
        const auto control_start = static_cast<std::size_t>(m_dim_controls * i);
        ta.set_time(m_tgrid[i]);
        std::copy(m_states.begin() + static_cast<std::ptrdiff_t>(state_start),
                  m_states.begin() + static_cast<std::ptrdiff_t>(state_start + m_dim_dynamics), ta.get_state_data());
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

        std::vector<std::vector<double>> segment_states;
        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? m_tgrid[i]
                                       : (m_tgrid[i] + (m_tgrid[i + 1u] - m_tgrid[i]) * static_cast<double>(k)
                                          / static_cast<double>(N - 1u));
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

    for (unsigned i = 0u; i < m_nseg_bck; ++i) {
        const auto segment = m_nseg - 1u - i;
        const auto state_start = static_cast<std::size_t>(m_dim_dynamics * (segment + 1u));
        const auto control_start = static_cast<std::size_t>(m_dim_controls * segment);
        ta.set_time(m_tgrid[segment + 1u]);
        std::copy(m_states.begin() + static_cast<std::ptrdiff_t>(state_start),
                  m_states.begin() + static_cast<std::ptrdiff_t>(state_start + m_dim_dynamics), ta.get_state_data());
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_start),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_start + m_dim_controls), ta.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + m_dim_controls);

        std::vector<std::vector<double>> segment_states;
        for (unsigned k = 0u; k < N; ++k) {
            const double t = (N == 1u) ? m_tgrid[segment + 1u]
                                       : (m_tgrid[segment + 1u]
                                          + (m_tgrid[segment] - m_tgrid[segment + 1u]) * static_cast<double>(k)
                                          / static_cast<double>(N - 1u));
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

void zoh_ms::update_nseg()
{
    m_nseg = static_cast<unsigned>(m_controls.size() / m_dim_controls);
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}

void zoh_ms::update_pars_no_control()
{
    const auto &pars = m_ta.get_pars();
    m_pars_no_control.assign(pars.begin() + static_cast<std::ptrdiff_t>(m_dim_controls), pars.end());
}

void zoh_ms::sanity_checks() const
{
    if (m_dim_dynamics == 0u || m_dim_controls == 0u) {
        throw std::logic_error("dim_dynamics and dim_controls must be positive.");
    }

    if (m_cut < 0. || m_cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh_ms leg must be in the [0, 1] interval.");
    }

    if (m_states.size() != (m_nseg + 1u) * m_dim_dynamics) {
        throw std::logic_error("states size must be (nseg + 1) * dim_dynamics.");
    }

    if (m_ta.get_dim() != m_dim_dynamics) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh_ms leg with a Taylor adaptive integrator "
                                           "state dimension of {}, while {} is required.",
                                           m_ta.get_dim(), m_dim_dynamics));
    }

    if (m_ta.get_pars().size() < m_dim_controls) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh_ms leg with a Taylor adaptive integrator "
                                           "parameters dimension of {}, while >= {} is required.",
                                           m_ta.get_pars().size(), m_dim_controls));
    }

    if ((m_controls.size() % m_dim_controls) != 0u) {
        throw std::logic_error("In a zoh_ms leg, controls size must be a multiple of dim_controls.");
    }

    if (m_tgrid.size() != m_nseg + 1u) {
        throw std::logic_error(
            "The tgrid and controls have incompatible sizes. They must be nseg + 1 and dim_controls * nseg.");
    }

    if (m_ta_var) {
        const auto expected = m_dim_dynamics + m_dim_dynamics * m_dim_dynamics + m_dim_dynamics * m_dim_controls;
        if (m_ta_var->get_dim() != expected) {
            throw std::logic_error(fmt::format("Attempting to construct a zoh_ms leg with a variational Taylor "
                                               "adaptive integrator state dimension of {}, while {} is required.",
                                               m_ta_var->get_dim(), expected));
        }
        if (m_ta_var->get_pars().size() != m_ta.get_pars().size()) {
            throw std::logic_error("The variational and nominal Taylor adaptive integrators for a zoh_ms leg must "
                                   "expose the same number of parameters.");
        }
    }
}

std::ostream &operator<<(std::ostream &s, const zoh_ms &leg)
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
    s << fmt::format("States: {}\n", leg.get_states());
    s << fmt::format("Time grid: {}\n", leg.get_tgrid());
    s << fmt::format("Control values: {}\n\n", leg.get_controls());
    s << fmt::format("Defects: {}\n", leg.compute_defects());
    return s;
}

} // namespace kep3::leg