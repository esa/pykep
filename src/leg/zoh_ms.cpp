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

#include <heyoka/expression.hpp>
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
    if (m_dim_dynamics == 0u || m_dim_controls == 0u) {
        throw std::logic_error("dim_dynamics and dim_controls must be positive.");
    }
    if (m_cut < 0. || m_cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh_ms leg must be in the [0, 1] interval.");
    }

    update_nseg();
    sanity_checks();
    initialize_ic_var();
    update_pars_no_control();

    const auto &sys = m_ta.get_sys();
    std::vector<heyoka::expression> dyn, vars;
    for (const auto &equation : sys) {
        vars.push_back(equation.first);
        dyn.push_back(equation.second);
    }
    m_dyn_cfunc = heyoka::cfunc<double>(dyn, vars, heyoka::kw::compact_mode = true);
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

const std::tuple<zoh_ms::sparsity_pattern, zoh_ms::sparsity_pattern, zoh_ms::sparsity_pattern> &
zoh_ms::defects_grad_sparsity() const
{
    return m_defects_grad_sparsity;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> zoh_ms::compute_defects_grad() const
{
    if (!m_ta_var) {
        throw std::logic_error("zoh_ms::compute_defects_grad() requires a variational integrator (ta_var)");
    }

    auto &ta_var = *m_ta_var;
    const auto dimension = static_cast<std::size_t>(m_dim_dynamics);
    const auto control_dimension = static_cast<std::size_t>(m_dim_controls);
    const auto sensitivity_dimension = dimension + control_dimension;
    std::vector<double> grad_states;
    std::vector<double> grad_controls;
    std::vector<double> grad_tgrid;
    grad_states.reserve(static_cast<std::size_t>(m_nseg) * (dimension * dimension + dimension));
    grad_controls.reserve(static_cast<std::size_t>(m_nseg) * dimension * control_dimension);
    grad_tgrid.reserve(static_cast<std::size_t>(m_nseg) * dimension * 2u);

    for (unsigned segment = 0u; segment < m_nseg; ++segment) {
        const bool forward = segment < m_nseg_fwd;
        const auto start_node = forward ? segment : segment + 1u;
        const auto end_node = forward ? segment + 1u : segment;
        const auto state_offset = static_cast<std::size_t>(m_dim_dynamics * start_node);
        const auto control_offset = static_cast<std::size_t>(m_dim_controls * segment);

        std::vector<double> initial_state(dimension);
        std::copy_n(m_states.begin() + static_cast<std::ptrdiff_t>(state_offset),
                    static_cast<std::ptrdiff_t>(dimension), initial_state.begin());
        std::vector<double> parameters;
        parameters.reserve(m_ta.get_pars().size());
        parameters.insert(parameters.end(), m_controls.begin() + static_cast<std::ptrdiff_t>(control_offset),
                          m_controls.begin() + static_cast<std::ptrdiff_t>(control_offset + control_dimension));
        parameters.insert(parameters.end(), m_pars_no_control.begin(), m_pars_no_control.end());

        // Each defect uses sensitivities initialized at its own starting node.
        ta_var.set_time(m_tgrid[start_node]);
        std::copy(initial_state.begin(), initial_state.end(), ta_var.get_state_data());
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + dimension);
        std::copy(parameters.begin(), parameters.end(), ta_var.get_pars_data());
        const bool success = propagate_until_safe_impl(ta_var, m_tgrid[end_node], m_max_steps);

        // A failed propagation restores [x0, I, 0], which is also the required fallback.
        std::vector<double> state_transition(dimension * dimension, 0.0);
        std::vector<double> control_sensitivity(dimension * control_dimension, 0.0);
        if (success) {
            const auto &variational_state = ta_var.get_state();
            for (std::size_t row = 0u; row < dimension; ++row) {
                for (std::size_t column = 0u; column < dimension; ++column) {
                    state_transition[row * dimension + column]
                        = variational_state[dimension + row * sensitivity_dimension + column];
                }
                for (std::size_t column = 0u; column < control_dimension; ++column) {
                    control_sensitivity[row * control_dimension + column]
                        = variational_state[dimension + row * sensitivity_dimension + dimension + column];
                }
            }
        } else {
            for (std::size_t component = 0u; component < dimension; ++component) {
                state_transition[component * dimension + component] = 1.0;
            }
        }

        std::vector<double> dynamics_start(dimension, 0.0);
        std::vector<double> dynamics_end(dimension, 0.0);
        if (success) {
            std::vector<double> final_state(dimension);
            std::copy_n(ta_var.get_state().begin(), static_cast<std::ptrdiff_t>(dimension), final_state.begin());
            m_dyn_cfunc(dynamics_start, initial_state, heyoka::kw::pars = parameters);
            m_dyn_cfunc(dynamics_end, final_state, heyoka::kw::pars = parameters);
        }

        // Moving the segment start changes the flow by -M f(x0); moving its end gives f(x1).
        std::vector<double> time_derivative_start(dimension, 0.0);
        if (success) {
            for (std::size_t row = 0u; row < dimension; ++row) {
                for (std::size_t column = 0u; column < dimension; ++column) {
                    time_derivative_start[row] -= state_transition[row * dimension + column] * dynamics_start[column];
                }
            }
        }

        // Append each row in the same order as the cached sparsity coordinates.
        for (std::size_t row = 0u; row < dimension; ++row) {
            if (forward) {
                for (std::size_t column = 0u; column < dimension; ++column) {
                    grad_states.push_back(state_transition[row * dimension + column]);
                }
                grad_states.push_back(-1.0);
            } else {
                grad_states.push_back(1.0);
                for (std::size_t column = 0u; column < dimension; ++column) {
                    grad_states.push_back(-state_transition[row * dimension + column]);
                }
            }

            const double control_sign = forward ? 1.0 : -1.0;
            for (std::size_t column = 0u; column < control_dimension; ++column) {
                grad_controls.push_back(control_sign * control_sensitivity[row * control_dimension + column]);
            }

            grad_tgrid.push_back(forward ? time_derivative_start[row] : -dynamics_end[row]);
            grad_tgrid.push_back(forward ? dynamics_end[row] : -time_derivative_start[row]);
        }
    }

    return {std::move(grad_states), std::move(grad_controls), std::move(grad_tgrid)};
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
    update_sparsity();
}

void zoh_ms::initialize_ic_var()
{
    const auto dimension = static_cast<std::size_t>(m_dim_dynamics);
    const auto sensitivity_dimension = dimension + static_cast<std::size_t>(m_dim_controls);
    m_ic_var.assign(dimension * sensitivity_dimension, 0.0);
    for (std::size_t component = 0u; component < dimension; ++component) {
        m_ic_var[component * sensitivity_dimension + component] = 1.0;
    }
}

void zoh_ms::update_sparsity()
{
    auto &[sp_states, sp_controls, sp_tgrid] = m_defects_grad_sparsity;
    sp_states.clear();
    sp_controls.clear();
    sp_tgrid.clear();

    const auto d = static_cast<std::size_t>(m_dim_dynamics);
    const auto c = static_cast<std::size_t>(m_dim_controls);
    const auto nseg = static_cast<std::size_t>(m_nseg);
    sp_states.reserve(nseg * (d * d + d));
    sp_controls.reserve(nseg * d * c);
    sp_tgrid.reserve(nseg * d * 2u);

    for (unsigned i = 0u; i < m_nseg; ++i) {
        for (unsigned k = 0u; k < m_dim_dynamics; ++k) {
            const auto row = static_cast<std::size_t>(m_dim_dynamics * i + k);
            if (i < m_nseg_fwd) {
                for (unsigned j = 0u; j < m_dim_dynamics; ++j) {
                    sp_states.emplace_back(row, static_cast<std::size_t>(m_dim_dynamics * i + j));
                }
                sp_states.emplace_back(row, static_cast<std::size_t>(m_dim_dynamics * (i + 1u) + k));
            } else {
                sp_states.emplace_back(row, static_cast<std::size_t>(m_dim_dynamics * i + k));
                for (unsigned j = 0u; j < m_dim_dynamics; ++j) {
                    sp_states.emplace_back(row, static_cast<std::size_t>(m_dim_dynamics * (i + 1u) + j));
                }
            }

            for (unsigned j = 0u; j < m_dim_controls; ++j) {
                sp_controls.emplace_back(row, static_cast<std::size_t>(m_dim_controls * i + j));
            }
            sp_tgrid.emplace_back(row, static_cast<std::size_t>(i));
            sp_tgrid.emplace_back(row, static_cast<std::size_t>(i + 1u));
        }
    }
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