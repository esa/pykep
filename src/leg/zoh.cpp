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
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>

#include <kep3/leg/zoh.hpp>
#include <kep3/linalg.hpp>

namespace kep3::leg

{

namespace
{

void propagate_until_impl(heyoka::taylor_adaptive<double> &ta, double t, const std::optional<unsigned> &max_steps)
{
    // Keep the propagation call-site uniform while preserving the optional max_steps setting.
    if (max_steps) {
        ta.propagate_until(t, heyoka::kw::max_steps = *max_steps);
    } else {
        ta.propagate_until(t);
    }
}

} // namespace

zoh::zoh(const std::array<double, 7> &state0, const std::vector<double> &controls, const std::array<double, 7> &state1,
         const std::vector<double> &tgrid, double cut, const heyoka::taylor_adaptive<double> &ta,
         std::optional<heyoka::taylor_adaptive<double>> ta_var, std::optional<unsigned> max_steps)
    : m_state0(state0), m_controls(controls), m_state1(state1), m_tgrid(tgrid), m_cut(cut), m_max_steps(max_steps),
      m_ta(ta), m_ta_var(std::move(ta_var))
{
    // The integrators are fixed at construction: the object invariants depend on their dimensions and parameter layout.
    update_nseg();
    update_ic_var();
    update_pars_no_control();
    // Construct cfunc from the system expressions and variables of the nominal integrator
    const auto &sys = m_ta.get_sys();
    std::vector<heyoka::expression> dyn, vars;
    for (const auto &pair : sys) {
        vars.push_back(pair.first);
        dyn.push_back(pair.second);
    }
    m_dyn_cfunc = heyoka::cfunc<double>(dyn, vars, heyoka::kw::compact_mode = true);
    sanity_checks();
}

void zoh::set_state0(const std::array<double, 7> &state0)
{
    m_state0 = state0;
}

void zoh::set_state1(const std::array<double, 7> &state1)
{
    m_state1 = state1;
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

void zoh::set(const std::array<double, 7> &state0, const std::vector<double> &controls,
              const std::array<double, 7> &state1, const std::vector<double> &tgrid, double cut,
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

const std::array<double, 7> &zoh::get_state0() const
{
    return m_state0;
}

const std::array<double, 7> &zoh::get_state1() const
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

std::vector<double> zoh::compute_throttle_constraints() const
{
    std::vector<double> retval(m_nseg, 0.);
    for (decltype(m_nseg) i = 0u; i < m_nseg; ++i) {
        retval[i] = m_controls[4u * i + 1u] * m_controls[4u * i + 1u]
                    + m_controls[4u * i + 2u] * m_controls[4u * i + 2u]
                    + m_controls[4u * i + 3u] * m_controls[4u * i + 3u] - 1.;
    }
    return retval;
}

std::array<double, 7> zoh::compute_mismatch_constraints() const
{
    // Intentionally alias the data member: propagation must operate on the stored integrator, not on a copy.
    auto &ta = m_ta;

    // Forward half-leg: propagate from the initial state up to the cut point.
    ta.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 4);

    for (decltype(m_nseg_fwd) i = 0u; i < m_nseg_fwd; ++i) {
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(4u * i),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(4u * i + 4u), ta.get_pars_data());
        propagate_until_impl(ta, m_tgrid[static_cast<std::size_t>(i) + 1u], m_max_steps);
    }

    std::array<double, 7> state_fwd{};
    std::copy(ta.get_state().begin(), ta.get_state().begin() + 7, state_fwd.begin());

    // Backward half-leg: restart from the final state and traverse the tail segments in reverse order.
    ta.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
    std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta.get_pars_data() + 4);

    for (decltype(m_nseg_bck) i = 0u; i < m_nseg_bck; ++i) {
        const auto control_idx = m_controls.size() - static_cast<std::size_t>(4u * (i + 1u));
        std::copy(m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx),
                  m_controls.begin() + static_cast<std::ptrdiff_t>(control_idx + 4u), ta.get_pars_data());
        propagate_until_impl(ta, m_tgrid[m_tgrid.size() - static_cast<std::size_t>(2u + i)], m_max_steps);
    }

    std::array<double, 7> state_bck{};
    std::copy(ta.get_state().begin(), ta.get_state().begin() + 7, state_bck.begin());

    // The leg is feasible when the forward and backward midpoint states coincide.
    std::array<double, 7> retval{};
    for (std::size_t i = 0u; i < retval.size(); ++i) {
        retval[i] = state_fwd[i] - state_bck[i];
    }
    return retval;
}

std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>, std::vector<double>>
zoh::compute_mc_grad() const
{
    using kep3::linalg::_dot;
    using kep3::linalg::mat71;
    using kep3::linalg::mat74;
    using kep3::linalg::mat77;
    using xt::adapt;
    using xt::all;
    using xt::range;
    using xt::view;

    if (!m_ta_var) {
        throw std::logic_error("zoh::compute_mc_grad() requires a variational integrator (ta_var)");
    }
    if (m_ta_var->get_dim() != 84u) {
        throw std::logic_error("zoh::compute_mc_grad() requires ta_var with state dim 84");
    }

    // Prepare output containers
    std::array<double, 49> dmc_dx0{};
    std::array<double, 49> dmc_dx1{};
    std::vector<double> dmc_dcontrols(7 * 4 * m_nseg, 0.0);
    std::vector<double> dmc_dtgrid(7 * (m_nseg + 1), 0.0);

    // Views for state, controls, tgrid, ic_var, pars_no_control
    auto x_state0 = adapt(m_state0);
    auto x_state1 = adapt(m_state1);
    auto x_controls = adapt(m_controls, {m_nseg * 4});
    auto x_tgrid = adapt(m_tgrid, {m_nseg + 1});
    auto x_ic_var = adapt(m_ic_var, {77});
    auto x_pars_no_control = adapt(m_pars_no_control);

    // Forward pass: propagate and collect STMs and control sensitivities
    std::vector<mat77> M_seg_fwd;
    std::vector<mat74> C_seg_fwd;
    std::vector<mat71> dyn_fwd;
    std::vector<mat77> M_fwd;

    auto &ta_var = *m_ta_var;
    ta_var.set_time(m_tgrid.front());
    std::copy(m_state0.begin(), m_state0.end(), ta_var.get_state_data());

    for (unsigned i = 0; i < m_nseg_fwd; ++i) {
        // Set variational IC
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + 7);
        // Set controls
        std::copy(m_controls.begin() + 4 * i, m_controls.begin() + 4 * i + 4, ta_var.get_pars_data());
        // Set non-control parameters
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 4);
        propagate_until_impl(ta_var, m_tgrid[i + 1], m_max_steps);
        // Extract STM and control sensitivity
        auto x_var = adapt(ta_var.get_state(), {84});
        mat77 M{};
        mat74 C{};
        for (unsigned r = 0; r < 7; ++r) {
            for (unsigned c = 0; c < 7; ++c) {
                M(r, c) = x_var(7 + r * 11 + c);
            }
            for (unsigned c = 0; c < 4; ++c) {
                C(r, c) = x_var(7 + r * 11 + 7 + c);
            }
        }
        M_seg_fwd.push_back(M);
        C_seg_fwd.push_back(C);
        // Compute dynamics at this point using m_dyn_cfunc
        mat71 dyn{};
        std::array<double, 7> state_arr{};
        std::copy_n(ta_var.get_state().begin(), 7, state_arr.begin());
        std::vector<double> pars_vec;
        pars_vec.reserve(4 + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.begin() + 4 * i, m_controls.begin() + 4 * i + 4);
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());
        std::array<double, 7> dyn_out{};
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        for (unsigned j = 0; j < 7; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_fwd.push_back(dyn);
    }
    // Compute M_fwd chain
    mat77 cur = xt::eye<double>(7);
    for (auto it = M_seg_fwd.rbegin(); it != M_seg_fwd.rend(); ++it) {
        cur = _dot<7, 7, 7>(cur, *it);
        M_fwd.push_back(cur);
    }
    std::reverse(M_fwd.begin(), M_fwd.end());
    M_fwd.push_back(xt::eye<double>(7));

    // 1. dmc/dx0
    std::copy(M_fwd[0].data(), M_fwd[0].data() + 49, dmc_dx0.begin());

    // 2. dmc/dcontrols (forward)
    xt::xarray<double> C_fwd = xt::zeros<double>({7u, 4u * m_nseg_fwd});
    for (unsigned i = 0; i < m_nseg_fwd; ++i) {
        auto prod = _dot<7, 7, 4>(M_fwd[i + 1], C_seg_fwd[i]);
        xt::view(C_fwd, xt::all(), xt::range(4 * i, 4 * i + 4)) = prod;
    }

    // 3. dmc/dtgrid (forward)
    xt::xarray<double> dmcdtgrid_fwd = xt::zeros<double>({7u, m_nseg + 1, 1u});
    if (m_nseg_fwd > 0) {
        xt::view(dmcdtgrid_fwd, xt::all(), 0, xt::all()) = xt::eval(-_dot<7, 7, 1>(M_fwd[1], dyn_fwd[0]));
        xt::view(dmcdtgrid_fwd, xt::all(), m_nseg_fwd, xt::all())
            = xt::eval(_dot<7, 7, 1>(M_fwd.back(), dyn_fwd.back()));
        for (unsigned i = 1; i < m_nseg_fwd; ++i) {
            auto tmp = _dot<7, 7, 1>(M_seg_fwd[i], dyn_fwd[i - 1]) - dyn_fwd[i];
            xt::view(dmcdtgrid_fwd, xt::all(), i, xt::all()) = xt::eval(_dot<7, 7, 1>(M_fwd[i + 1], tmp));
        }
    }

    // Backward pass: propagate and collect STMs and control sensitivities
    std::vector<mat77> M_seg_bck;
    std::vector<mat74> C_seg_bck;
    std::vector<mat71> dyn_bck;
    std::vector<mat77> M_bck;

    ta_var.set_time(m_tgrid.back());
    std::copy(m_state1.begin(), m_state1.end(), ta_var.get_state_data());

    for (unsigned i = 0; i < m_nseg_bck; ++i) {
        std::copy(m_ic_var.begin(), m_ic_var.end(), ta_var.get_state_data() + 7);
        // Controls in reverse order
        unsigned control_idx = m_controls.size() - 4 * (i + 1);
        std::copy(m_controls.begin() + control_idx, m_controls.begin() + control_idx + 4, ta_var.get_pars_data());
        std::copy(m_pars_no_control.begin(), m_pars_no_control.end(), ta_var.get_pars_data() + 4);
        propagate_until_impl(ta_var, m_tgrid[m_tgrid.size() - 2 - i], m_max_steps);
        auto x_var = adapt(ta_var.get_state(), {84});
        mat77 M{};
        mat74 C{};
        for (unsigned r = 0; r < 7; ++r) {
            for (unsigned c = 0; c < 7; ++c) {
                M(r, c) = x_var(7 + r * 11 + c);
            }
            for (unsigned c = 0; c < 4; ++c) {
                C(r, c) = x_var(7 + r * 11 + 7 + c);
            }
        }
        M_seg_bck.push_back(M);
        C_seg_bck.push_back(C);
        mat71 dyn{};
        std::array<double, 7> state_arr{};
        std::copy_n(ta_var.get_state().begin(), 7, state_arr.begin());
        std::vector<double> pars_vec;
        pars_vec.reserve(4 + m_pars_no_control.size());
        pars_vec.insert(pars_vec.end(), m_controls.end() - 4 * (i + 1), m_controls.end() - 4 * i);
        pars_vec.insert(pars_vec.end(), m_pars_no_control.begin(), m_pars_no_control.end());
        std::array<double, 7> dyn_out{};
        m_dyn_cfunc(dyn_out, state_arr, heyoka::kw::pars = pars_vec);
        for (unsigned j = 0; j < 7; ++j) {
            dyn(j, 0) = dyn_out[j];
        }
        dyn_bck.push_back(dyn);
    }
    // Compute M_bck chain
    cur = xt::eye<double>(7);
    for (auto it = M_seg_bck.rbegin(); it != M_seg_bck.rend(); ++it) {
        cur = _dot<7, 7, 7>(cur, *it);
        M_bck.push_back(cur);
    }
    std::reverse(M_bck.begin(), M_bck.end());
    M_bck.push_back(xt::eye<double>(7));

    // 4. dmc/dx1
    for (unsigned i = 0; i < 49; ++i)
        dmc_dx1[i] = -M_bck[0].data()[i];

    // 5. dmc/dcontrols (backward)
    // Fill C_bck from right to left so that the first backward segment's gradient appears in the last block
    // (Python-style)
    xt::xarray<double> C_bck = xt::zeros<double>({7u, 4u * m_nseg_bck});
    for (unsigned i = 0; i < m_nseg_bck; ++i) {
        auto prod = _dot<7, 7, 4>(M_bck[i + 1], C_seg_bck[i]);
        unsigned block = m_nseg_bck - 1 - i;
        xt::view(C_bck, xt::all(), xt::range(4 * block, 4 * block + 4)) = prod;
    }

    // 6. dmc/dtgrid (backward)
    xt::xarray<double> dmcdtgrid_bck = xt::zeros<double>({7u, m_nseg + 1, 1u});
    if (m_nseg_bck > 0) {
        xt::view(dmcdtgrid_bck, xt::all(), m_nseg, xt::all()) = xt::eval(_dot<7, 7, 1>(M_bck[1], dyn_bck[0]));
        xt::view(dmcdtgrid_bck, xt::all(), m_nseg_fwd, xt::all())
            -= xt::eval(_dot<7, 7, 1>(M_bck.back(), dyn_bck.back()));
        for (unsigned i = 1; i < m_nseg_bck; ++i) {
            auto tmp = _dot<7, 7, 1>(M_seg_bck[i], dyn_bck[i - 1]) - dyn_bck[i];
            xt::view(dmcdtgrid_bck, xt::all(), m_nseg - i, xt::all()) = xt::eval(-_dot<7, 7, 1>(M_bck[i + 1], tmp));
        }
    }

    // 7. Assemble dmc/dcontrols (7 x 4*nseg): [C_fwd, -C_bck]
    auto x_dmc_dcontrols = adapt(dmc_dcontrols, {7u, 4u * m_nseg});
    if (m_nseg_fwd > 0)
        xt::view(x_dmc_dcontrols, xt::all(), xt::range(0, 4 * m_nseg_fwd)) = C_fwd;
    if (m_nseg_bck > 0)
        xt::view(x_dmc_dcontrols, xt::all(), xt::range(4 * m_nseg_fwd, 4 * m_nseg)) = -C_bck;

    // 8. Assemble dmc/dtgrid (7 x (nseg+1,1)): sum of fwd and bck, then flatten last dim
    auto x_dmc_dtgrid = adapt(dmc_dtgrid, {7u, m_nseg + 1, 1u});
    x_dmc_dtgrid = dmcdtgrid_fwd + dmcdtgrid_bck;

    return {dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid};
}

std::vector<double> zoh::compute_tc_grad() const
{
    // Returns a flattened row-major matrix of shape nseg x (4 * nseg).
    std::vector<double> retval(m_nseg * 4 * m_nseg, 0.0);
    for (unsigned i = 0; i < m_nseg; ++i) {
        retval[i * 4 * m_nseg + 4 * i + 1] = 2.0 * m_controls[4 * i + 1];
        retval[i * 4 * m_nseg + 4 * i + 2] = 2.0 * m_controls[4 * i + 2];
        retval[i * 4 * m_nseg + 4 * i + 3] = 2.0 * m_controls[4 * i + 3];
    }
    return retval;
}

void zoh::update_nseg()
{
    // The control layout is fixed to [T, ix, iy, iz] repeated once per segment.
    m_nseg = static_cast<unsigned>(m_controls.size() / 4u);
    m_nseg_fwd = static_cast<unsigned>(static_cast<double>(m_nseg) * m_cut);
    m_nseg_bck = m_nseg - m_nseg_fwd;
}

void zoh::update_ic_var()
{
    // This is the flattened [I_7 | 0_(7x4)] initial conditions expected by the first-order variational propagators.
    m_ic_var.fill(0.);
    for (std::size_t i = 0u; i < 7u; ++i) {
        m_ic_var[i * 11u + i] = 1.;
    }
}

void zoh::update_pars_no_control()
{
    auto const &pars = m_ta.get_pars();
    // The first four parameters are reserved for the segment controls; everything else is copied through unchanged.
    m_pars_no_control.assign(pars.begin() + 4, pars.end());
}

void zoh::sanity_checks() const
{
    if (m_cut < 0. || m_cut > 1.) {
        throw std::logic_error("The cut parameter of a zoh leg must be in the [0, 1] interval.");
    }

    if (m_ta.get_dim() != 7u) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a Taylor adaptive integrator state "
                                           "dimension of {}, while 7 is required.",
                                           m_ta.get_dim()));
    }

    if (m_ta.get_pars().size() <= 4u) {
        throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a Taylor adaptive integrator "
                                           "parameters dimension of {}, while >4 is required.",
                                           m_ta.get_pars().size()));
    }

    // The public API accepts a flat control vector; segment count is inferred entirely from its size.
    if ((m_controls.size() % 4u) != 0u) {
        throw std::logic_error("In a zoh leg the controls size must be a multiple of 4: [T, ix, iy, iz] * nseg.");
    }

    // The grid defines one boundary per segment plus the initial time.
    if (m_tgrid.size() != m_nseg + 1u) {
        throw std::logic_error("The tgrid and controls have incompatible sizes. They must be nseg + 1 and 4 * nseg.");
    }

    if (m_ta_var) {
        // Gradients assume a 7 + 7x7 + 7x4 variational layout.
        if (m_ta_var->get_dim() != 84u) {
            throw std::logic_error(fmt::format("Attempting to construct a zoh leg with a variational Taylor adaptive "
                                               "integrator state dimension of {}, while 84 is required.",
                                               m_ta_var->get_dim()));
        }
        if (m_ta_var->get_pars().size() != m_ta.get_pars().size()) {
            throw std::logic_error("The variational and nominal Taylor adaptive integrators for a zoh leg must expose "
                                   "the same number of parameters.");
        }
    }
}

std::pair<std::vector<std::vector<std::array<double, 7>>>, std::vector<std::vector<std::array<double, 7>>>>
zoh::get_state_info(unsigned N) const
{
	// Forward segments
	std::vector<std::vector<std::array<double, 7>>> state_fwd;
	auto &ta = m_ta;
	ta.set_time(m_tgrid.front());
	std::copy(m_state0.begin(), m_state0.end(), ta.get_state_data());
	for (unsigned i = 0u; i < m_nseg_fwd; ++i) {
		std::copy(m_controls.begin() + 4 * i, m_controls.begin() + 4 * i + 4, ta.get_pars_data());
		std::vector<std::array<double, 7>> segment_states;
		// Uniform grid for this segment
		double t0 = m_tgrid[i];
		double t1 = m_tgrid[i + 1];
		for (unsigned k = 0u; k < N; ++k) {
			double t = t0 + (t1 - t0) * k / (N - 1);
			propagate_until_impl(ta, t, m_max_steps);
			std::array<double, 7> state;
			std::copy(ta.get_state().begin(), ta.get_state().begin() + 7, state.begin());
			segment_states.push_back(state);
		}
		state_fwd.push_back(segment_states);
	}

	// Backward segments
	std::vector<std::vector<std::array<double, 7>>> state_bck;
	ta.set_time(m_tgrid.back());
	std::copy(m_state1.begin(), m_state1.end(), ta.get_state_data());
	for (unsigned i = 0u; i < m_nseg_bck; ++i) {
		unsigned control_idx = m_controls.size() - 4 * (i + 1);
		std::copy(m_controls.begin() + control_idx, m_controls.begin() + control_idx + 4, ta.get_pars_data());
		std::vector<std::array<double, 7>> segment_states;
		double t0 = m_tgrid[m_tgrid.size() - 1 - i];
		double t1 = m_tgrid[m_tgrid.size() - 2 - i];
		for (unsigned k = 0u; k < N; ++k) {
			double t = t0 + (t1 - t0) * k / (N - 1);
			propagate_until_impl(ta, t, m_max_steps);
			std::array<double, 7> state;
			std::copy(ta.get_state().begin(), ta.get_state().begin() + 7, state.begin());
			segment_states.push_back(state);
		}
		state_bck.push_back(segment_states);
	}
	return {state_fwd, state_bck};
}

std::ostream &operator<<(std::ostream &s, const zoh &leg)
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
    s << fmt::format("Throttle constraints: {}\n\n", leg.compute_throttle_constraints());
    return s;
}



} // namespace kep3::leg
