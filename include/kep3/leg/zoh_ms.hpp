// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_ZOH_MS_H
#define kep3_LEG_ZOH_MS_H

#include <cstddef>
#include <optional>
#include <ostream>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/ostream.h>

#include <heyoka/taylor.hpp>

#include <kep3/detail/visibility.hpp>

namespace kep3::leg
{
/**
 * @brief Generic zero-order-hold trajectory leg with forward/backward multiple shooting.
 *
 * The flat vector states contains nseg + 1 mesh nodes of dimension dim_dynamics.
 * The flat vector controls contains nseg piecewise-constant controls of dimension dim_controls,
 * and tgrid contains the nseg + 1 node times. Segments may have non-uniform durations.
 * Each segment is propagated independently from its starting mesh node, rather than from the
 * result of the preceding segment. With nseg_fwd = floor(nseg * cut), the first nseg_fwd
 * segments are propagated forward and the remaining ones backward.
 *
 * The dynamics are provided by a user-supplied pair of compatible Taylor-adaptive integrators.
 * The nominal integrator must have state dimension ``dim_dynamics`` and at least
 * ``dim_controls`` parameters, with the first ``dim_controls`` parameters representing the
 * segment controls. When provided, the variational integrator must implement the same dynamics
 * and expose the corresponding first-order variations with respect to state and controls.
 *
 * A transfer is feasible when all nseg * dim_dynamics defects vanish. For invertible flows,
 * the feasible trajectories do not depend on cut, but the off-feasibility defects and their
 * Jacobian conditioning do. Any
 * additional constraints on controls (e.g. throttle constraints) are intentionally left to the
 * caller (typically a UDP).
 *
 */

class kep3_DLL_PUBLIC zoh_ms
{

public:
    using sparsity_pattern = std::vector<std::pair<std::size_t, std::size_t>>;

    /// Construct an empty object for deserialization; not a propagatable leg.
    zoh_ms() = default;

    /**
     * @brief Construct a multiple-shooting leg, copying the supplied data and integrators.
     * @param states Flat mesh nodes, of length (nseg + 1) * dim_dynamics.
     * @param controls Flat segment controls, of length nseg * dim_controls.
     * @param tgrid Node times, of length nseg + 1, in the integrators' time units.
     * @param cut Fraction of segments propagated forward, in [0, 1].
     * @param tas Nominal integrator and optional first-order variational integrator with respect
     * to states and controls. Control parameters precede non-control parameters.
     * @param max_steps Maximum integration steps per segment; nullopt means no limit.
     * @param dim_dynamics State dimension. Default is 7.
     * @param dim_controls Control dimension. Default is 4.
    * @note Sparsity patterns are computed after validation and segment-count initialization.
     */
    zoh_ms(const std::vector<double> &states, const std::vector<double> &controls,
        const std::vector<double> &tgrid, double cut,
        const std::pair<heyoka::taylor_adaptive<double>, std::optional<heyoka::taylor_adaptive<double>>> &tas,
        std::optional<unsigned> max_steps = std::nullopt, unsigned dim_dynamics = 7u, unsigned dim_controls = 4u);

    /**
     * @brief Replace mesh data or propagation settings; integrators and dimensions remain unchanged.
     * @note Setters that change nseg or nseg_fwd refresh the cached sparsity patterns after
     * validation and segment-count updates. Numerical-value changes alone preserve the cache.
     */
    void set_states(const std::vector<double> &states);
    void set_controls(const std::vector<double> &controls);
    void set_tgrid(const std::vector<double> &tgrid);
    void set_cut(double cut);
    void set_max_steps(std::optional<unsigned> max_steps);
    /// Replace the mesh, preserving the current cut when cut is nullopt.
    void set(const std::vector<double> &states, const std::vector<double> &controls, const std::vector<double> &tgrid,
             std::optional<double> cut = std::nullopt,
             std::optional<unsigned> max_steps = std::nullopt);

    /// Access mesh data, propagation settings, integrators and derived segment counts.
    [[nodiscard]] const std::vector<double> &get_states() const;
    [[nodiscard]] const std::vector<double> &get_controls() const;
    [[nodiscard]] const std::vector<double> &get_tgrid() const;
    [[nodiscard]] double get_cut() const;
    [[nodiscard]] unsigned get_dim_dynamics() const;
    [[nodiscard]] unsigned get_dim_controls() const;
    [[nodiscard]] std::optional<unsigned> get_max_steps() const;
    [[nodiscard]] const heyoka::taylor_adaptive<double> &get_ta() const;
    [[nodiscard]] const std::optional<heyoka::taylor_adaptive<double>> &get_ta_var() const;
    [[nodiscard]] bool has_ta_var() const;
    [[nodiscard]] unsigned get_nseg() const;
    [[nodiscard]] unsigned get_nseg_fwd() const;
    [[nodiscard]] unsigned get_nseg_bck() const;

    /**
     * @brief Propagate each segment independently and compute its forward or backward defect.
     * @return Flat vector of length nseg * dim_dynamics, ordered by segment then state component.
     * @note A failed propagation restores the segment starting state, which replaces the flow
     * in the defect. This method modifies the internal nominal integrator.
     */
    [[nodiscard]] std::vector<double> compute_defects() const;

    /**
     * @brief Compute sparse defect derivatives with respect to states, controls and tgrid.
     * @return Tuple (grad_states, grad_controls, grad_tgrid) of flat values in exactly the
     * coordinate order returned by defects_grad_sparsity(). Structural entries are retained
     * even when their numerical values are zero.
     * @note With d = dim_dynamics, c = dim_controls and n = nseg, the value counts are
     * n * (d * d + d), n * d * c and 2 * n * d, respectively. The propagated node contributes
     * a full state-transition matrix, and the other node contributes only the diagonal of
     * the signed identity. Backward flow sensitivities enter with a minus sign.
     * @note A variational integrator is required. On propagation failure the flow sensitivities
     * reduce to the identity for states and zero for controls and times. This method modifies
     * the internal variational integrator.
     */
    [[nodiscard]] std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    compute_defects_grad() const;

    /**
     * @brief Replace interior nodes by propagating inward from the current endpoint states.
     * @param ballistic If true, use zero controls without changing the stored controls.
     * Default is false.
     * @note Endpoints and tgrid are preserved. With stored controls and successful propagation,
     * only the segment where the forward and backward meshes meet may have a defect.
     * Ballistic propagation may leave defects on every segment when stored controls are nonzero.
     * Failed propagation restores the starting state. The nominal integrator is modified.
     */
    void set_initial_guess(bool ballistic = false);

    /**
    * @brief Return the cached structural coordinates of the three defect Jacobian blocks.
    * @return Const reference to the tuple (sp_states, sp_controls, sp_tgrid) of zero-based
    * (row, column) pairs, owned by this leg. Structural setters replace its contents;
    * copy the patterns if a snapshot is required.
     * Rows index the flat defects; columns index each block's own states, controls or tgrid
     * vector, without offsets into a combined decision vector.
     * @note For n = nseg, d = dim_dynamics and c = dim_controls, the block shapes are
     * (n * d, (n + 1) * d), (n * d, n * c) and (n * d, n + 1). The coordinate counts are
     * n * (d * d + d), n * d * c and 2 * n * d. Entries are in row-major order within each
     * block, matching compute_defects_grad(). The signed identity contributes only its diagonal.
     * No variational integrator is needed to construct the patterns.
     */
    [[nodiscard]] const std::tuple<sparsity_pattern, sparsity_pattern, sparsity_pattern> &defects_grad_sparsity() const;

    /**
    * @brief Sample each segment independently from its corresponding mesh node.
    * Forward segments are ordered from node 0 toward the cut and sampled in increasing time.
    * Backward segments are ordered from the final node toward the cut and sampled in decreasing time.
     *
    * @param N Number of sampling points per segment (including endpoints). Default is 5.
    * @return Tuple (state_fwd, state_bck, success). Each entry is an N x dim_dynamics sequence,
    * possibly shorter on failure. Each direction stops at its first incomplete segment;
    * success is true only if all requested segments return all N samples.
     * @note This method modifies the internal state of the nominal integrator.
     */
    std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<std::vector<double>>>, bool>
    get_state_info(unsigned N = 5) const;

private:
    void update_nseg();
    void initialize_ic_var();
    void update_pars_no_control();
    void update_sparsity();
    void sanity_checks() const;

    std::vector<double> m_states;
    std::vector<double> m_controls;
    std::vector<double> m_tgrid;
    double m_cut = 0.5;
    std::optional<unsigned> m_max_steps;
    unsigned m_dim_dynamics = 7u;
    unsigned m_dim_controls = 4u;

    mutable heyoka::taylor_adaptive<double> m_ta;
    mutable std::optional<heyoka::taylor_adaptive<double>> m_ta_var;

    std::vector<double> m_pars_no_control;
    std::vector<double> m_ic_var;

    unsigned m_nseg = 0u;
    unsigned m_nseg_fwd = 0u;
    unsigned m_nseg_bck = 0u;

    std::tuple<sparsity_pattern, sparsity_pattern, sparsity_pattern> m_defects_grad_sparsity;

    heyoka::cfunc<double> m_dyn_cfunc;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & m_states;
        ar & m_controls;
        ar & m_tgrid;
        ar & m_cut;
        ar & m_max_steps;
        ar & m_dim_dynamics;
        ar & m_dim_controls;
        ar & m_ta;
        ar & m_ta_var;
        ar & m_pars_no_control;
        ar & m_ic_var;
        ar & m_nseg;
        ar & m_nseg_fwd;
        ar & m_nseg_bck;
        ar & m_dyn_cfunc;
        if constexpr (Archive::is_loading::value) {
            update_sparsity();
        }
    }
};

/// Stream a summary of the multiple-shooting leg.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const zoh_ms &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::zoh_ms> : fmt::ostream_formatter {
};

#endif // kep3_LEG_ZOH_MS_H
