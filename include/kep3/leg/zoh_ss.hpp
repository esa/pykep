// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_ZOH_SS_H
#define kep3_LEG_ZOH_SS_H

#include <array>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/ostream.h>

#include <heyoka/taylor.hpp>

#include <kep3/detail/visibility.hpp>

namespace kep3::leg
{
/// The Zero-Order-Hold (ZOH) solar-sail leg model.
/**
 * This class implements a zero-order-hold solar-sail leg between a starting and final six-dimensional state.
 * The transfer is modelled as a sequence of non-uniform segments with controls [alpha, beta] fixed on each segment.
 *
 * The dynamics are provided by a user-supplied pair of compatible Taylor-adaptive integrators.
 * The nominal integrator must have state dimension 6 and at least 2 parameters, with the first
 * two parameters representing [alpha, beta]. When provided, the variational integrator must expose
 * first-order variations w.r.t. the 6 state variables and these 2 controls.
 */
class kep3_DLL_PUBLIC zoh_ss
{

public:
    // Default Constructor.
    zoh_ss() = default;

    // Constructor
    zoh_ss(const std::array<double, 6> &state0, const std::vector<double> &controls,
           const std::array<double, 6> &state1, const std::vector<double> &tgrid, double cut,
           const std::pair<heyoka::taylor_adaptive<double>, std::optional<heyoka::taylor_adaptive<double>>> &tas,
           std::optional<unsigned> max_steps = std::nullopt);

    // Setters
    void set_state0(const std::array<double, 6> &state0);
    void set_state1(const std::array<double, 6> &state1);
    void set_controls(const std::vector<double> &controls);
    void set_tgrid(const std::vector<double> &tgrid);
    void set_cut(double cut);
    void set_max_steps(std::optional<unsigned> max_steps);
    void set(const std::array<double, 6> &state0, const std::vector<double> &controls,
             const std::array<double, 6> &state1, const std::vector<double> &tgrid, double cut,
             std::optional<unsigned> max_steps = std::nullopt);

    // Getters
    [[nodiscard]] const std::array<double, 6> &get_state0() const;
    [[nodiscard]] const std::array<double, 6> &get_state1() const;
    [[nodiscard]] const std::vector<double> &get_controls() const;
    [[nodiscard]] const std::vector<double> &get_tgrid() const;
    [[nodiscard]] double get_cut() const;
    [[nodiscard]] std::optional<unsigned> get_max_steps() const;
    [[nodiscard]] const heyoka::taylor_adaptive<double> &get_ta() const;
    [[nodiscard]] const std::optional<heyoka::taylor_adaptive<double>> &get_ta_var() const;
    [[nodiscard]] bool has_ta_var() const;
    [[nodiscard]] unsigned get_nseg() const;
    [[nodiscard]] unsigned get_nseg_fwd() const;
    [[nodiscard]] unsigned get_nseg_bck() const;

    // Compute constraints
    [[nodiscard]] std::array<double, 6> compute_mismatch_constraints() const;

    // Compute mismatch constraint gradients.
    // Returns flattened row-major matrices with shapes:
    // dmc/dx0: 6 x 6, dmc/dx1: 6 x 6, dmc/du: 6 x (2 * nseg), dmc/dtgrid: 6 x (nseg + 1)
    [[nodiscard]] std::tuple<std::array<double, 36>, std::array<double, 36>, std::vector<double>, std::vector<double>>
    compute_mc_grad() const;

    /**
     * Returns state histories sampled along each ZOH segment, for both forward and backward propagation.
     * On propagation failure, the integration state is restored and the method returns partial trajectories
     * together with success = false.
     *
     * @param N Number of sampling points per segment (including endpoints). Default is 2.
     * @return Tuple (state_fwd, state_bck, success).
     */
    std::tuple<std::vector<std::vector<std::array<double, 6>>>, std::vector<std::vector<std::array<double, 6>>>, bool>
    get_state_info(unsigned N = 2) const;

private:
    void update_nseg();
    void update_ic_var();
    void update_pars_no_control();
    void sanity_checks() const;

    // Constructor args storage
    std::array<double, 6> m_state0{{1., 0., 0., 0., 1., 0.}};
    std::vector<double> m_controls{0.2, 0.1, 0.1, 0.2, 0.3, 0.4};
    std::array<double, 6> m_state1{{1.1, 0.1, 0., 0., 0.9, 0.1}};
    std::vector<double> m_tgrid{0., 0.4, 0.9, 1.3};
    double m_cut = 0.5;
    std::optional<unsigned> m_max_steps;

    // Taylor-adaptive integrators
    mutable heyoka::taylor_adaptive<double> m_ta;
    mutable std::optional<heyoka::taylor_adaptive<double>> m_ta_var;

    // Derived quantities
    std::vector<double> m_pars_no_control;
    std::array<double, 48> m_ic_var{};

    // Cached segment counts
    unsigned m_nseg = 3u;
    unsigned m_nseg_fwd = 1u;
    unsigned m_nseg_bck = 2u;

    // Compiled function for dynamics evaluation.
    heyoka::cfunc<double> m_dyn_cfunc;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar & m_state0;
        ar & m_state1;
        ar & m_controls;
        ar & m_tgrid;
        ar & m_cut;
        ar & m_ta;
        ar & m_ta_var;
        ar & m_pars_no_control;
        ar & m_ic_var;
        ar & m_nseg;
        ar & m_nseg_fwd;
        ar & m_nseg_bck;
        ar & m_dyn_cfunc;
    }
};

// Streaming operator for the class kep3::leg::zoh_ss.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const zoh_ss &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::zoh_ss> : fmt::ostream_formatter {
};

#endif // kep3_LEG_ZOH_SS_H
