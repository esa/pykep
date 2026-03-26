// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_LEG_SIMS_FLANAGAN_HF_H
#define kep3_LEG_SIMS_FLANAGAN_HF_H

#include <array>
#include <fmt/ostream.h>
#include <optional>
#include <utility>
#include <vector>

#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/epoch.hpp>
#include <kep3/ta/zero_hold_kep.hpp>

namespace kep3::leg
{
/// The High-fidelity leg model expanding upon Sims-Flanagan
/**
 * This class represents, generically, a low-thrust leg as a sequence of successive
 * fixed low-thrust segments.
 * The leg achieves to transfer a given spacecraft from an initial to a final state in the
 * time given (and can be considered as feasible) whenever the method evaluate_mismatch
 * returns all zeros and the method get_throttles_con returns all values less than zero.
 *
 * The dynamics is by default Keplerian, but the user can construct his own and pass it to the
 * leg. The user constructed dynamics will need to meet the requirements of having 7 variables
 * with the mass being the last one (zero_hold tas) and five parameters that represent mu, veff and the thrust direction.
 */
class kep3_DLL_PUBLIC sims_flanagan_hf
{
public:
    // Constructors
    // Default Constructor.
    sims_flanagan_hf(); // = default;
    // Main constructor from rvm states
    sims_flanagan_hf(
        const std::array<double, 7> &rvms, const std::vector<double> &throttles, const std::array<double, 7> &rvmf,
        double tof, double max_thrust, double veff, double mu, double cut, double tol = 1e-16,
        std::optional<std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &>>
        = std::nullopt);
    // Backwards-compatible constructor with rv and m states separately
    sims_flanagan_hf(
        const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
        const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double veff,
        double mu, double cut = 0.5, double tol = 1e-16,
        std::optional<std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &>>
        = std::nullopt);

    // Setters
    void set_tof(double tof);
    void set_rvs(const std::array<std::array<double, 3>, 2> &rv);
    void set_ms(double mass);
    void set_throttles(const std::vector<double> &throttles);
    void set_throttles(const std::vector<double>::const_iterator &it1, const std::vector<double>::const_iterator &it2);
    void set_rvf(const std::array<std::array<double, 3>, 2> &rv);
    void set_mf(double mass);
    void set_max_thrust(double max_thrust);
    void set_veff(double veff);
    void set_mu(double mu);
    void set_cut(double cut);
    void set_tol(double tol);
    void set_rvms(const std::array<double, 7> &rvms);
    void set_rvmf(const std::array<double, 7> &rvmf);
    // void set_tas(const heyoka::taylor_adaptive<double> &tas);
    // void set_tas_var(const heyoka::taylor_adaptive<double> &tas_var);
    // Backwards-compatible setting function with rv and m states separately
    void set(const std::array<std::array<double, 3>, 2> &rvs, double ms, const std::vector<double> &throttles,
             const std::array<std::array<double, 3>, 2> &rvf, double mf, double tof, double max_thrust, double veff,
             double mu, double cut = 0.5, double tol = 1e-16);
    // Setting function with rvm states
    void set(const std::array<double, 7> &rvms, const std::vector<double> &throttles, const std::array<double, 7> &rvmf,
             double tof, double max_thrust, double veff, double mu, double cut = 0.5, double tol = 1e-16);
    void set(const std::array<double, 7> &rvms, const std::vector<double> &throttles, const std::array<double, 7> &rvmf,
             double time_of_flight);

    // Getters
    [[nodiscard]] double get_tof() const;
    [[nodiscard]] std::array<std::array<double, 3>, 2> get_rvs() const;
    [[nodiscard]] std::array<std::array<double, 3>, 2> get_rvf() const;
    [[nodiscard]] double get_ms() const;
    [[nodiscard]] const std::vector<double> &get_throttles() const;
    [[nodiscard]] double get_mf() const;
    [[nodiscard]] double get_max_thrust() const;
    [[nodiscard]] double get_veff() const;
    [[nodiscard]] double get_mu() const;
    [[nodiscard]] double get_cut() const;
    [[nodiscard]] double get_tol() const;
    [[nodiscard]] unsigned get_nseg() const;
    [[nodiscard]] unsigned get_nseg_fwd() const;
    [[nodiscard]] unsigned get_nseg_bck() const;
    [[nodiscard]] const heyoka::taylor_adaptive<double> &get_tas() const;
    [[nodiscard]] const heyoka::taylor_adaptive<double> &get_tas_var() const;
    [[nodiscard]] const std::array<double, 7> &get_rvms() const;
    [[nodiscard]] const std::array<double, 7> &get_rvmf() const;

    // Compute constraints
    [[nodiscard]] std::array<double, 7> compute_mismatch_constraints() const;
    [[nodiscard]] std::vector<double> compute_throttle_constraints() const;
    [[nodiscard]] std::vector<double> compute_constraints() const;
    // [[nodiscard]] std::vector<double> set_and_compute_constraints(const std::vector<double> &chromosome);

    // Get state derivative
    [[nodiscard]] std::array<double, 7> get_state_derivative(const std::array<double, 7> &state,
                                                             const std::array<double, 3> &throttles) const;

    // Compute all gradients w.r.t. all legs
    [[nodiscard]]
    std::tuple<std::vector<std::array<double, 49u>>, std::vector<std::array<double, 21u>>,
               std::vector<std::array<double, 7u>>> compute_all_gradients() const;

    // Process all gradients to retrieve relevant gradients (w.r.t. initial and final rvm state as well as w.r.t.
    // throttles and tof)
    [[nodiscard]] std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
    get_relevant_gradients(const std::vector<std::array<double, 49u>> &dxdx_per_seg,
                           const std::vector<std::array<double, 21u>> &dxdu_per_seg,
                           const std::vector<std::array<double, 7u>> &dxdtof_per_seg) const;

    // Compute mismatch constraint gradients (w.r.t. initial and final rvm state as well as w.r.t. throttles and
    // tof)
    [[nodiscard]] std::tuple<std::array<double, 49>, std::array<double, 49>, std::vector<double>>
    compute_mc_grad() const;

    // Compute throttle constraint gradients
    [[nodiscard]] std::vector<double> compute_tc_grad() const;

    // Retrieve the state history of the sims flanagan leg
    [[nodiscard]] std::vector<std::vector<double>> get_state_history(unsigned grid_points_per_segment) const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned)
    {
        ar & m_rvms;
        ar & m_vars;
        ar & m_throttles;
        ar & m_thrusts;
        ar & m_tof;
        ar & m_rvmf;
        ar & m_max_thrust;
        ar & m_veff;
        ar & m_mu;
        ar & m_cut;
        ar & m_tol;
        ar & m_nseg;
        ar & m_nseg_fwd;
        ar & m_nseg_bck;
        ar & m_ta;
        ar & m_ta_var;
        ar & m_cf_dyn;
    }

    // Initial rvm state
    std::array<double, 7> m_rvms{1., 0., 0., 0., 1., 0., 1.};
    // Initial variational state
    std::array<double, 70> m_vars{};
    // Sequence of throttles.
    std::vector<double> m_throttles{0., 0., 0., 0., 0., 0.};
    // Sequence of thrusts.
    std::vector<double> m_thrusts{0., 0., 0., 0., 0., 0.};
    // Final rvm state
    std::array<double, 7> m_rvmf{0., 1., 0., -1., 0., 0., 1.};
    // Time of flight (defaults to 1/4 of the period)
    double m_tof = kep3::pi / 2;
    // Spacecraft propulsion system maximum thrust.
    double m_max_thrust{1};
    // Spacecraft effective velocity Isp*G0 (in units of the dynamics)
    double m_veff{1.};
    // Spacecraft gravitational parameter.
    double m_mu{1.};
    // The cut parameter
    double m_cut = 0.5;
    // The tolerance
    double m_tol = 1e-16;
    // Segment sizes
    unsigned m_nseg = 2u;
    unsigned m_nseg_fwd = 1u;
    unsigned m_nseg_bck = 1u;
    // Taylor-adaptive integrator
    // m_ta needs to be mutable because the heyoka integrator needs to be modifiable
    mutable heyoka::taylor_adaptive<double> m_ta{};
    // Variational Taylor-adaptive integrator
    // m_ta_var needs to be mutable because the heyoka integrator needs to be modifiable
    mutable heyoka::taylor_adaptive<double> m_ta_var{};
    // Dynamics c_func needed for the derivatives w.r.t. tof
    mutable heyoka::cfunc<double> m_cf_dyn{};
};

// Streaming operator for the class kep3::leg::sims_flanagan.
kep3_DLL_PUBLIC std::ostream &operator<<(std::ostream &, const sims_flanagan_hf &);

} // namespace kep3::leg

template <>
struct fmt::formatter<kep3::leg::sims_flanagan_hf> : fmt::ostream_formatter {
};

#endif