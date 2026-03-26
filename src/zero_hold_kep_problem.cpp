// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <heyoka/taylor.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/zero_hold_kep_problem.hpp>
#include <kep3/ta/zero_hold_kep.hpp>

namespace kep3
{

using heyoka::taylor_outcome;
using ta::get_ta_zero_hold_kep;
using ta::get_ta_zero_hold_kep_var;

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
zero_hold_kep_problem::zero_hold_kep_problem(double mu, double veff, double tol)
    : m_mu(mu), m_veff(veff), m_tol(tol), m_ta(get_ta_zero_hold_kep(tol)), m_ta_var(get_ta_zero_hold_kep_var(tol)), m_var_ic(70)
{
    // We set mu and veff for the non variational
    *(m_ta.get_pars_data()) = m_mu;
    *(m_ta.get_pars_data() + 1) = m_veff;
    // ... and variational version of the integrator
    *(m_ta_var.get_pars_data()) = m_mu;
    *(m_ta_var.get_pars_data() + 1) = m_veff;
    // We copy the initial conditions for the variational equations
    std::copy(m_ta_var.get_state().begin()+7, m_ta_var.get_state().end(), m_var_ic.begin());
};

std::array<double, 7> zero_hold_kep_problem::propagate(const std::array<double, 7> &rvm_state, std::array<double, 3> thrust, double tof)
{
    if (rvm_state[6]==0) {
        throw std::domain_error("zero_hold_kep_problem: initial mass is zero!");
    }
    std::array<double, 7> retval{};
    // We set the thrust magnitude.
    std::copy(thrust.begin(), thrust.end(), m_ta.get_pars_data() + 2);
    // Set the Taylor Integration initial conditions
    m_ta.set_time(0.);
    std::copy(rvm_state.begin(), rvm_state.end(), m_ta.get_state_data());
    // ... and integrate
    auto out = m_ta.propagate_until(tof);
    if (std::get<0>(out) != taylor_outcome::time_limit) {
        throw std::domain_error("zero_hold_kep_problem: failiure to reach the final time requested during a propagation.");
    }
    std::copy(m_ta.get_state().begin(), m_ta.get_state().end(), retval.begin());
    return retval;
}

std::tuple<std::array<double, 7>, std::array<double, 49>, std::array<double, 21>>
zero_hold_kep_problem::propagate_var(const std::array<double, 7> &rvm_state, std::array<double, 3> thrust, double tof)
{
    if (rvm_state[6]==0) {
        throw std::domain_error("zero_hold_kep_problem: initial mass is zero!");
    }
    std::array<double, 7> retval{};
    std::array<double, 49> dxdx{};
    std::array<double, 21> dxdu{};
    // We set the thrust magnitude.
    std::copy(thrust.begin(), thrust.end(), m_ta_var.get_pars_data() + 2);
    // Set the Taylor Integration initial conditions
    m_ta_var.set_time(0.);
    std::copy(rvm_state.begin(), rvm_state.end(), m_ta_var.get_state_data());
    std::copy(m_var_ic.begin(), m_var_ic.end(), m_ta_var.get_state_data()+7);
    // ... and integrate
    auto out = m_ta_var.propagate_until(tof);
    if (std::get<0>(out) != taylor_outcome::time_limit) {
        fmt::print("State: {}", m_ta_var.get_state());
        throw std::domain_error("zero_hold_kep_problem: failiure to reach the final time requested during a variational propagation.");
    }
    // We now copy the result into the various return values
    std::copy(m_ta_var.get_state().begin(), m_ta_var.get_state().begin() + 7, retval.begin());
    for (auto i = 0; i < 7; ++i) {
        std::copy(m_ta_var.get_state().begin() + 7 + 10l * i, m_ta_var.get_state().begin() + 7 + 10l * i + 7,
                  dxdx.begin() + 7l * i);
        std::copy(m_ta_var.get_state().begin() + 14 + 10l * i, m_ta_var.get_state().begin() + 14 + 10l * i + 3,
                  dxdu.begin() + 3l * i);
    }
    return {retval, dxdx, dxdu};
}

double zero_hold_kep_problem::get_mu() const
{
    return m_mu;
}
double zero_hold_kep_problem::get_veff() const
{
    return m_veff;
}
double zero_hold_kep_problem::get_tol() const
{
    return m_tol;
}
heyoka::taylor_adaptive<double> zero_hold_kep_problem::get_ta() const
{
    return m_ta;
}
heyoka::taylor_adaptive<double> zero_hold_kep_problem::get_ta_var() const
{
    return m_ta_var;
}

// Setters
void zero_hold_kep_problem::set_mu(double mu_in)
{
    m_mu = mu_in;
    *(m_ta.get_pars_data()) = m_mu;
    *(m_ta_var.get_pars_data()) = m_mu;
};
void zero_hold_kep_problem::set_veff(double veff_in)
{
    m_veff = veff_in;
    *(m_ta.get_pars_data()+1) = m_veff;
    *(m_ta_var.get_pars_data()+1) = m_veff;
};
//Streaming operator
std::ostream &operator<<(std::ostream &os, const zero_hold_kep_problem &p)
{
    os << "zero_hold_kep Problem:\n";
    os << fmt::format("\tmu central body: {}\n", p.get_mu());
    os << fmt::format("\tpropulsion veff (Isp g0): {}\n", p.get_veff());
    os << fmt::format("\ttolerance requested: {}\n", p.get_tol());
    return os;
}

} // namespace kep3