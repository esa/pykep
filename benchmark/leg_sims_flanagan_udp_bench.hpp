// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_TEST_LEG_SIMS_FLANAGAN_UDP_BENCH_H
#define kep3_TEST_LEG_SIMS_FLANAGAN_UDP_BENCH_H

#include <array>
#include <vector>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>

#include <xtensor/views/xview.hpp>

#include <pagmo/utils/gradients_and_hessians.hpp>

#include <kep3/core_astro/constants.hpp>
#include <kep3/leg/sims_flanagan.hpp>

struct sf_bench_udp {
    sf_bench_udp() = default;
    sf_bench_udp(std::array<std::array<double, 3>, 2> rvs, double ms, std::array<std::array<double, 3>, 2> rvf,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                double max_thrust, double isp, unsigned nseg, bool analytical)
        : m_rvs(rvs), m_rvf(rvf), m_ms(ms), m_max_thrust(max_thrust), m_isp(isp), m_nseg(nseg), m_analytical(analytical)
    {
    }

    [[nodiscard]] std::vector<double> fitness(const std::vector<double> &x) const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 3];// * kep3::DAY2SEC; // in s
        double mf = x[m_nseg * 3 + 1];              // in kg
        kep3::leg::sims_flanagan leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, mf, tof, m_max_thrust,
                                     m_isp*kep3::G0, 1);

        // We set the throttles
        leg.set_throttles(x.begin(), x.end() - 2);

        std::vector<double> retval(1 + 7 + m_nseg, 0.);
        // Fitness
        retval[0] = -mf;
        // Equality Constraints
        auto eq_con = leg.compute_mismatch_constraints();
        retval[1] = eq_con[0]; // / kep3::AU;
        retval[2] = eq_con[1]; // / kep3::AU;
        retval[3] = eq_con[2]; // / kep3::AU;
        retval[4] = eq_con[3]; // / kep3::EARTH_VELOCITY;
        retval[5] = eq_con[4]; // / kep3::EARTH_VELOCITY;
        retval[6] = eq_con[5]; // / kep3::EARTH_VELOCITY;
        retval[7] = eq_con[6]; // / 1e8; //
        //  Inequality Constraints
        auto ineq_con = leg.compute_throttle_constraints();
        std::copy(ineq_con.begin(), ineq_con.end(), retval.begin() + 8);
        return retval;
    }

     [[nodiscard]] std::vector<double> gradient(const std::vector<double> &x) const
     {
        if (m_analytical) {
            return _gradient_analytical(x);
        } else {
            return _gradient_numerical(x);
        }
     }

    [[nodiscard]] std::vector<double> _gradient_numerical(const std::vector<double> &x) const
    {
        return pagmo::estimate_gradient([this](const std::vector<double> &x) { return this->fitness(x); }, x);
    }

    [[nodiscard]] std::vector<double> _gradient_analytical(const std::vector<double> &x) const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // We set the leg (avoiding the allocation for the throttles is possible but requires mutable data members.)
        double tof = x[m_nseg * 3]; // * kep3::DAY2SEC; // in s
        double mf = x[m_nseg * 3 + 1];              // in kg
        kep3::leg::sims_flanagan leg(m_rvs, m_ms, std::vector<double>(m_nseg * 3, 0.), m_rvf, mf, tof, m_max_thrust,
                                     m_isp*kep3::G0, 1);
        // We set the throttles
        leg.set_throttles(x.begin(), x.end() - 2);

        // We compute the gradients
        auto grad_mc_all = (leg.compute_mc_grad());
        auto grad_tc = leg.compute_tc_grad();
        auto grad_mc_xf = std::get<1>(grad_mc_all);
        auto grad_mc = std::move(std::get<2>(grad_mc_all));

        // We allocate the final gradient containing all
        std::vector<double> gradient((1u + 7u + m_nseg) * (m_nseg * 3u + 2u), 0);
        // Create the various xtensor objects adapting the std containers
        auto xgradient
            = xt::adapt(gradient, {1u + 7u + static_cast<unsigned>(m_nseg), static_cast<unsigned>(m_nseg) * 3u + 2u});
        auto xgrad_mc = xt::adapt(grad_mc, {7u, static_cast<unsigned>(m_nseg) * 3u + 1u});
        auto xgrad_mc_xf = xt::adapt(grad_mc_xf, {7u, 7u});
        auto xgrad_tc = xt::adapt(grad_tc, {static_cast<unsigned>(m_nseg), static_cast<unsigned>(m_nseg) * 3u});

        // Row 1 - fitness gradient
        xgradient(0, m_nseg * 3 + 1) = -1.; // fitness gradient - obj fun
        // [1:4,:-1] - fitness gradient - position mismatch
        xt::view(xgradient, xt::range(1u, 4u), xt::range(0, m_nseg * 3u + 1u))
            = xt::view(xgrad_mc, xt::range(0u, 3u), xt::all()); // / kep3::AU; // throttles, tof
        // [4:7,:-1] - fitness gradient - velocity mismatch
        xt::view(xgradient, xt::range(4u, 7u), xt::range(0, m_nseg * 3u + 1u))
            = xt::view(xgrad_mc, xt::range(3u, 6u), xt::all()); // / kep3::EARTH_VELOCITY; // throttles, tof
        // [7:8,:-1] - fitness gradient - mass mismatch
        xt::view(xgradient, xt::range(7u, 8u), xt::range(0, static_cast<unsigned>(m_nseg) * 3u + 1))
            = xt::view(xgrad_mc, xt::range(6u, 7u), xt::all()); // / 1e8; // throttles, tof
        // [8:,:-2] - fitness gradient - throttle constraints
        xt::view(xgradient, xt::range(8u, 8u + static_cast<unsigned>(m_nseg)),
                 xt::range(0, static_cast<unsigned>(m_nseg) * 3u))
            = xgrad_tc;

        // [1:4,-1] - fitness gradient, position mismatch w.r.t. mf
        xt::view(xgradient, xt::range(1u, 4u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_mc_xf, xt::range(0u, 3u), xt::range(6u, 7u)); // / kep3::AU; // mf
        // [4:7,-1] - fitness gradient - velocity mismatch w.r.t. mf
        xt::view(xgradient, xt::range(4u, 7u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_mc_xf, xt::range(3u, 6u), xt::range(6u, 7u)); // / kep3::EARTH_VELOCITY; // mf
        // [7:8,-1] - fitness gradient - mass mismatch w.r.t. mf
        xt::view(xgradient, xt::range(7u, 8u), xt::range(m_nseg * 3u + 1u, m_nseg * 3u + 2u))
            = xt::view(xgrad_mc_xf, xt::range(6u, 7u), xt::range(6u, 7u)); // / 1e8; // mf

        // Units for the tof
        xt::view(xgradient, xt::all(), xt::range(m_nseg * 3u, m_nseg * 3u + 1u)); // *= kep3::DAY2SEC;
        return gradient;
    }

    [[nodiscard]] std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        // x = [throttles, tof (in days), mf (in kg)]
        // std::vector<double> lb(m_nseg * 3 + 2, -1.);
        // std::vector<double> ub(m_nseg * 3 + 2, +1.);
        // lb[m_nseg * 3] = 1.;            // days
        // ub[m_nseg * 3] = 2500.;         // days
        // lb[m_nseg * 3 + 1] = m_ms / 2.; // kg
        // ub[m_nseg * 3 + 1] = m_ms;      // kg
        // return {lb, ub};
        std::vector<double> lb(m_nseg * 3 + 2, -1.);
        std::vector<double> ub(m_nseg * 3 + 2, +1.);
        lb[m_nseg * 3] = kep3::pi / 12;     // days
        ub[m_nseg * 3] = 2 * kep3::pi;     // days
        lb[m_nseg * 3 + 1] = 0.5; // kg
        ub[m_nseg * 3 + 1] = 1;   // kg
        return {lb, ub};
    }

    [[nodiscard]] static std::vector<double>::size_type get_nec()
    {
        return 7u;
    }

    [[nodiscard]] std::vector<double>::size_type get_nic() const
    {
        return m_nseg;
    }

    std::array<std::array<double, 3>, 2> m_rvs{};
    std::array<std::array<double, 3>, 2> m_rvf{};
    double m_ms{};
    double m_max_thrust{};
    double m_isp{};
    std::size_t m_nseg{};
    bool m_analytical{};
};

#endif