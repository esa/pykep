// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef kep3_SF_CHECKS_H
#define kep3_SF_CHECKS_H

#include <kep3/detail/visibility.hpp>
#include <vector>

// These checks are used for the low- and high-fidelity legs (in sims_flanagan.cpp and sims_flanagan_hf.cpp)
namespace kep3::leg {

kep3_DLL_PUBLIC void _check_tof(double tof);
kep3_DLL_PUBLIC void _check_throttles(const std::vector<double> &throttles);
kep3_DLL_PUBLIC void _check_talphas(const std::vector<double> &talphas,unsigned nseg);
kep3_DLL_PUBLIC void _check_max_thrust(double max_thrust);
kep3_DLL_PUBLIC void _check_veff(double veff);
kep3_DLL_PUBLIC void _check_isp(double isp);
kep3_DLL_PUBLIC void _check_mu(double mu);
kep3_DLL_PUBLIC void _check_cut(double cut);
kep3_DLL_PUBLIC void _check_tol(double tol);
kep3_DLL_PUBLIC void _check_nseg(unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
kep3_DLL_PUBLIC void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
kep3_DLL_PUBLIC void _sanity_checks(const std::vector<double> &throttles, double tof, double max_thrust, double isp, double mu,
                    double cut, double tol, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
kep3_DLL_PUBLIC void _sanity_checks_alpha(const std::vector<double> &throttles, const std::vector<double> &talphas, double tof, double max_thrust, double isp, double mu,
                    double cut, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
kep3_DLL_PUBLIC void _sanity_checks_alpha(const std::vector<double> &throttles, const std::vector<double> &talphas, double tof, double max_thrust, double isp, double mu,
                    double cut, double tol, unsigned nseg, unsigned nseg_fwd, unsigned nseg_bck);
    

} // namespace kep3::leg

#endif
