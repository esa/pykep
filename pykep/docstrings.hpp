// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef PYKEP_DOCSTRINGS_HPP
#define PYKEP_DOCSTRINGS_HPP

#include <string>

namespace pykep
{
// Modules
std::string core_module_doc();

// Anomaly conversions
std::string m2e_doc();
std::string e2m_doc();
std::string m2f_doc();
std::string f2m_doc();
std::string e2f_doc();
std::string f2e_doc();

std::string n2h_doc();
std::string h2n_doc();
std::string n2f_doc();
std::string f2n_doc();
std::string h2f_doc();
std::string f2h_doc();

std::string zeta2f_doc();
std::string f2zeta_doc();

// Anomaly conversions (vectorized)
std::string m2e_v_doc();
std::string e2m_v_doc();
std::string m2f_v_doc();
std::string f2m_v_doc();
std::string e2f_v_doc();
std::string f2e_v_doc();

std::string n2h_v_doc();
std::string h2n_v_doc();
std::string n2f_v_doc();
std::string f2n_v_doc();
std::string h2f_v_doc();
std::string f2h_v_doc();

std::string zeta2f_v_doc();
std::string f2zeta_v_doc();

// Elements conversions
std::string ic2par_doc();
std::string par2ic_doc();
std::string ic2mee_doc();
std::string mee2ic_doc();
std::string par2mee_doc();
std::string mee2par_doc();

// MIMA
std::string mima_doc();
std::string mima_from_hop_doc();
std::string mima2_doc();
std::string mima2_from_hop_doc();
std::string hohmann_doc();
std::string bielliptic_doc();

// Encodings
std::string alpha2direct_doc();
std::string direct2alpha_doc();
std::string eta2direct_doc();
std::string direct2eta_doc();

// Epoch
std::string epoch_from_float_doc();
std::string epoch_from_datetime_doc();
std::string epoch_from_string_doc();

// Planet
std::string planet_docstring();
std::string planet_get_name_docstring();
std::string planet_get_extra_info_docstring();
std::string planet_get_mu_central_body_docstring();
std::string planet_get_mu_self_docstring();
std::string planet_get_radius_docstring();
std::string planet_get_safe_radius_docstring();
std::string planet_eph_docstring();
std::string planet_eph_v_docstring();
std::string planet_acc_docstring();
std::string planet_acc_v_docstring();
std::string planet_period_docstring();
std::string planet_elements_docstring();

// UDPLAS
std::string udpla_keplerian_from_elem_docstring();
std::string udpla_keplerian_from_posvel_docstring();
std::string udpla_jpl_lp_docstring();
std::string udpla_vsop2013_docstring();

// Taylor Adaptive propagators
// basic
std::string get_kep_docstring();
std::string get_kep_var_docstring();
std::string kep_dyn_docstring();
std::string get_cr3bp_docstring();
std::string get_cr3bp_var_docstring();
std::string cr3bp_dyn_docstring();
std::string cr3bp_jacobi_C_docstring();
std::string cr3bp_effective_potential_U_docstring();
std::string get_bcp_docstring();
std::string get_bcp_var_docstring();
std::string bcp_dyn_docstring();
// zero holds
std::string get_zoh_kep_docstring();
std::string get_zoh_kep_var_docstring();
std::string zoh_kep_dyn_docstring();
std::string get_zoh_eq_docstring();
std::string get_zoh_eq_var_docstring();
std::string zoh_eq_dyn_docstring();
std::string get_zoh_cr3bp_docstring();
std::string get_zoh_cr3bp_var_docstring();
std::string zoh_cr3bp_dyn_docstring();
std::string get_zoh_ss_docstring();
std::string get_zoh_ss_var_docstring();
std::string zoh_ss_dyn_docstring();
// TPBVPs
std::string get_pc_docstring();
std::string get_pc_var_docstring();
std::string pc_dyn_docstring();
std::string get_pc_H_cfunc_docstring();
inline std::string get_pc_SF_cfunc_docstring(){return "The Switching Function along an optimal trajectory.";};
inline std::string get_pc_u_cfunc_docstring(){return "The optimal throttle.";};
inline std::string get_pc_i_vers_cfunc_docstring(){return "The optimal thrust direction.";};
inline std::string get_pc_dyn_cfunc_docstring(){return "The augmented dynamics.";};

std::string get_peq_docstring();
std::string get_peq_var_docstring();
std::string peq_dyn_docstring();
inline std::string get_peq_H_cfunc_docstring(){return "The Hamiltonian along an optimal trajectory.";};
inline std::string get_peq_SF_cfunc_docstring(){return "The Switching Function along an optimal trajectory.";};
inline std::string get_peq_u_cfunc_docstring(){return "The optimal throttle.";};
inline std::string get_peq_i_vers_cfunc_docstring(){return "The optimal thrust direction.";};
inline std::string get_peq_dyn_cfunc_docstring(){return "The augmented dynamics.";};


// Lambert Problem
std::string lambert_problem_docstring();

// Flybys
std::string fb_con_docstring();
std::string fb_dv_docstring();
std::string fb_vout_docstring();

// Propagators
std::string propagate_lagrangian_docstring();
std::string propagate_lagrangian_grid_docstring();

// LEG
// Sims Flanagan
std::string leg_sf_docstring();
std::string leg_sf_rvs_docstring();
std::string leg_sf_ms_docstring();
std::string leg_sf_throttles_docstring();
std::string leg_sf_rvf_docstring();
std::string leg_sf_mf_docstring();
std::string leg_sf_tof_docstring();
std::string leg_sf_max_thrust_docstring();
std::string leg_sf_veff_docstring();
std::string leg_sf_mu_docstring();
std::string leg_sf_cut_docstring();
std::string leg_sf_mc_docstring();
std::string leg_sf_tc_docstring();
std::string leg_sf_mc_grad_docstring();

std::string leg_sf_tc_grad_docstring();
std::string leg_sf_nseg_docstring();
std::string leg_sf_nseg_fwd_docstring();
std::string leg_sf_nseg_bck_docstring();
// Alpha
std::string leg_sf_alpha_docstring();
std::string leg_sf_talphas_docstring();
// Zoh
std::string leg_zoh_docstring();
std::string leg_zoh_state0_docstring();
std::string leg_zoh_controls_docstring();
std::string leg_zoh_state1_docstring();
std::string leg_zoh_tgrid_docstring();
std::string leg_zoh_cut_docstring();
std::string leg_zoh_max_steps_docstring();
std::string leg_zoh_nseg_docstring();
std::string leg_zoh_nseg_fwd_docstring();
std::string leg_zoh_nseg_bck_docstring();
std::string leg_zoh_mc_docstring();
std::string leg_zoh_tc_docstring();
std::string leg_zoh_mc_grad_docstring();
std::string leg_zoh_tc_grad_docstring();
std::string leg_zoh_get_state_info_docstring();

} // namespace pykep

#endif