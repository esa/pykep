// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <vector>

#include <kep3/core_astro/constants.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/zoh.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/zoh_kep.hpp>
#include <kep3/udpla/jpl_lp.hpp>
#include <kep3/detail/s11n.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

namespace
{

struct zoh_reference_case {
    std::array<double, 7> state0;
    std::vector<double> controls;
    std::array<double, 7> state1;
    std::vector<double> tgrid;
    double cut;
    heyoka::taylor_adaptive<double> ta;
    heyoka::taylor_adaptive<double> ta_var;
};

zoh_reference_case make_reference_case()
{
    constexpr double tol = 1e-14;
    constexpr double t0 = 1234.;
    constexpr double t1 = 3456.;
    constexpr double m0 = 1000.;
    constexpr double m1 = 1000.;
    constexpr unsigned nseg = 5u;

    kep3::planet pl0{kep3::udpla::jpl_lp{"Venus"}};
    kep3::planet pl1{kep3::udpla::jpl_lp{"Earth"}};

    auto const rv0 = pl0.eph(t0);
    auto const rv1 = pl1.eph(t1);
    kep3::lambert_problem lp{rv0[0], rv1[0], (t1 - t0) * kep3::DAY2SEC, kep3::MU_SUN};

    const double L = kep3::AU;
    const double MU = kep3::MU_SUN;
    const double TIME = std::sqrt(L * L * L / MU);
    const double V = L / TIME;
    const double ACC = V / TIME;
    const double MASS = 1000.;
    const double F = MASS * ACC;
    const double veff = 6000. * kep3::G0;
    const double veff_nd = veff / V;

    zoh_reference_case retval{.state0 = {rv0[0][0] / L, rv0[0][1] / L, rv0[0][2] / L, lp.get_v0()[0][0] / V,
                                         lp.get_v0()[0][1] / V, lp.get_v0()[0][2] / V, m0 / MASS},
                              .controls = {0.03 / F, 1.0, 0.0,      0.0, 0.02 / F, 1.0, 0.0,       0.0, 0.003 / F, 0.0,
                                           1.0,      0.0, 0.03 / F, 1.0, 0.0,      0.0, 0.001 / F, 0.0, 0.0,       1.0},
                              .state1 = {rv1[0][0] / L, rv1[0][1] / L, rv1[0][2] / L, lp.get_v1()[0][0] / V,
                                         lp.get_v1()[0][1] / V, lp.get_v1()[0][2] / V, m1 / MASS},
                              .tgrid = {},
                              .cut = 0.5,
                              .ta = kep3::ta::get_ta_zoh_kep(tol),
                              .ta_var = kep3::ta::get_ta_zoh_kep_var(tol)};

    retval.tgrid.resize(nseg + 1u);
    for (unsigned i = 0u; i <= nseg; ++i) {
        retval.tgrid[i] = (t0 * kep3::DAY2SEC / TIME)
                          + (t1 * kep3::DAY2SEC / TIME - t0 * kep3::DAY2SEC / TIME) * static_cast<double>(i)
                                / static_cast<double>(nseg);
    }

    *(retval.ta.get_pars_data() + 4) = 1. / veff_nd;
    *(retval.ta_var.get_pars_data() + 4) = 1. / veff_nd;
    return retval;
}

} // namespace

TEST_CASE("compute_throttle_constraints")
{
    auto data = make_reference_case();
    kep3::leg::zoh leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    auto const tc = leg.compute_throttle_constraints();
    REQUIRE(tc.size() == 5u);

    std::vector<double> const expected(5u, 0.);
    REQUIRE(kep3_tests::L_infinity_norm_rel(tc, expected) < 1e-15);
}

TEST_CASE("compute_mismatch_constraints")
{
    auto data = make_reference_case();
    kep3::leg::zoh leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    auto const mc = leg.compute_mismatch_constraints();
    std::array<double, 7> const expected = {
        0.02176525074416702,  -0.10717348832206364,  -0.014802731334029318, 0.14380114520228338,
        -0.04229446021128315, -0.004488343342410163, -0.05481446161533288,
    };

    REQUIRE(kep3_tests::L_infinity_norm_rel(mc, expected) < 1e-12);
}

TEST_CASE("compute_mc_grad")
{
    auto data = make_reference_case();

    kep3::leg::zoh leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, data.ta_var}};

    auto [dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid] = leg.compute_mc_grad();

    // Reference values (replace with actual Python output):
    std::array<double, 49> ref_dmc_dx0 = {
        -28.89197973155973, 4.92732369493274, 1.27710619508696, -1.47434073353856, -25.17586592639953, 0.08699425041571, -0.56895096700586,
        71.88229905095362, 5.95649343260037, -4.1544831346279, 18.45194188646853, 56.40711924262009, -0.87956814211855, 0.15061589926087,
        1.56837668910072, -0.2404131408908, -7.60291499393688, 0.25719954093257, 1.71782453076697, 1.60174924383033, 0.01042118644378,
        -6.59447701324815, -0.00890340331431, 0.3797264695942, -1.26874642854301, -5.22371035996367, 0.04855251593817, -0.08467912206302,
        6.32922892709135, 0.06743746385021, -0.37097492202881, 1.26374683453272, 5.17927272463923, -0.05427468947683, 0.0335335538595,
        0.41800461037587, 0.00101825121685, -0.16386117735228, 0.0758472362242, 0.33916202615184, -0.10645617178342, 0.00270148679708,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
    };
    std::array<double, 49> ref_dmc_dx1 = {
        -8.40021864775383, -75.6986205546155, 0.53570236082703, 87.74846729377357, 52.35998454108677, -4.70937868766582, 0.53928327630283,
        -5.7422125022256, -2.75678549222092, 0.31428325431877, -1.17173999566234, 8.71332274748447, 0.08954544977415, -0.03970469749386,
        0.80483950133332, 4.52446209727957, 5.47446861225021, -5.0146978804342, -3.13065263119156, 4.91356441276588, 0.00411232706043,
        0.57984028714607, 7.01341498016066, -0.00561154278267, -8.28205244864012, -4.61432897110805, 0.48428753319856, -0.06278043232123,
        -0.27223624676809, -2.78047563854411, 0.00458735583423, 3.06887753408114, 1.89728985782157, -0.17567800818225, 0.01332964256027,
        -0.03712794795976, -0.45035017428401, 0.0012550540849, 0.54217957161615, 0.29552497006717, 0.1513906071779, 0.00235480052286,
        -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0
    };
    std::vector<double> ref_dmc_dcontrols = {
        91.73476557268101, 0.45763200657108, -0.28135725774724, -0.01031907317353, 31.09418929737732, 0.10440792080321, -0.00533602707709, -0.00045457651689,
        0.87488440350236, -0.01505254086631, 0.00044277737647, 0.00006261294551, -107.24705388291584, -0.54698571687815, 0.05053123459141, 0.0103081991757,
        16.79209386337708, -0.05601179996783, -0.00413638005807, 0.00247869629428, -28.72309170977026, -0.14404629946369, 0.70835140389999, 0.01013584022919,
        -1.57360553816639, -0.00529368909715, 0.1014276654063, 0.00034130505513, -28.26664951403545, 0.00044565526668, -0.01430914720517, -0.00002930931036,
        10.74724997022017, 0.05447548651581, -0.40260086656745, -0.00347096298556, -2.16537828966123, 0.00546872494741, -0.01561391615775, -0.00033912604726,
        -1.96983381789672, -0.00987971969853, 0.01846270259939, 0.30411337468551, -0.13517457261022, -0.00045468343358, 0.00034411567169, 0.09738982749413,
        -0.05758336979501, 0.00006266124429, -0.00002914289658, -0.01410234734889, 2.03676539391028, 0.01037027183734, -0.00324142530578, -0.38428189608695,
        -85.31665463779889, 0.00255640436003, 0.00023628798338, -0.01438546063671, 11.08043299892456, 0.05492934754674, -0.05849793618098, -0.00250935691183,
        8.48706381330209, 0.02843898539217, -0.002592095189, -0.00024169665146, -0.51512924450038, 0.00408431229459, -0.00026073045229, -0.00003362758322,
        12.53383288666958, 0.06395502800351, -0.01539836025854, -0.00252392862746, -2.17639178017643, 0.00552115654423, -0.00012445972877, -0.00032591005929,
        -6.12153308575359, -0.0306588627639, 0.07089021529296, 0.00205764362375, -0.7605391268277, -0.00255674650435, 0.02650220945787, 0.00016538098315,
        7.33908456461765, -0.0002635557282, 0.00371640230435, 0.00001755508815, -3.39462214990647, -0.01723306474664, 0.03203446613325, 0.0011183595583,
        0.7755171292174, -0.00200017859979, 0.00025122641828, 0.0001220302752, -0.48604398023409, -0.00243297443363, 0.00378558235203, 0.01741142675005,
        -0.07193614639063, -0.00024180730376, 0.00016774468822, 0.02471982711497, 0.03436097037797, -0.0000336749338, 0.00001739168815, 0.00358066335047,
        -0.50112854757231, -0.00255366494369, 0.00101179970477, 0.02496849363433, 0.9666097425291, -0.00034556512332, 0.00000895037566, 0.0001614547959,
        -3.86969447003105, 0.0, 0.0, 0.0, -3.86969447003105, 0.0, 0.0, 0.0, -3.86969447003105, -0.0, -0.0, -0.0,
        -3.86969447003106, -0.0, -0.0, -0.0, -3.86969447003105, -0.0, -0.0, -0.0
    };
    std::vector<double> ref_dmc_dtgrid
        = {0.07525339735098,  0.01412034689889,  0.14380114520228,  0.04121585590418,  -0.11047763236567,
           -0.16391311299066, 0.15640064005516,  -0.00150248544324, -0.04229446021128, -0.00591789505836,
           0.01413121620724,  -0.12081701554951, -0.00701113187738, -0.00012100836737, -0.00448834334241,
           -0.00034239368632, 0.0051859488985,   0.00677692837499,  -0.03602776627764, 0.00213308282678,
           0.00459258866154,  -0.00626861634219, 0.01154642401326,  0.02402428711826,  0.01851178663655,
           -0.00053649624179, -0.00054693897454, 0.00146647941931,  -0.00389679114243, -0.0149980396971,
           0.0020298538426,   -0.00004823755368, 0.00001333345595,  0.00013904588755,  -0.00068527212892,
           -0.0014487235035,  0.00256083519394,  -0.00085361173131, -0.00145113994323, 0.00230475167454,
           -0.00247547402081, -0.00008536117313};

    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx0, ref_dmc_dx0) < 1e-10);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx1, ref_dmc_dx1) < 1e-10);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dcontrols, ref_dmc_dcontrols) < 1e-10);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dtgrid, ref_dmc_dtgrid) < 1e-10);
}

TEST_CASE("get_state_info")
{
    auto data = make_reference_case();
    unsigned nseg = static_cast<unsigned>(data.tgrid.size() - 1);
    double cut = data.cut;
    unsigned N = 5;
    unsigned nseg_fwd = static_cast<unsigned>(nseg * cut);
    unsigned nseg_bck = nseg - nseg_fwd;

    kep3::leg::zoh leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    auto [state_fwd, state_bck] = leg.get_state_info(N);

    // Check sizes
    REQUIRE(state_fwd.size() == nseg_fwd);
    REQUIRE(state_bck.size() == nseg_bck);
    for (const auto &seg : state_fwd) {
        REQUIRE(seg.size() == N);
        for (const auto &s : seg)
            REQUIRE(s.size() == 7);
    }
    for (const auto &seg : state_bck) {
        REQUIRE(seg.size() == N);
        for (const auto &s : seg)
            REQUIRE(s.size() == 7);
    }

    // Check mismatch: last state of last fwd seg vs last state of last bck seg
    const auto &xfwd = state_fwd.empty() ? std::array<double, 7>{} : state_fwd.back().back();
    const auto &xbck = state_bck.empty() ? std::array<double, 7>{} : state_bck.back().back();
    std::array<double, 7> mismatch{};
    for (size_t i = 0; i < 7; ++i) {
        mismatch[i] = xfwd[i] - xbck[i];
    }
    auto constraint = leg.compute_mismatch_constraints();
    for (size_t i = 0; i < 7; ++i) {
        REQUIRE(std::abs(mismatch[i] - constraint[i]) < 1e-8);
    }
}

TEST_CASE("serialization")
{
    // Create a reference zoh leg
    auto data = make_reference_case();
    kep3::leg::zoh zoh1{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    // Store the string representation
    std::stringstream ss;
    auto before = boost::lexical_cast<std::string>(zoh1);
    // Serialize
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << zoh1;
    }
    // Deserialize
    kep3::leg::zoh zoh2{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> zoh2;
    }
    auto after = boost::lexical_cast<std::string>(zoh2);
    // Compare the string representations
    REQUIRE(before == after);
}