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

#include <kep3/leg/zoh_ss.hpp>
#include <kep3/ta/zoh_ss.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

namespace
{

struct zoh_ss_reference_case {
    std::array<double, 6> state0;
    std::vector<double> controls;
    std::array<double, 6> state1;
    std::vector<double> tgrid;
    double cut;
    heyoka::taylor_adaptive<double> ta;
    heyoka::taylor_adaptive<double> ta_var;
};

zoh_ss_reference_case make_reference_case()
{
    constexpr double tol = 1e-12;

    auto ta = kep3::ta::get_ta_zoh_ss(tol);
    auto ta_var = kep3::ta::get_ta_zoh_ss_var(tol);
    ta.get_pars_data()[2] = 0.01;
    ta_var.get_pars_data()[2] = 0.01;

    return {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0},
            {0.2, 0.1, 0.25, 0.2, 0.15, 0.3, 0.1, 0.4},
            {0.8, 0.2, 0.0, -0.2, 1.1, 0.1},
            {0.0, 0.3, 0.6, 0.9, 1.2},
            0.5,
            std::move(ta),
            std::move(ta_var)};
}

} // namespace

TEST_CASE("compute_mismatch_constraints")
{
    auto data = make_reference_case();

    kep3::leg::zoh_ss leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    const auto mc = leg.compute_mismatch_constraints();
    const std::array<double, 6> expected = {
        0.1714872409652637,
        1.0157421502952377,
        0.05352046338689242,
        -1.2277053650770107,
        -0.11577459801172241,
        -0.06533123154780118,
    };

    REQUIRE(kep3_tests::L_infinity_norm_rel(mc, expected) < 1e-12);
}

TEST_CASE("compute_mc_grad")
{
    auto data = make_reference_case();

    kep3::leg::zoh_ss leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, data.ta_var}};

    const auto [dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid] = leg.compute_mc_grad();

    const std::array<double, 36> ref_dmc_dx0 = {
        1.344552046764875, 0.0975315679037739, -0.000348427626414218, 0.6625903947593111, 0.029043296928082486,
        4.279450003309203e-21, 0.10831715289054136, 0.8571717833466151, 7.521400808045772e-21,
        0.03012803945055804, 0.5759445549709923, -0.0003484276264142181, -0.0006868162678417112,
        -0.0001341736460265959, 0.8270437438960572, -0.00013417364602659587, -4.1151956922114735e-05,
        0.5650588268555373, 1.1053783907674073, 0.45753529913779695, -0.0011831852975620592,
        1.2850187652019927, 0.18383478917284035, -2.54329439102585e-20, 0.5456295604786638,
        -0.3643583928719417, -3.4827062579528104e-20, 0.1945739530948616, 0.9163982663891408,
        -0.0011831852975620592, -0.0022955305782998213, -0.0006511832704030957, -0.5589323459668033,
        -0.0006511832704030957, -0.00026613721619460107, 0.8274834660641959,
    };

    const std::array<double, 36> ref_dmc_dx1 = {
        -1.6876565785774478, 0.011044175269084128, 0.023060607174478016, 0.7302516579249552,
        -0.031134221734963833, -0.007023221462174532, 0.062209783592802395, -0.7174059197003063,
        -0.004038506637931591, -0.03653325746972896, 0.5468819683434383, 0.0021126021831381174,
        0.02785942297501654, -0.0030030851328812767, -0.6867096997891968, -0.007574581344878912,
        0.001788433890292376, 0.5347310061696088, 2.346037625044413, -0.3838929234821028,
        -0.10738275462733564, -1.612802237423542, 0.27251266647708927, 0.044298015975614585,
        -0.7996130712422133, -0.7888663074955842, 0.035938197940533746, 0.3259858743471086,
        -0.8049093720981736, -0.01924014541141937, -0.14480125344545086, 0.02743048878950537,
        -0.9991429018090692, 0.04945104273126944, -0.017457526418509655, -0.6774320330508642,
    };

    const std::vector<double> ref_dmc_dcontrols = {
        -0.0008251225281877982, -2.659641284759342e-05, -0.0003186140523337648, -3.920582013169406e-05,
        0.0002267985796451559, -2.8581083847973313e-05, 0.00072542590099947, 2.358194503976297e-05,
        -1.6919871530211055e-05, 0.000246457900643507, -5.976481346142649e-05, 9.34231147119829e-05,
        -0.00022378343262285908, -9.344086816599986e-05, -0.0005709802294965314, -0.00018007016125270238,
        0.001118021284326077, -2.471938377631667e-05, 0.0003457748001760872, -2.0501977024680522e-05,
        -0.0006336828696665956, 2.240165234090839e-05, -0.0017503629786729095, 5.6792694887203023e-05,
        -0.0019924341467552867, -2.890799909395033e-05, -0.002120474510497177, -0.00028765108087785995,
        -0.0014411486481230925, 0.0002257700340366845, -0.0021152154689840865, -0.00010794569345045138,
        -0.00021513852304499016, 0.0005290225740308213, -0.0005235188329563173, 0.0006046634902169105,
        0.0016324795828550087, 0.0006023161234615025, 0.0015021655639689076, 0.00038440118847714674,
        0.002317147046164722, -5.125700195792929e-05, 0.002287200226179161, -0.000135546480068724,
        0.004191037830336969, -0.00014978914896056507, 0.00342950879671592, -0.00010501621276691841,
    };

    const std::vector<double> ref_dmc_dtgrid = {
        0.5588157749327002, 0.00011657103410311898, -1.2277053650770127, 7.986920371516604e-05,
        0.6686931499064942, -0.8274364260480092, -4.7040016186716826e-05, -0.11577459801171774,
        -1.8802849376275965e-05, 0.9432768669252899, -0.001071632293093212, -0.00011155300446884708,
        -0.06533123154780074, -0.00022174739743079097, 0.06673616424279359, 0.8153514847274255,
        0.00039663765089503755, 0.46472026617271456, -0.0002969667010208621, -1.2801714218500142,
        0.5569283216976959, -0.00014073158163330213, -1.4377541692444584, 8.426676441919412e-05,
        0.8808823123639766, -0.001565000136218636, -0.00036082044183270864, -0.1041944034123343,
        0.0006987413912366292, 0.10542148259914902,
    };

    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx0, ref_dmc_dx0) < 1e-12);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx1, ref_dmc_dx1) < 1e-12);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dcontrols, ref_dmc_dcontrols) < 1e-12);
    REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dtgrid, ref_dmc_dtgrid) < 1e-12);
}

TEST_CASE("get_state_info")
{
    auto data = make_reference_case();
    constexpr unsigned N = 5u;

    kep3::leg::zoh_ss leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, std::nullopt}};

    const auto [state_fwd, state_bck, success] = leg.get_state_info(N);

    REQUIRE(success == true);
    REQUIRE(state_fwd.size() == 2u);
    REQUIRE(state_bck.size() == 2u);

    for (const auto &seg : state_fwd) {
        REQUIRE(seg.size() == N);
        for (const auto &s : seg) {
            REQUIRE(s.size() == 6u);
        }
    }

    for (const auto &seg : state_bck) {
        REQUIRE(seg.size() == N);
        for (const auto &s : seg) {
            REQUIRE(s.size() == 6u);
        }
    }

    std::array<double, 6> mismatch{};
    const auto &xfwd = state_fwd.back().back();
    const auto &xbck = state_bck.back().back();
    for (std::size_t i = 0u; i < 6u; ++i) {
        mismatch[i] = xfwd[i] - xbck[i];
    }

    const auto mc = leg.compute_mismatch_constraints();
    REQUIRE(kep3_tests::L_infinity_norm_rel(mismatch, mc) < 1e-12);
}

TEST_CASE("compute_mc_grad_corner_cuts")
{
    constexpr double eps = 1e-7;

    auto make_eye6 = []() {
        std::array<double, 36> I{};
        for (unsigned i = 0u; i < 6u; ++i) {
            I[6u * i + i] = 1.;
        }
        return I;
    };

    auto make_minus_eye6 = []() {
        std::array<double, 36> I{};
        for (unsigned i = 0u; i < 6u; ++i) {
            I[6u * i + i] = -1.;
        }
        return I;
    };

    SECTION("cut_zero")
    {
        auto data = make_reference_case();
        data.cut = 0.;

        kep3::leg::zoh_ss leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, data.ta_var}};
        const auto [dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid] = leg.compute_mc_grad();

        CAPTURE(dmc_dx0[0], dmc_dx0[7], dmc_dx0[14], dmc_dx0[21], dmc_dx0[28], dmc_dx0[35]);

        const auto I = make_eye6();
        REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx0, I) < 1e-13);

        auto data_p = make_reference_case();
        data_p.cut = 0.;
        data_p.state0[0] += eps;
        kep3::leg::zoh_ss leg_p{data_p.state0, data_p.controls, data_p.state1, data_p.tgrid, data_p.cut,
                                {data_p.ta, std::nullopt}};

        auto data_m = make_reference_case();
        data_m.cut = 0.;
        data_m.state0[0] -= eps;
        kep3::leg::zoh_ss leg_m{data_m.state0, data_m.controls, data_m.state1, data_m.tgrid, data_m.cut,
                                {data_m.ta, std::nullopt}};

        const auto mc_p = leg_p.compute_mismatch_constraints();
        const auto mc_m = leg_m.compute_mismatch_constraints();
        const double fd = (mc_p[0] - mc_m[0]) / (2. * eps);
        REQUIRE(std::abs(fd - dmc_dx0[0]) < 1e-6);
        (void)dmc_dx1;
        (void)dmc_dcontrols;
        (void)dmc_dtgrid;
    }

    SECTION("cut_one")
    {
        auto data = make_reference_case();
        data.cut = 1.;

        kep3::leg::zoh_ss leg{data.state0, data.controls, data.state1, data.tgrid, data.cut, {data.ta, data.ta_var}};
        const auto [dmc_dx0, dmc_dx1, dmc_dcontrols, dmc_dtgrid] = leg.compute_mc_grad();

        CAPTURE(dmc_dx1[0], dmc_dx1[7], dmc_dx1[14], dmc_dx1[21], dmc_dx1[28], dmc_dx1[35]);

        const auto minus_I = make_minus_eye6();
        REQUIRE(kep3_tests::L_infinity_norm_rel(dmc_dx1, minus_I) < 1e-13);

        auto data_p = make_reference_case();
        data_p.cut = 1.;
        data_p.state1[0] += eps;
        kep3::leg::zoh_ss leg_p{data_p.state0, data_p.controls, data_p.state1, data_p.tgrid, data_p.cut,
                                {data_p.ta, std::nullopt}};

        auto data_m = make_reference_case();
        data_m.cut = 1.;
        data_m.state1[0] -= eps;
        kep3::leg::zoh_ss leg_m{data_m.state0, data_m.controls, data_m.state1, data_m.tgrid, data_m.cut,
                                {data_m.ta, std::nullopt}};

        const auto mc_p = leg_p.compute_mismatch_constraints();
        const auto mc_m = leg_m.compute_mismatch_constraints();
        const double fd = (mc_p[0] - mc_m[0]) / (2. * eps);
        REQUIRE(std::abs(fd - dmc_dx1[0]) < 1e-6);
        (void)dmc_dx0;
        (void)dmc_dcontrols;
        (void)dmc_dtgrid;
    }
}
