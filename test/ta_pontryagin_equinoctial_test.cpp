// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <heyoka/taylor.hpp>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/ta/pontryagin_equinoctial.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_peq;
using kep3::ta::get_ta_peq_cache_dim;
using kep3::ta::get_ta_peq_var;
using kep3::ta::get_ta_peq_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    {
        // MASS
        //  The non variational one.
        REQUIRE(get_ta_peq_cache_dim() == 0u);
        auto ta_cached = get_ta_peq(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_cache_dim() == 1u);
        ta_cached = get_ta_peq(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_cache_dim() == 1u);
        ta_cached = get_ta_peq(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_cache_dim() == 2u);

        // The variational integrator.
        REQUIRE(get_ta_peq_var_cache_dim() == 0u);
        ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_var_cache_dim() == 1u);
        ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_var_cache_dim() == 1u);
        ta_cached = get_ta_peq_var(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_peq_var_cache_dim() == 2u);
    }
    {
        // TIME
        //  The non variational one. (no cache is not empty)
        REQUIRE(get_ta_peq_cache_dim() == 2u);
        auto ta_cached = get_ta_peq(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_cache_dim() == 3u);
        ta_cached = get_ta_peq(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_cache_dim() == 3u);
        ta_cached = get_ta_peq(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_cache_dim() == 4u);

        // The variational integrator.
        REQUIRE(get_ta_peq_var_cache_dim() == 2u);
        ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_var_cache_dim() == 3u);
        ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_var_cache_dim() == 3u);
        ta_cached = get_ta_peq_var(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_peq_var_cache_dim() == 4u);
    }
}

TEST_CASE("dynamics_mass")
{
    auto ta_cached = get_ta_peq(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 5); // [mu, c1, c2, eps, l0]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        std::vector<double> pars = {1, 1e-4, 1, 1, 1e-4};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.0);
        std::vector<double> const ground_truth
            = {9.9989683610656938e-02, 2.0000464728286413e-01, 3.0001445961286149e-01, 3.9999965108792562e-01,
               4.9999907790556286e-01, 2.6893895051671013e+01, 6.9990000134984875e-01, 4.2293276260195773e+02,
               1.9220340493191951e+01, 2.8262413631404467e+01, 5.0000094914060389e-01, 5.0000022800154342e-01,
               5.6376609370642772e-01, 4.9711717148509021e-01};

        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
    {
        // We test a case provided by Zhong Zhang tools (Tshinghua)
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1.48102e+11, 0.1, -1.70469e-10, 0.0500417, 3.49378e-18, 4.95038, 1500,
                                  0.0,         0.1, 0.1,          0.1,       0.1,         0.1,     0.0};
        std::vector<double> pars = {1.32712440018e20, 0.6, 3000.0 * 9.80665, 1e-2, 1.0};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(100.0 * 86400.0);
        std::vector<double> const ground_truth
            = {1.4803130514268405e+11,  9.9027557237587743e-02, -1.1380120593980601e-03, 5.0011064027660791e-02,
               -2.1795658780245426e-05, 6.9815960043535021e+00, 1.4982552310395881e+03,  1.8561671535926119e-12,
               -1.6580839961454405e-01, 1.9590460144511043e-01, 1.0000530638532458e-01,  9.9993276760264954e-02,
               9.0789957994717868e-02,  -2.1460001306482017e-07};

        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_mass")
{
    auto ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 5);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        std::vector<double> pars = {1., 1e-4, 1., 1., 1e-4};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.);
        std::vector<double> const ground_truth
            = {9.9989683610656882e-02,  2.0000464728286399e-01,  3.0001445961286144e-01,  3.9999965108792579e-01,
               4.9999907790556286e-01,  2.6893895051671077e+01,  6.9990000134984798e-01,  4.2293276260195813e+02,
               1.9220340493191923e+01,  2.8262413631404364e+01,  5.0000094914060433e-01,  5.0000022800154409e-01,
               5.6376609370643949e-01,  4.9711717148508550e-01,  -4.2696391041655561e-08, 2.6650144402136371e-07,
               -1.7709917226851863e-07, -3.0544847853732319e-09, 5.7865420617795647e-08,  -1.0171097252531012e-07,
               -5.7257305158580533e-11, 1.2570663438515129e-06,  2.6654822896278566e-07,  -3.4050038164526772e-06,
               8.1445597585252041e-07,  -1.6635312594786876e-07, 5.1473831420242596e-08,  2.4393346653839308e-06,
               2.2357607088308811e-10,  -3.3966780568721268e-06, -1.7724752078944122e-07, 8.1474279598518956e-07,
               -3.5184071628851020e-06, -2.2561636169664630e-08, 3.0372146124014144e-07,  2.5998944339553773e-06,
               -2.7237233803075127e-10, 6.5000565283105387e-07,  -3.0565820421894761e-09, -1.6634240531075940e-07,
               -2.2531915430232931e-08, -1.3712480382138708e-06, -1.4244062693192912e-07, 1.7056091670616464e-06,
               7.3749151630759741e-12,  1.5129759351269834e-08,  5.7881730635940692e-08,  5.1460963631626650e-08,
               3.0358353727975374e-07,  -1.4243755235524716e-07, -1.1354976085189271e-06, 8.6496227926880447e-07,
               -3.4873951254771097e-11, 4.0762004111110064e-07,  3.1780002406352363e-05,  -1.2458320659895615e-04,
               2.8122504345663425e-04,  1.0436403895190930e-05,  -5.9257700343678319e-05, -1.3950884840425339e-04,
               4.7250382119895283e-08,  -6.9472394130844358e-04, -5.7259901044543953e-11, 2.2356854990096814e-10,
               -2.7227069905699890e-10, 7.3750240228374919e-12,  -3.4873905963511755e-11, -1.8469276592578567e-09,
               -7.1930737857764503e-10, 1.3498478780937949e-05,  1.0003910324077412e+00,  -2.2469869203812763e-03,
               -3.8627575895986668e-04, -3.9438593659106938e-05, -2.6693115781559678e-04, 8.4486807536367814e+02,
               3.5362734197323346e-07,  -9.5668258433216948e-03, 2.6196684856772314e-05,  9.9958942536081208e-01,
               -7.7002988984296251e-05, -2.0687524525161881e-05, 2.1076199997692550e-05,  3.7441142038009907e+01,
               1.4234818331019140e-08,  -3.6796505828102312e-04, -5.0671367828331118e-05, 1.5872218382573613e-04,
               9.9907783153857410e-01,  -2.0535848295828776e-05, 1.2118214300153330e-04,  5.5525540669911535e+01,
               -8.9924200151759468e-08, 7.7086105588059634e-04,  -5.2499064102052514e-08, -3.9263784227973904e-07,
               -1.1920617547975809e-08, 1.0000007182965609e+00,  1.6104969611651453e-06,  2.6591172234684664e-08,
               3.0061138008928027e-11,  -3.8011236242396707e-07, -2.1392541961485479e-08, 8.3035396014744179e-08,
               -2.1302795279899263e-07, -9.1041703122098208e-08, 1.0000010162173245e+00,  -3.1777384303532454e-07,
               1.4458819917101324e-11,  -1.4025390059975511e-07, 7.0666160376392904e-06,  -2.8691260724364882e-05,
               6.6287518360265706e-05,  2.2564823797330643e-06,  -1.3279828871113185e-05, 1.1274985683988947e+00,
               1.0690143341927000e-08,  -1.5601669031123345e-04, -1.4741101890795997e-05, 6.6501636667184275e-06,
               2.0657650951993110e-05,  -4.9839643622054181e-07, -1.3157312204282862e-06, -5.7764348665161232e-03,
               9.9999999857762290e-01,  1.3336991356451000e-04};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-8);
    }
    {
        // We test a case provided by Zhong Zhang tools (Tshinghua)
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1.48102e+11, 0.1, -1.70469e-10, 0.0500417, 3.49378e-18, 4.95038, 1500,
                                  0.0,         0.1, 0.1,          0.1,       0.1,         0.1,     0.0};
        std::vector<double> pars = {1.32712440018e20, 0.6, 3000.0 * 9.80665, 1e-2, 1.0};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(100.0 * 86400.0);
        std::vector<double> const ground_truth
            = {1.4803130514268405e+11,  9.9027557237587716e-02,  -1.1380120593980599e-03, 5.0011064027660777e-02,
               -2.1795658780245416e-05, 6.9815960043535004e+00,  1.4982552310395879e+03,  1.8561671535926107e-12,
               -1.6580839961454394e-01, 1.9590460144511065e-01,  1.0000530638532454e-01,  9.9993276760264940e-02,
               9.0789957994717813e-02,  -2.1460001306481999e-07, -2.6170735790319098e+20, -1.6228701132131331e+09,
               1.1082801657982795e+09,  2.6671027380313352e+07,  -1.5423690932088519e+07, 5.0310065378288114e+08,
               -7.0129895953906015e+07, 2.4195718374771244e+04,  -1.6225659115616734e+09, -1.0320064345628710e-02,
               6.6711108730635512e-03,  1.3346184564593064e-04,  2.1862835689676258e-05,  3.4915444428363892e-03,
               -9.6369401368235339e-04, 2.0843483931630457e-07,  1.1096783698155313e+09,  6.6802856780675121e-03,
               -4.9715731422605158e-03, 3.7178665742099573e-05,  1.4706130630363277e-04,  -1.8950729045978131e-03,
               -1.1268457521318763e-03, 2.1203967450856599e-07,  2.6673471640976019e+07,  1.3343675170792509e-04,
               3.7285663422010896e-05,  -5.9557777459996481e-04, 2.6491823819434440e-04,  1.5985736986089781e-04,
               -3.0416427044510564e-05, 7.9751414786776469e-09,  -1.5418486771723585e+07, 2.1951745502598216e-05,
               1.4697403423530480e-04,  2.6499630180242038e-04,  -4.2958582256626105e-04, -4.3750696535857400e-06,
               -2.1576867212768622e-05, 3.8810679523438389e-09,  -1.9864651240199022e+07, -2.5643583439404334e-04,
               4.7249092865853151e-05,  -1.7679618637666754e-05, 2.1921568347899124e-04,  5.8284576754193576e-06,
               -1.3479522866864891e-03, 1.8222190114457774e-07,  -7.0072180687067866e+07, -9.6270443996213160e-04,
               -1.1260126342321054e-03, -3.0296057452733661e-05, -2.1579926156808232e-05, -1.0428837036874235e-03,
               -1.7274654263109441e+00, 3.1834767614912022e-04,  1.0071574078993777e+00,  3.0949981232259192e-14,
               -1.6524405741410746e-14, 3.7522750092186654e-16,  -7.3314873131499181e-16, 1.8547611486588712e-11,
               3.0447170533970449e-15,  -7.6049230594778460e-19, -9.0525301975856388e+08, 9.9668061797762153e-01,
               2.1853911951405232e-03,  -1.3239467279006515e-04, 1.5916051031241970e-04,  -2.6569777954442619e+00,
               -4.7756444023765550e-04, 1.0242885383480540e-07,  2.5359522922517404e+08,  1.6193090406002280e-03,
               9.9862587617275111e-01,  -3.7495070293689949e-05, -1.2059664916215307e-05, 9.5885032999914110e-01,
               -3.9492002319824895e-06, 5.3973824208154713e-09,  3.2713235065784389e+06,  -6.3677959109269158e-06,
               -2.9744047058766255e-05, 9.9994820149972463e-01,  9.8092583617070478e-05,  4.2891312956319523e-05,
               5.2538725716624614e-06,  -9.7000827887502489e-10, 5.8536373525987770e+06,  2.9283846346380420e-05,
               2.0713649471093712e-06,  -1.3070310264186490e-04, 1.0000581351254894e+00,  -2.6037133739433378e-05,
               -6.6752193845636886e-06, 1.7502247714103865e-09,  1.0528949015429577e+08,  -3.9310186012440885e-04,
               6.3432608112206897e-04,  4.2350268245565008e-05,  -6.9378387029547624e-05, 9.0768615523697116e-01,
               3.7778449348655749e-04,  -7.7139200676634864e-08, -4.7246292918636427e+04, -6.4910517049364865e-07,
               -7.5921600488605057e-07, -2.0427170482627680e-08, -1.4550303253412321e-08, -7.0316618589944210e-07,
               9.9999978725865057e-01,  4.6470436698100695e-11};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("dynamics_time")
{
    auto ta_cached = get_ta_peq(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 3); // [mu, c1, c2]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        std::vector<double> pars = {1, 1e-4, 1};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.0);
        std::vector<double> const ground_truth
            = {9.9989683484950270e-02, 2.0000464762253184e-01, 3.0001445954786066e-01, 3.9999965108641278e-01,
               4.9999907786480108e-01, 2.6893895121143434e+01, 6.9990000000000008e-01, 4.2293276355864123e+02,
               1.9220340529988459e+01, 2.8262413554318336e+01, 5.0000094917861571e-01, 5.0000022801556965e-01,
               5.6376610930810478e-01, 4.9711715814807972e-01};

        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_time")
{
    auto ta_cached = get_ta_peq_var(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 3);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        std::vector<double> pars = {1, 1e-4, 1};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.);
        std::vector<double> const ground_truth
            = {9.9989683484950270e-02,  2.0000464762253178e-01,  3.0001445954786088e-01,  3.9999965108641289e-01,
               4.9999907786480102e-01,  2.6893895121143455e+01,  6.9989999999999997e-01,  4.2293276355864072e+02,
               1.9220340529988434e+01,  2.8262413554318286e+01,  5.0000094917861537e-01,  5.0000022801556887e-01,
               5.6376610930810933e-01,  4.9711715814808011e-01,  -4.2694045211548559e-08, 2.6649311770193512e-07,
               -1.7710270863845420e-07, -3.0564873982975259e-09, 5.7873012660844764e-08,  -1.0151288911447955e-07,
               0.0000000000000000e+00,  0.0000000000000000e+00,  2.6653990473536968e-07,  -3.4049857199962175e-06,
               8.1449685269341686e-07,  -1.6634802440339336e-07, 5.1468893034985907e-08,  2.4388280939358384e-06,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -1.7725106099626137e-07, 8.1478367379906788e-07,
               -3.5184103615851497e-06, -2.2576911510907171e-08, 3.0378511650369313e-07,  2.5996695437895554e-06,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -3.0585850539702357e-09, -1.6633730414221070e-07,
               -2.2547185370889866e-08, -1.3713184666270294e-06, -1.4242763090318087e-07, 1.7056891720972812e-06,
               0.0000000000000000e+00,  0.0000000000000000e+00,  5.7889324631465610e-08,  5.1456023985411342e-08,
               3.0364717160399241e-07,  -1.4242455595025263e-07, -1.1355471101040182e-06, 8.6497914583340125e-07,
               0.0000000000000000e+00,  0.0000000000000000e+00,  3.1778870542065043e-05,  -1.2458002867030360e-04,
               2.8122628721927932e-04,  1.0438557933502884e-05,  -5.9266343096809520e-05, -1.3959734392773420e-04,
               0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
               0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
               0.0000000000000000e+00,  0.0000000000000000e+00,  1.0003910173991175e+00,  -2.2470440567257570e-03,
               -3.8633062259343964e-04, -3.9435128950313492e-05, -2.6696600612478231e-04, 8.4486807587569695e+02,
               0.0000000000000000e+00,  0.0000000000000000e+00,  2.6195808702406448e-05,  9.9958942252195493e-01,
               -7.7002001872287745e-05, -2.0687852639967772e-05, 2.1078252277565067e-05,  3.7441142053248441e+01,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -5.0670603442863455e-05, 1.5871640599372683e-04,
               9.9907782449614047e-01,  -2.0540856128558779e-05, 1.2120175807010855e-04,  5.5525540577435912e+01,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -5.2505497434348618e-08, -3.9265208069234857e-07,
               -1.1963391937008132e-08, 1.0000007183054169e+00,  1.6105620415885631e-06,  2.6610743446749823e-08,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -2.1395803096787098e-08, 8.3039238628006595e-08,
               -2.1305527926752284e-07, -9.1065218453256924e-08, 1.0000010162674751e+00,  -3.1775927366261674e-07,
               0.0000000000000000e+00,  0.0000000000000000e+00,  7.0663704326338798e-06,  -2.8690761529955843e-05,
               6.6287695848795967e-05,  2.2569625301650765e-06,  -1.3281790498939194e-05, 1.1274985801394359e+00,
               0.0000000000000000e+00,  0.0000000000000000e+00,  -1.4741142302228122e-05, 6.6503221592830074e-06,
               2.0657458995708366e-05,  -4.9839125053452309e-07, -1.3157556034061730e-06, -5.7764361958391889e-03,
               1.0000000000000000e+00,  0.0000000000000000e+00};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        // This test is for highly unphysical conditions and the numerics of high derivatives becomes unstable
        // to the point of losing 8 digits of precision on different machines. So we ask 1e-8
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-8);
    }
}

TEST_CASE("various_cfunc")
{
    REQUIRE_NOTHROW(kep3::ta::get_peq_dyn_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_peq_dyn_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_peq_i_vers_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_peq_i_vers_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_peq_u_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_peq_u_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_peq_SF_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_peq_SF_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_peq_H_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_peq_H_cfunc(kep3::optimality_type::MASS));
}