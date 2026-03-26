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
#include <kep3/ta/pontryagin_cartesian.hpp>

#include "catch.hpp"
#include "test_helpers.hpp"

using heyoka::taylor_adaptive;
using heyoka::taylor_outcome;

using kep3::ta::get_ta_pc;
using kep3::ta::get_ta_pc_cache_dim;
using kep3::ta::get_ta_pc_var;
using kep3::ta::get_ta_pc_var_cache_dim;

// This needs to be the first test as cache dimension will be assumed to be zero here.
TEST_CASE("caches")
{
    {
        // MASS
        //  The non variational one.
        REQUIRE(get_ta_pc_cache_dim() == 0u);
        auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 1u);
        ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 1u);
        ta_cached = get_ta_pc(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_cache_dim() == 2u);

        // The variational integrator.
        REQUIRE(get_ta_pc_var_cache_dim() == 0u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 1u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 1u);
        ta_cached = get_ta_pc_var(1e-8, kep3::optimality_type::MASS);
        REQUIRE(get_ta_pc_var_cache_dim() == 2u);
    }
    {
        // TIME
        //  The non variational one. (no cache is not empty)
        REQUIRE(get_ta_pc_cache_dim() == 2u);
        auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 3u);
        ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 3u);
        ta_cached = get_ta_pc(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_cache_dim() == 4u);

        // The variational integrator.
        REQUIRE(get_ta_pc_var_cache_dim() == 2u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 3u);
        ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 3u);
        ta_cached = get_ta_pc_var(1e-8, kep3::optimality_type::TIME);
        REQUIRE(get_ta_pc_var_cache_dim() == 4u);
    }
}

TEST_CASE("dynamics_mass")
{
    auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 5); // [mu, c1, c2, eps, l0]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1., 0.5, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {0.3296405122183833,  0.9437361993515112, -0.000126441649019,  -0.9446333908544454, 0.3295140980732467,
               -0.0000302961826291, 9.993495810584642,  -0.2930516116019541, 0.1612382337718589,  1.2739955068582864,
               0.9472048639975497,  0.0186280856006421, -0.6140454396473527, 0.9999296119431648};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_mass")
{
    auto ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::MASS);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 5);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1., 0.5, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {0.3296405122183832,  0.9437361993515109,  -0.000126441649019,  -0.9446333908544454, 0.3295140980732463,
               -0.0000302961826292, 9.99349581058464,    -0.2930516116019541, 0.1612382337718588,  1.273995506858286,
               0.9472048639975497,  0.0186280856006422,  -0.6140454396473525, 0.9999296119431648,  0.0000510654227239,
               -0.0000298189349962, -0.0000217842414248, -0.0001709127492249, 0.0000606364462851,  0.0000914445088186,
               -0.0001647862264403, 0.0001841557742587,  -0.0000292258835123, 0.0001409721990404,  -0.0000221621982323,
               0.0000615571481262,  -0.0002458392872784, 0.0000801765457717,  -0.0001148195469165, 0.0001293410230013,
               -0.0000131533128852, -0.0000133007447032, 0.000128807585954,   0.0000627595124607,  0.0000491882508609,
               -0.0002220482039032, -0.0000585156386462, 0.0000662625508623,  0.0000952627093301,  -0.0000136994934963,
               0.0000101498938056,  -0.0002217814547051, -0.0000067094844097, 0.0001031363758976,  -0.0003052447368777,
               0.0003388861904555,  -0.0000103695140484, 0.0003988136559688,  -0.0000397571410607, -0.0000344582704363,
               -0.0004852344245896, 0.0001438034675844,  -0.0002232837605096, 0.0002504859870914,  0.0000564539416515,
               0.0000039580624579,  0.0002869315036968,  -0.0000528897225851, 0.0000099634335878,  -0.0003066917694243,
               -0.0000136899077258, 0.0000159644583412,  0.0001984065648323,  0.0000662507174415,  -0.0000063482502781,
               -0.0003813109734398, -0.0001467019621364, -0.0000597339367116, -0.0030581645051647, 0.0033876023454567,
               1.4420740730782136,  0.920504461192634,   -0.0000866599738527, -1.5430230109206016, -1.1126263439809247,
               0.0001026034647977,  -0.0000121654716071, 0.0000154310093872,  1.2656777129337764,  1.6024311541672625,
               -0.0001601302495506, -2.3852206033936696, -0.321767110272336,  0.0002446240120997,  -0.0002483963509898,
               0.0002809829252663,  -0.0001377758644031, -0.0001171273685256, 0.3295958024617905,  0.0002174314434472,
               0.0000497860337846,  0.9443892336975448,  0.0000083315281977,  -0.0000101750735497, -1.5764972150288585,
               -0.3739582379234117, 0.0000404235543384,  2.265054465978205,   0.6326590269354981,  -0.0000879455572974,
               0.0000383485449313,  -0.0000440025058561, -0.4490943638822759, -1.2889700541945686, 0.000046069423557,
               0.9777263320447432,  0.7790169660888409,  -0.0000881434138879, 0.0000662209863897,  -0.0000749414521565,
               0.000040643607706,   0.0000309397180723,  -0.9438792979197026, -0.0000718529150407, -0.0000198552061604,
               0.3298502486562563,  -0.0000253407256842, 0.0000290751372009,  0.0000437812740011,  0.0000146715804732,
               -0.000001352215403,  -0.0000845890189984, -0.0000329083076638, -0.0000137446167106, 0.999967014454261,
               0.0000367387932054};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("dynamics_time")
{
    auto ta_cached = get_ta_pc(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == false);
    REQUIRE(ta_cached.get_dim() == 14);
    REQUIRE(ta_cached.get_pars().size() == 3); // [mu, c1, c2]

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.29326772571214e-01,  9.43518279694245e-01,  -2.37243691077672e-04, -9.45216052816242e-01,
               3.29089569819722e-01,  -5.58634286200389e-05, 9.98765500000000e+00,  -2.93073464383497e-01,
               1.60767429870367e-01,  1.27401079872040e+00,  9.47277110771498e-01,  1.87535673100874e-02,
               -6.14093196354726e-01, 9.99866710251057e-01};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-13);
    }
}

TEST_CASE("variational_dynamics_time")
{
    auto ta_cached = get_ta_pc_var(1e-16, kep3::optimality_type::TIME);
    REQUIRE(ta_cached.is_variational() == true);
    REQUIRE(ta_cached.get_dim() == 14 + 14 * 8);
    REQUIRE(ta_cached.get_pars().size() == 3);
    REQUIRE(ta_cached.get_vorder() == 1);
    REQUIRE(ta_cached.get_vargs().size() == 8);

    {
        // We test a meaningless case.
        taylor_adaptive<double> ta(ta_cached); // making a copy as to be able to modify the object.
        std::vector<double> ic = {1., 0., 0., 0., 1., 0., 10., 1., 1., 1., 1., 1., 1., 1.};
        std::vector<double> pars = {1., 0.01, 1.};
        ta.set_time(0.);
        std::copy(ic.begin(), ic.end(), ta.get_state_data());
        std::copy(pars.begin(), pars.end(), ta.get_pars_data());

        auto out = ta.propagate_until(1.2345);
        std::vector<double> const ground_truth
            = {3.29326772571214e-01,  9.43518279694245e-01,  -2.37243691077678e-04, -9.45216052816242e-01,
               3.29089569819723e-01,  -5.58634286200446e-05, 9.98765500000000e+00,  -2.93073464383497e-01,
               1.60767429870367e-01,  1.27401079872040e+00,  9.47277110771498e-01,  1.87535673100874e-02,
               -6.14093196354726e-01, 9.99866710251056e-01,  8.41133971750316e-05,  -6.23805339365836e-05,
               -4.35576346536962e-05, -2.91230168218841e-04, 1.30236152736551e-04,  1.82818786897539e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -6.15912904833284e-05, 2.66032245545687e-04,
               -4.43225236984674e-05, 1.34507269522455e-04,  -4.54930989076605e-04, 1.60305288190260e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -2.62904019617088e-05, -2.65922554162050e-05,
               2.44305289111464e-04,  1.25451303614909e-04,  9.83231934912956e-05,  -4.15197128839756e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.45747427272735e-04,  -3.84863004523546e-05,
               2.03218553225493e-05,  -3.49785695279272e-04, 1.60038710912614e-05,  2.06198842045081e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -3.42763969439676e-05, 7.55063066434171e-04,
               -7.95411628295219e-05, -2.67105446050315e-05, -9.02168112273193e-04, 2.87633150217543e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.12967971828034e-04,  7.93937073268113e-06,
               5.44943821249505e-04,  -1.05916110655104e-04, 1.98653792978839e-05,  -5.79800432453000e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,
               0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,  0.00000000000000e+00,
               0.00000000000000e+00,  0.00000000000000e+00,  1.44196561049032e+00,  9.20877930716147e-01,
               -1.62136878165076e-04, -1.54295185938881e+00, -1.11299681577074e+00, 1.93806447760147e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  1.26630644933279e+00,  1.60313599146158e+00,
               -3.05364032176103e-04, -2.38639680762912e+00, -3.22452855839271e-01, 4.80016576560558e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -2.59565205932807e-04, -2.18973016400636e-04,
               3.29244288341037e-01,  4.08237064294930e-04,  9.12183025923601e-05,  9.44745593234806e-01,
               0.00000000000000e+00,  0.00000000000000e+00,  -1.57652935813679e+00, -3.74079622523328e-01,
               7.63176133867316e-05,  2.26516174045298e+00,  6.32816141123617e-01,  -1.68107758370038e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  -4.49254170142732e-01, -1.28907236682424e+00,
               8.78481189862216e-05,  9.78050639895916e-01,  7.79112671238256e-01,  -1.71054976102050e-04,
               0.00000000000000e+00,  0.00000000000000e+00,  7.56303455884195e-05,  5.74555043756629e-05,
               -9.43788098497564e-01, -1.31402239732434e-04, -3.44548398801908e-05, 3.29727673372486e-01,
               0.00000000000000e+00,  0.00000000000000e+00,  8.00334644311539e-05,  2.67404597714927e-05,
               -2.55263719046308e-06, -1.53938341730639e-04, -5.93271855549472e-05, -2.42455086700587e-05,
               1.00000000000000e+00,  0.00000000000000e+00};
        REQUIRE(std::get<0>(out) == taylor_outcome::time_limit);
        REQUIRE(kep3_tests::L_infinity_norm_rel(ta.get_state(), ground_truth) <= 1e-12);
    }
}

TEST_CASE("various_cfunc")
{
    REQUIRE_NOTHROW(kep3::ta::get_pc_dyn_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_pc_dyn_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_pc_i_vers_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_pc_i_vers_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_pc_u_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_pc_u_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_pc_SF_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_pc_SF_cfunc(kep3::optimality_type::MASS));
    REQUIRE_NOTHROW(kep3::ta::get_pc_H_cfunc(kep3::optimality_type::TIME));
    REQUIRE_NOTHROW(kep3::ta::get_pc_H_cfunc(kep3::optimality_type::MASS));
}