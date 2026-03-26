// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/mima.hpp>
#include <kep3/lambert_problem.hpp>

#include "catch.hpp"
#include "kep3/udpla/keplerian.hpp"

TEST_CASE("mima")
{
    // We take the first item from the zeonodo database https://zenodo.org/records/11502524 containing
    // independent cases reporting the mima values.
    {
        std::array<double, 3> rs = {3.574644002632926178e+10, -5.688222150272903442e+10, -1.304897435568400574e+10};
        std::array<double, 3> vs = {4.666425901145393436e+04, 2.375697019573154466e+04, 1.165422004315219965e+04};
        std::array<double, 3> rt = {7.672399994418635559e+10, -1.093562401274179382e+11, 4.796635567684053421e+09};
        std::array<double, 3> vt = {2.725105661271001009e+04, 1.599599495457483499e+04, 6.818440757625087826e+03};
        double tof = 3.311380772794449854e+02 * kep3::DAY2SEC;
        double Tmax = 0.6;
        double veff = kep3::G0 * 4000;
        auto lp = kep3::lambert_problem(rs, rt, tof, kep3::MU_SUN);
        std::array<double, 3> dv1 = {lp.get_v0()[0][0] - vs[0], lp.get_v0()[0][1] - vs[1], lp.get_v0()[0][2] - vs[2]};
        std::array<double, 3> dv2 = {vt[0] - lp.get_v1()[0][0], vt[1] - lp.get_v1()[0][1], vt[2] - lp.get_v1()[0][2]};
        auto mima_res = kep3::mima(dv1, dv2, tof, Tmax, veff);
        double mima_from_zenodo_db = 1.711341975126993020e+02;
        REQUIRE(mima_res.first == Approx(mima_from_zenodo_db).epsilon(1e-8));
    }
    // We take a second item from the zeonodo database https://zenodo.org/records/11502524
    {
        std::array<double, 3> rs = {-4.054163103119991455e+11, -3.112051036509102173e+11, 8.852823556219964600e+10};
        std::array<double, 3> vs = {1.053626627938712227e+04, -9.040187399709659076e+03, 6.140398196916326924e+03};
        std::array<double, 3> rt = {-4.377345423917691040e+10, -4.913367837977642822e+11, 3.030052465928871918e+10};
        std::array<double, 3> vt = {1.549115173657366176e+04, -3.341220883615214916e+02, -3.245198147308494299e+03};
        double tof = 3.219932820383384069e+02 * kep3::DAY2SEC;
        double Tmax = 0.6;
        double veff = kep3::G0 * 4000;
        auto lp = kep3::lambert_problem(rs, rt, tof, kep3::MU_SUN);
        std::array<double, 3> dv1 = {lp.get_v0()[0][0] - vs[0], lp.get_v0()[0][1] - vs[1], lp.get_v0()[0][2] - vs[2]};
        std::array<double, 3> dv2 = {vt[0] - lp.get_v1()[0][0], vt[1] - lp.get_v1()[0][1], vt[2] - lp.get_v1()[0][2]};
        auto mima_res = kep3::mima(dv1, dv2, tof, Tmax, veff);
        double mima_from_zenodo_db = 1.101217159182178875e+03;
        REQUIRE(mima_res.first == Approx(mima_from_zenodo_db).epsilon(1e-8));
    }
}

TEST_CASE("mima_from_hop")
{
    kep3::epoch when{64328.0, kep3::epoch::julian_type::MJD};
    std::array<double, 6> el_s
        = {286031128778.39996, 0.0763, 0.3153809958353754, 4.70313873534912, 5.27246513734967, 2.9460074243530188};
    std::array<double, 6> el_f
        = {270472950225.6, 0.066, 0.39968039870670147, 4.255287249287375, 3.9732420421650914, 2.3813298840517};
    kep3::planet pl_s{kep3::udpla::keplerian(when, el_s, kep3::MU_SUN)};
    kep3::planet pl_f{kep3::udpla::keplerian(when, el_f, kep3::MU_SUN)};
    double Tmax = 0.6;
    double veff = kep3::G0 * 4000;
    kep3::epoch when_s{5500.0, kep3::epoch::julian_type::MJD2000};
    kep3::epoch when_f{5700.0, kep3::epoch::julian_type::MJD2000};
    auto mima_res = kep3::mima_from_hop(pl_s, pl_f, when_s, when_f, Tmax, veff);
    double ground_truth = 1270.16102850;
    REQUIRE(mima_res.first == Approx(ground_truth).epsilon(1e-8));
}

TEST_CASE("mima2")
{
    // We take the first item from the zeonodo database https://zenodo.org/records/11502524 containing
    // independent cases reporting the mima2 values.
    {
        std::array<double, 3> rs = {3.574644002632926178e+10, -5.688222150272903442e+10, -1.304897435568400574e+10};
        std::array<double, 3> vs = {4.666425901145393436e+04, 2.375697019573154466e+04, 1.165422004315219965e+04};
        std::array<double, 3> rt = {7.672399994418635559e+10, -1.093562401274179382e+11, 4.796635567684053421e+09};
        std::array<double, 3> vt = {2.725105661271001009e+04, 1.599599495457483499e+04, 6.818440757625087826e+03};
        double tof = 3.311380772794449854e+02 * kep3::DAY2SEC;
        double Tmax = 0.6;
        double veff = kep3::G0 * 4000.;
        auto lp = kep3::lambert_problem(rs, rt, tof, kep3::MU_SUN);
        std::array<double, 3> dv1 = {lp.get_v0()[0][0] - vs[0], lp.get_v0()[0][1] - vs[1], lp.get_v0()[0][2] - vs[2]};
        std::array<double, 3> dv2 = {vt[0] - lp.get_v1()[0][0], vt[1] - lp.get_v1()[0][1], vt[2] - lp.get_v1()[0][2]};
        auto mima2_res = kep3::mima2({rs, lp.get_v0()[0]}, dv1, dv2, tof, Tmax, veff, kep3::MU_SUN);
        double mima2_from_zenodo_db = 1.397851641912264995e+02;
        REQUIRE(mima2_res.first == Approx(mima2_from_zenodo_db).epsilon(1e-8));
    }
    { // This second case is also from the Zenodo db
        std::array<double, 3> rs = {8.899464427764886475e+10, -4.581927411496286621e+11, 2.048886307096130981e+11};
        std::array<double, 3> vs = {1.385571222713435418e+04, 6.481857970028194359e+03, 2.812533151527441078e+03};
        std::array<double, 3> rt = {3.668961862639051514e+11, -2.093798042740150452e+11, 4.217417200463520050e+10};
        std::array<double, 3> vt = {5.115925672273307100e+03, 1.436244285179517283e+04, -7.297450556937028523e+03};
        double tof = 3.024673748374755746e+02 * kep3::DAY2SEC;
        double Tmax = 0.6;
        double veff = kep3::G0 * 4000;
        auto lp = kep3::lambert_problem(rs, rt, tof, kep3::MU_SUN);
        std::array<double, 3> dv1 = {lp.get_v0()[0][0] - vs[0], lp.get_v0()[0][1] - vs[1], lp.get_v0()[0][2] - vs[2]};
        std::array<double, 3> dv2 = {vt[0] - lp.get_v1()[0][0], vt[1] - lp.get_v1()[0][1], vt[2] - lp.get_v1()[0][2]};
        vs[0] += dv1[0];
        vs[1] += dv1[1];
        vs[2] += dv1[2];
        auto mima2_res = kep3::mima2({rs, lp.get_v0()[0]}, dv1, dv2, tof, Tmax, veff, kep3::MU_SUN);
        double mima2_from_zenodo_db = 1092.1862621801;
        REQUIRE(mima2_res.first == Approx(mima2_from_zenodo_db).epsilon(1e-8));
    }
}

TEST_CASE("mima2_from_hop")
{
    kep3::epoch when{64328.0, kep3::epoch::julian_type::MJD};
    std::array<double, 6> el_s
        = {286031128778.39996, 0.0763, 0.3153809958353754, 4.70313873534912, 5.27246513734967, 2.9460074243530188};
    std::array<double, 6> el_f
        = {270472950225.6, 0.066, 0.39968039870670147, 4.255287249287375, 3.9732420421650914, 2.3813298840517};
    kep3::planet pl_s{kep3::udpla::keplerian(when, el_s, kep3::MU_SUN)};
    kep3::planet pl_f{kep3::udpla::keplerian(when, el_f, kep3::MU_SUN)};
    double Tmax = 0.6;
    double veff = kep3::G0 * 4000;
    kep3::epoch when_s{5500.0, kep3::epoch::julian_type::MJD2000};
    kep3::epoch when_f{5700.0, kep3::epoch::julian_type::MJD2000};
    auto mima2_res = kep3::mima2_from_hop(pl_s, pl_f, when_s, when_f, Tmax, veff);
    double ground_truth = 1336.53752329;
    REQUIRE(mima2_res.first == Approx(ground_truth).epsilon(1e-8));
}