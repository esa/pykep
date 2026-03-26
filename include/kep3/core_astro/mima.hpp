// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_MIMA_H
#define kep3_MIMA_H

#include "kep3/epoch.hpp"
#include "kep3/planet.hpp"
#include <array>
#include <utility>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

// mima (https://ieeexplore.ieee.org/abstract/document/7850107)
// Hennes, D., Izzo, D., & Landau, D. (2016, December). Fast approximators for optimal low-thrust hops between main belt
// asteroids. In 2016 IEEE Symposium Series on Computational Intelligence (SSCI) (pp. 1-7). IEEE.
kep3_DLL_PUBLIC std::pair<double, double> mima(const std::array<double, 3> &dv1, const std::array<double, 3> &dv2,
                                               double tof, double Tmax, double veff);

kep3_DLL_PUBLIC std::pair<double, double> mima_from_hop(const kep3::planet &pl_s, const kep3::planet &pl_f,
                                                        const kep3::epoch &when_s, const kep3::epoch &when_f,
                                                        double Tmax, double veff);

// mima2 (https://arxiv.org/pdf/2410.20839)
// Izzo, D., ... & Yam, C. H. (2025). Asteroid mining: ACT&Friends’ results for the GTOC12 problem. Astrodynamics, 9(1),
// 19-40.
kep3_DLL_PUBLIC std::pair<double, double> mima2(const std::array<std::array<double, 3>, 2> &posvel1,
                                                const std::array<double, 3> &dv1, const std::array<double, 3> &dv2,
                                                double tof, double Tmax, double veff, double mu);

kep3_DLL_PUBLIC std::pair<double, double> mima2_from_hop(const kep3::planet &pl_s, const kep3::planet &pl_f,
                                                        const kep3::epoch &when_s, const kep3::epoch &when_f,
                                                        double Tmax, double veff);

std::pair<double, double> _mima_compute_transfer(double x, const std::array<std::array<double, 3>, 2> &posvel1,
                                                         double tof, const std::array<double, 3> &dv1,
                                                         const std::array<double, 3> &dv2, double mu);

} // namespace kep3
#endif // kep3_MIMA_H