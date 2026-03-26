// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef kep3_ENCODING_CONVERSIONS_H
#define kep3_ENCODING_CONVERSIONS_H

#include <vector>

#include <kep3/detail/visibility.hpp>

namespace kep3
{
kep3_DLL_PUBLIC std::vector<double> alpha2direct(const std::vector<double> &alphas, double tof);
kep3_DLL_PUBLIC std::pair<std::vector<double>, double> direct2alpha(const std::vector<double> &tofs);
kep3_DLL_PUBLIC std::vector<double> eta2direct(const std::vector<double> &etas, double max_tof);
kep3_DLL_PUBLIC std::vector<double> direct2eta(const std::vector<double> &tofs, double max_tof);
} // namespace kep3

#endif