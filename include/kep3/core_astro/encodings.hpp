// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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