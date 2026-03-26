// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_mee2PAR2mee_H
#define kep3_mee2PAR2mee_H

#include <array>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::array<double, 6> mee2par(const std::array<double, 6> &eq, bool retrogade = false);

kep3_DLL_PUBLIC std::array<double, 6> par2mee(const std::array<double, 6> &par, bool retrogade = false);

} // namespace kep3
#endif // kep3_mee2PAR2mee_H