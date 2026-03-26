// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_IC2mee2IC_H
#define kep3_IC2mee2IC_H

#include <array>

#include <kep3/detail/visibility.hpp>

namespace kep3
{

kep3_DLL_PUBLIC std::array<double, 6> ic2mee(const std::array<std::array<double, 3>, 2> &pos_vel, double mu,
                                            bool retrogade = false);

kep3_DLL_PUBLIC std::array<std::array<double, 3>, 2> mee2ic(const std::array<double, 6> &eq, double mu,
                                                           bool retrogade = false);

} // namespace kep3
#endif // kep3_IC2mee2IC_H