// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xadapt.hpp>

#include <kep3/linalg.hpp>


namespace kep3::linalg
{
mat33 _skew(const mat31 &v)
{
    return {{0., -v(2, 0), v(1, 0)}, {v(2, 0), 0., -v(0, 0)}, {-v(1, 0), v(0, 0), 0.}};
}

mat31 _cross(const mat31 &v1, const mat31 &v2)
{
    return {{v1(1, 0) * v2(2, 0) - v2(1, 0) * v1(2, 0), v2(0, 0) * v1(2, 0) - v1(0, 0) * v2(2, 0),
             v1(0, 0) * v2(1, 0) - v2(0, 0) * v1(1, 0)}};
}
} // namespace kep3::linalg
