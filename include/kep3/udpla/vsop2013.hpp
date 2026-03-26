// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef kep3_UDPLA_VSOP2013_H
#define kep3_UDPLA_VSOP2013_H

#include <array>
#include <memory>
#include <string>

#include <kep3/detail/s11n.hpp>
#include <kep3/detail/visibility.hpp>
#include <kep3/planet.hpp>

namespace kep3::udpla
{

class kep3_DLL_PUBLIC vsop2013
{
    struct impl;

    std::unique_ptr<impl> m_impl;

    friend class boost::serialization::access;
    void save(boost::archive::binary_oarchive &, unsigned) const;
    void load(boost::archive::binary_iarchive &, unsigned);
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    vsop2013();
    explicit vsop2013(std::string, double = 1e-5);
    vsop2013(vsop2013 &&) noexcept;
    vsop2013(const vsop2013 &);
    vsop2013 &operator=(vsop2013 &&) noexcept;
    vsop2013 &operator=(const vsop2013 &);
    ~vsop2013();

    // Mandatory UDPLA methods.
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double) const;

    // Optional UDPLA methods
    [[nodiscard]] std::string get_name() const;
};

} // namespace kep3::udpla

KEP3_S11N_EXPORT_KEY_AND_EXTERN_TEMPLATES(kep3::udpla::vsop2013, kep3::detail::planet_iface)

#endif
