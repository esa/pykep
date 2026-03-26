// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef PYKEP_PLANET_HPP
#define PYKEP_PLANET_HPP

#include <array>
#include <string>

#include <kep3/core_astro/constants.hpp>
#include <kep3/epoch.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pykep
{
namespace py = pybind11;
struct python_udpla {
    py::object m_obj;

    python_udpla();
    explicit python_udpla(py::object obj);

    // Mandatory methods
    [[nodiscard]] std::array<std::array<double, 3>, 2> eph(double) const;

    // Optional methods
    [[nodiscard]] std::vector<double> eph_v(const std::vector<double> &) const;
    [[nodiscard]] std::array<double, 3> acc(double) const;
    [[nodiscard]] std::vector<double> acc_v(const std::vector<double> &) const;
    [[nodiscard]] std::string get_name() const;
    [[nodiscard]] std::string get_extra_info() const;
    [[nodiscard]] double get_mu_central_body() const;
    [[nodiscard]] double get_mu_self() const;
    [[nodiscard]] double get_radius() const;
    [[nodiscard]] double get_safe_radius() const;
    [[nodiscard]] double period(double) const;
    [[nodiscard]] std::array<double, 6> elements(double, kep3::elements_type el_type = kep3::elements_type::KEP_F) const;
};
} // namespace pykep

#endif