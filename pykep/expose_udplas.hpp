// Copyright (c) 2023-2026 Dario Izzo (dario.izzo@gmail.com)
//                          Advanced Concepts Team, European Space Agency (ESA)
//
// This file is part of the kep3 library.
//
// SPDX-License-Identifier: MPL-2.0
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef PYKEP_EXPOSE_UDPLAS_HPP
#define PYKEP_EXPOSE_UDPLAS_HPP

#include <pybind11/pybind11.h>

#include <kep3/planet.hpp>

#include "common_utils.hpp"

namespace pykep
{
namespace py = pybind11;

// C++ UDPLA exposition function.
template <typename Udpla>
inline py::class_<Udpla> expose_one_udpla(py::module &p_module, py::class_<kep3::planet> &pla, // NOLINT
                                          const char *name, const char *descr)
{
    py::class_<Udpla> c(p_module, name, descr);

    // We require all udplas to be def-ctible at the bare minimum.
    c.def(py::init<>());

    // We mark it as a C++ udpla.
    c.attr("_pykep_cpp_udpla") = true;

    // Expose the udpa constructor from planet.
    pla.def(py::init<const Udpla &>(), py::arg("udpla"));

    // Expose extract of the planet class.
    pla.def("_cpp_extract", &generic_cpp_extract<kep3::planet, Udpla>, py::return_value_policy::reference_internal);

    // Add the udpla to the udpla submodule.
    p_module.attr(name) = c;

    return c;
}

void expose_all_udplas(py::module &, py::class_<kep3::planet> &);

} // namespace pykep

#endif