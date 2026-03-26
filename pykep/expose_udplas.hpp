// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

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