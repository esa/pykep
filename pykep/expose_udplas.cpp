// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <string>

#include <pybind11/pybind11.h>

#include <kep3/planet.hpp>
#include <kep3/udpla/jpl_lp.hpp>
#include <kep3/udpla/keplerian.hpp>
#include <kep3/udpla/vsop2013.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "expose_udplas.hpp"

namespace pykep
{
namespace py = pybind11;

// Split algorithm exposition functions.
void expose_all_udplas(py::module &udpla_module, py::class_<kep3::planet> &planet_class) // NOLINT
{
    // null udpla
    auto null_udpla = pykep::expose_one_udpla<kep3::detail::null_udpla>(
        udpla_module, planet_class, "_null_udpla", "A moot udpla used as default to construct a planet.");
    // Constructor.
    null_udpla.def(py::init<>());

    // keplerian udpla
    auto keplerian_udpla
        = pykep::expose_one_udpla<kep3::udpla::keplerian>(udpla_module, planet_class, "_keplerian", "keplerian udpla");
    // Constructors.
    keplerian_udpla
        .def(py::init<const kep3::epoch &, const std::array<double, 6> &, double, std::string, std::array<double, 3>,
                      kep3::elements_type>(),
             py::arg("when"), py::arg("elem"), py::arg("mu_central_body"), py::arg("name") = "unknown UDPLA",
             py::arg("added_params") = std::array<double, 3>({-1, -1, -1}),
             py::arg("el_type") = kep3::elements_type::KEP_F, pykep::udpla_keplerian_from_elem_docstring().c_str())
        .def(py::init<const kep3::epoch &, const std::array<std::array<double, 3>, 2> &, double, std::string,
                      std::array<double, 3>>(),
             py::arg("when"), py::arg("posvel"), py::arg("mu_central_body"), py::arg("name") = "unknown UDPLA",
             py::arg("added_params") = std::array<double, 3>({-1, -1, -1}),
             pykep::udpla_keplerian_from_posvel_docstring().c_str())
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::udpla::keplerian>)
        // other methods
        .def_property_readonly("ref_epoch", &kep3::udpla::keplerian::get_ref_epoch);

    // jpl_lp_udpla udpla
    auto jpl_lp_udpla = pykep::expose_one_udpla<kep3::udpla::jpl_lp>(
        udpla_module, planet_class, "_jpl_lp", "Low precision ephemerides from JPL for the solar system planets");
    // Constructors.
    jpl_lp_udpla
        .def(py::init<std::string>(), py::arg("body") = "earth", pykep::udpla_jpl_lp_docstring().c_str())
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::udpla::jpl_lp>)
        .def_property("safe_radius", &kep3::udpla::jpl_lp::get_safe_radius, &kep3::udpla::jpl_lp::set_safe_radius, "The planet's safe radius. Mostly used to constraint fly-by distances.")
        .def("set_safe_radius", &kep3::udpla::jpl_lp::set_safe_radius, "Set the planet's safe radius. Mostly used to constraint fly-by distances.")
        .def("get_radius", &kep3::udpla::jpl_lp::get_radius, "The planet's radius.")
        .def_property_readonly("radius", &kep3::udpla::jpl_lp::get_radius, "The planet's radius.")
        .def("get_radius", &kep3::udpla::jpl_lp::get_radius, "The planet's radius.");
        
    // vsop2013 udpla.
    auto vsop2013_udpla = pykep::expose_one_udpla<kep3::udpla::vsop2013>(
        udpla_module, planet_class, "_vsop2013", "Analytical planetary ephemerides from the VSOP2013 theory");
    // Constructors.
    vsop2013_udpla.def(py::init<std::string, double>(), py::arg("body") = "mercury", py::arg("thresh") = 1e-5,
                       pykep::udpla_vsop2013_docstring().c_str());
}

} // namespace pykep
