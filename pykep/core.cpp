// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.
#include <pybind11/pytypes.h>
#include <string>

#include <fmt/chrono.h>

#include <kep3/core_astro/basic_transfers.hpp>
#include <kep3/core_astro/constants.hpp>
#include <kep3/core_astro/convert_anomalies.hpp>
#include <kep3/core_astro/encodings.hpp>
#include <kep3/core_astro/flyby.hpp>
#include <kep3/core_astro/ic2mee2ic.hpp>
#include <kep3/core_astro/ic2par2ic.hpp>
#include <kep3/core_astro/mee2par2mee.hpp>
#include <kep3/core_astro/mima.hpp>
#include <kep3/core_astro/propagate_lagrangian.hpp>
#include <kep3/epoch.hpp>
#include <kep3/lambert_problem.hpp>
#include <kep3/leg/sims_flanagan.hpp>
#include <kep3/leg/sims_flanagan_alpha.hpp>
#include <kep3/planet.hpp>
#include <kep3/ta/bcp.hpp>
#include <kep3/ta/cr3bp.hpp>
#include <kep3/ta/kep.hpp>
#include <kep3/ta/pontryagin_cartesian.hpp>
#include <kep3/ta/pontryagin_equinoctial.hpp>
#include <kep3/ta/zero_hold_cr3bp.hpp>
#include <kep3/ta/zero_hold_kep.hpp>
#include <kep3/udpla/jpl_lp.hpp>
#include <kep3/udpla/keplerian.hpp>
#include <kep3/zero_hold_kep_problem.hpp>

#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "expose_udplas.hpp"

#include "python_udpla.hpp"

namespace py = pybind11;
namespace pk = pykep;

PYBIND11_MODULE(core, m) // NOLINT
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = pk::core_module_doc();

    // We expose various global constants:
    m.attr("AU") = py::float_(kep3::AU);
    m.attr("CAVENDISH") = py::float_(kep3::CAVENDISH);
    m.attr("MU_SUN") = py::float_(kep3::MU_SUN);
    m.attr("MU_EARTH") = py::float_(kep3::MU_EARTH);
    m.attr("MU_MOON") = py::float_(kep3::MU_MOON);
    m.attr("EARTH_VELOCITY") = py::float_(kep3::EARTH_VELOCITY);
    m.attr("EARTH_J2") = py::float_(kep3::EARTH_J2);
    m.attr("EARTH_RADIUS") = py::float_(kep3::EARTH_RADIUS);
    m.attr("RAD2DEG") = py::float_(kep3::RAD2DEG);
    m.attr("DEG2RAD") = py::float_(kep3::DEG2RAD);
    m.attr("DAY2SEC") = py::float_(kep3::DAY2SEC);
    m.attr("SEC2DAY") = py::float_(kep3::SEC2DAY);
    m.attr("YEAR2DAY") = py::float_(kep3::YEAR2DAY);
    m.attr("DAY2YEAR") = py::float_(kep3::DAY2YEAR);
    m.attr("G0") = py::float_(kep3::G0);
    m.attr("CR3BP_MU_EARTH_MOON") = py::float_(kep3::CR3BP_MU_EARTH_MOON);
    m.attr("BCP_MU_EARTH_MOON") = py::float_(kep3::BCP_MU_EARTH_MOON);
    m.attr("BCP_MU_S") = py::float_(kep3::BCP_MU_S);
    m.attr("BCP_RHO_S") = py::float_(kep3::BCP_RHO_S);
    m.attr("BCP_OMEGA_S") = py::float_(kep3::BCP_OMEGA_S);

    // We expose here global enums:
    py::enum_<kep3::elements_type>(m, "el_type", "")
        .value("KEP_M", kep3::KEP_M, "Keplerian Elements :math:`[a,e,i,\\Omega,\\omega,M]` (Mean anomaly)")
        .value("KEP_F", kep3::KEP_F, "Keplerian Elements :math:`[a,e,i,\\Omega,\\omega,f]` (True anomaly)")
        .value("MEE", kep3::MEE, "Modified Equinoctial Elements :math:`[p,f,g,h,k,L]` (Mean Longitude)")
        .value("MEE_R", kep3::MEE_R,
               "Modified Equinoctial Elements (retrograde) :math:`[p,f,g,h,k,L]` (Mean "
               "Longitude)")
        .value("POSVEL", kep3::POSVEL, "Position and Velocity")
        .export_values();

    py::enum_<kep3::optimality_type>(m, "optimality_type", "")
        .value("MASS", kep3::optimality_type::MASS, "Mass optimality")
        .value("TIME", kep3::optimality_type::TIME, "Time optimality")
        .export_values();

    // We expose the various anomaly conversions
    m.def("m2e", &kep3::m2e, pk::m2e_doc().c_str());
    m.def("e2m", &kep3::e2m, pk::e2m_doc().c_str());
    m.def("m2f", &kep3::m2f, pk::m2f_doc().c_str());
    m.def("f2m", &kep3::f2m, pk::f2m_doc().c_str());
    m.def("e2f", &kep3::e2f, pk::e2f_doc().c_str());
    m.def("f2e", &kep3::f2e, pk::f2e_doc().c_str());
    m.def("n2h", &kep3::n2h, pk::n2h_doc().c_str());
    m.def("h2n", &kep3::h2n, pk::h2n_doc().c_str());
    m.def("n2f", &kep3::n2f, pk::n2f_doc().c_str());
    m.def("f2n", &kep3::f2n, pk::f2n_doc().c_str());
    m.def("h2f", &kep3::h2f, pk::h2f_doc().c_str());
    m.def("f2h", &kep3::f2h, pk::f2h_doc().c_str());
    m.def("zeta2f", &kep3::zeta2f, pk::zeta2f_doc().c_str());
    m.def("f2zeta", &kep3::f2zeta, pk::f2zeta_doc().c_str());

    // And their vectorized versions
    m.def("m2e_v", py::vectorize(kep3::m2e), pk::m2e_v_doc().c_str());
    m.def("e2m_v", py::vectorize(kep3::e2m), pk::e2m_v_doc().c_str());
    m.def("m2f_v", py::vectorize(kep3::m2f), pk::m2f_v_doc().c_str());
    m.def("f2m_v", py::vectorize(kep3::f2m), pk::f2m_v_doc().c_str());
    m.def("e2f_v", py::vectorize(kep3::e2f), pk::e2f_v_doc().c_str());
    m.def("f2e_v", py::vectorize(kep3::f2e), pk::f2e_v_doc().c_str());
    m.def("n2h_v", py::vectorize(kep3::n2h), pk::n2h_v_doc().c_str());
    m.def("h2n_v", py::vectorize(kep3::h2n), pk::h2n_v_doc().c_str());
    m.def("n2f_v", py::vectorize(kep3::n2f), pk::n2f_v_doc().c_str());
    m.def("f2n_v", py::vectorize(kep3::f2n), pk::f2n_v_doc().c_str());
    m.def("h2f_v", py::vectorize(kep3::h2f), pk::h2f_v_doc().c_str());
    m.def("f2h_v", py::vectorize(kep3::f2h), pk::f2h_v_doc().c_str());
    m.def("zeta2f_v", py::vectorize(kep3::zeta2f), pk::zeta2f_v_doc().c_str());
    m.def("f2zeta_v", py::vectorize(kep3::f2zeta), pk::f2zeta_v_doc().c_str());

    // Exposing element conversions
    m.def("ic2par", &kep3::ic2par, py::arg("posvel"), py::arg("mu"), pk::ic2par_doc().c_str());
    m.def("par2ic", &kep3::par2ic, py::arg("elem"), py::arg("mu"), pk::par2ic_doc().c_str());
    m.def("ic2mee", &kep3::ic2mee, py::arg("posvel"), py::arg("mu"), py::arg("retrogde") = false,
          pk::ic2mee_doc().c_str());
    m.def("mee2ic", &kep3::mee2ic, py::arg("eq_elem"), py::arg("mu"), py::arg("retrogde") = false,
          pk::mee2ic_doc().c_str());
    m.def("par2mee", &kep3::par2mee, py::arg("elem"), py::arg("retrogde") = false, pk::par2mee_doc().c_str());
    m.def("mee2par", &kep3::mee2par, py::arg("eq_elem"), py::arg("retrogde") = false, pk::mee2par_doc().c_str());

    // Exposing mima functions and basic transfer functionalities
    m.def("mima", &kep3::mima, py::arg("dv1"), py::arg("dv2"), py::arg("tof"), py::arg("Tmax"), py::arg("veff"),
          pk::mima_doc().c_str());
    m.def("mima_from_hop", &kep3::mima_from_hop, py::arg("pl_s"), py::arg("pl_f"), py::arg("when_s"), py::arg("when_f"),
          py::arg("Tmax"), py::arg("veff"), pk::mima_from_hop_doc().c_str());
    m.def("mima2", &kep3::mima2, py::arg("posvel1"), py::arg("dv1"), py::arg("dv2"), py::arg("tof"), py::arg("Tmax"),
          py::arg("veff"), py::arg("mu"), pk::mima2_doc().c_str());
    m.def("mima2_from_hop", &kep3::mima2_from_hop, py::arg("pl_s"), py::arg("pl_f"), py::arg("when_s"),
          py::arg("when_f"), py::arg("Tmax"), py::arg("veff"), pk::mima2_from_hop_doc().c_str());
    m.def("hohmann", &kep3::hohmann, py::arg("r1"), py::arg("r2"), py::arg("mu"), pk::hohmann_doc().c_str());
    m.def("bielliptic", &kep3::bielliptic, py::arg("r1"), py::arg("r2"), py::arg("rb"), py::arg("mu"),
          pk::bielliptic_doc().c_str());

    // Exposing encoding conversions
    m.def("alpha2direct", &kep3::alpha2direct, py::arg("alphas"), py::arg("tof"), pk::alpha2direct_doc().c_str());
    m.def("direct2alpha", &kep3::direct2alpha, py::arg("tofs"), pk::direct2alpha_doc().c_str());
    m.def("eta2direct", &kep3::eta2direct, py::arg("etas"), py::arg("max_tof"), pk::eta2direct_doc().c_str());
    m.def("direct2eta", &kep3::direct2eta, py::arg("tofs"), py::arg("max_tof"), pk::direct2eta_doc().c_str());

    // Class epoch
    py::class_<kep3::epoch> epoch_class(m, "epoch", "Represents a specific point in time.");

    py::enum_<kep3::epoch::julian_type>(epoch_class, "julian_type")
        .value("MJD2000", kep3::epoch::julian_type::MJD2000, "Modified Julian Date 2000.")
        .value("MJD", kep3::epoch::julian_type::MJD, "Modified Julian Date.")
        .value("JD", kep3::epoch::julian_type::JD, "Julian Date.");

    py::enum_<kep3::epoch::string_format>(epoch_class, "string_format")
        .value("ISO", kep3::epoch::string_format::ISO, "ISO 8601 format for dates.");

    // This must go after the enum class registration
    epoch_class
        // Construtor from julian floats/int
        .def(py::init<double, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000, pk::epoch_from_float_doc().c_str())
        .def(py::init<int, kep3::epoch::julian_type>(), py::arg("when"),
             py::arg("julian_type") = kep3::epoch::julian_type::MJD2000)
        // Constructor from string
        .def(py::init<std::string, kep3::epoch::string_format>(), py::arg("when"),
             py::arg("string_format") = kep3::epoch::string_format::ISO, pk::epoch_from_string_doc().c_str())
        // Constructor from datetime py::object
        .def(py::init([](const py::object &in) {
                 // We check that `in` is a datetimeobject
                 const py::object Datetime = py::module_::import("datetime").attr("datetime");
                 if (!py::isinstance(in, Datetime)) {
                     pykep::py_throw(PyExc_TypeError, ("it seems you are trying to construct kep3::epoch object from a "
                                                       "python object that is not of type datetime"));
                 }
                 // We collect its info
                 const int y = in.attr("year").cast<int>();
                 auto m = in.attr("month").cast<unsigned>();
                 auto d = in.attr("day").cast<unsigned>();
                 const int h = in.attr("hour").cast<int>();
                 const int min = in.attr("minute").cast<int>();
                 const int s = in.attr("second").cast<int>();
                 const int us = in.attr("microsecond").cast<int>();
                 return kep3::epoch(y, m, d, h, min, s, 0, us);
             }),
             py::arg("when"), pk::epoch_from_datetime_doc().c_str())
        // repr()
        .def("__repr__", &pykep::ostream_repr<kep3::epoch>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::epoch>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::epoch>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::epoch>, &pykep::pickle_setstate_wrapper<kep3::epoch>))
        // now().
        .def_static("now", &kep3::epoch::now, "Returns a pykep.epoch with the current UTC date.")
        // julian dates
        .def_property_readonly("mjd2000", &kep3::epoch::mjd2000, "The Modified Julian Date 2000")
        .def_property_readonly("mjd", &kep3::epoch::mjd, "The Modified Julian Date")
        .def_property_readonly("jd", &kep3::epoch::jd, "The Julian Date")
        // comparison operators
        .def("__lt__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 < ep2; })
        .def("__gt__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 > ep2; })
        .def("__le__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 <= ep2; })
        .def("__ge__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 >= ep2; })
        .def("__eq__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 == ep2; })
        .def("__ne__", [](const kep3::epoch &ep1, const kep3::epoch &ep2) { return ep1 != ep2; })
        // math
        .def("__add__",
             [](kep3::epoch ep, double dt) { return ep + std::chrono::duration<double, std::ratio<86400>>(dt); })
        .def("__add__", [](kep3::epoch ep, std::chrono::duration<double, std::ratio<1>> dt) { return ep + dt; })
        .def("__sub__",
             [](kep3::epoch ep, double dt) { return ep - std::chrono::duration<double, std::ratio<86400>>(dt); })
        .def("__sub__", [](kep3::epoch ep, std::chrono::duration<double, std::ratio<1>> dt) { return ep - dt; });

    // Class planet (type erasure machinery here)
    py::class_<kep3::planet> planet_class(m, "planet", py::dynamic_attr{}, pykep::planet_docstring().c_str());
    // Expose extract.
    planet_class.def("_cpp_extract", &pykep::generic_cpp_extract<kep3::planet, kep3::udpla::keplerian>,
                     py::return_value_policy::reference_internal);
    // repr().
    planet_class.def("__repr__", [](const kep3::planet &pl) { return pl.get_name(); });
    // Full info
    planet_class.def("info", &pykep::ostream_repr<kep3::planet>);
    // Copy and deepcopy.
    planet_class.def("__copy__", &pykep::generic_copy_wrapper<kep3::planet>);
    planet_class.def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::planet>);
    // UDPLA extraction for python stuff.
    planet_class.def("_py_extract", &pykep::generic_py_extract<kep3::planet>);
    // Pickle support.
    planet_class.def(
        py::pickle(&pykep::pickle_getstate_wrapper<kep3::planet>, &pykep::pickle_setstate_wrapper<kep3::planet>));
    // Planet methods.
    planet_class.def(
        "eph",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when) {
            return std::visit([&](const auto &v) { return pl.eph(v); }, when);
        },
        py::arg("when"), pykep::planet_eph_docstring().c_str());
    planet_class.def(
        "acc",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when) {
            return std::visit([&](const auto &v) { return pl.acc(v); }, when);
        },
        py::arg("when"), pykep::planet_acc_docstring().c_str());
    // Vectorized versions. Note that the udpla method flattens everything but planet returns a non flat array.
    planet_class.def(
        "eph_v",
        [](const kep3::planet &pl, const std::vector<double> &eps) {
            std::vector<double> res = pl.eph_v(eps);
            // We create a capsule for the py::array_t to manage ownership change.
            auto vec_ptr = std::make_unique<std::vector<double>>(std::move(res));

            py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                const std::unique_ptr<std::vector<double>> vptr(static_cast<std::vector<double> *>(ptr));
            });

            // NOTE: at this point, the capsule has been created successfully (including
            // the registration of the destructor). We can thus release ownership from vec_ptr,
            // as now the capsule is responsible for destroying its contents. If the capsule constructor
            // throws, the destructor function is not registered/invoked, and the destructor
            // of vec_ptr will take care of cleaning up.
            auto *ptr = vec_ptr.release();

            return py::array_t<double>(py::array::ShapeContainer{boost::numeric_cast<py::ssize_t>(eps.size()),
                                                                 static_cast<py::ssize_t>(6)}, // shape
                                       ptr->data(), std::move(vec_caps));
        },
        py::arg("when"), pykep::planet_eph_v_docstring().c_str());

    planet_class.def(
        "acc_v",
        [](const kep3::planet &pl, const std::vector<double> &mjd2000s) {
            std::vector<double> res = pl.acc_v(mjd2000s);
            // We create a capsule for the py::array_t to manage ownership change.
            auto vec_ptr = std::make_unique<std::vector<double>>(std::move(res));

            py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                const std::unique_ptr<std::vector<double>> vptr(static_cast<std::vector<double> *>(ptr));
            });

            // NOTE: at this point, the capsule has been created successfully (including
            // the registration of the destructor). We can thus release ownership from vec_ptr,
            // as now the capsule is responsible for destroying its contents. If the capsule constructor
            // throws, the destructor function is not registered/invoked, and the destructor
            // of vec_ptr will take care of cleaning up.
            auto *ptr = vec_ptr.release();

            return py::array_t<double>(py::array::ShapeContainer{boost::numeric_cast<py::ssize_t>(mjd2000s.size()),
                                                                 static_cast<py::ssize_t>(3)}, // shape
                                       ptr->data(), std::move(vec_caps));
        },
        py::arg("when"), pykep::planet_acc_v_docstring().c_str());

#define PYKEP3_EXPOSE_PLANET_GETTER(name)                                                                              \
    planet_class.def(                                                                                                  \
        "get_" #name, [](const kep3::planet &pl) { return pl.get_##name(); },                                          \
        pykep::planet_get_##name##_docstring().c_str());

    PYKEP3_EXPOSE_PLANET_GETTER(name);
    PYKEP3_EXPOSE_PLANET_GETTER(extra_info);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_central_body);
    PYKEP3_EXPOSE_PLANET_GETTER(mu_self);
    PYKEP3_EXPOSE_PLANET_GETTER(radius);
    PYKEP3_EXPOSE_PLANET_GETTER(safe_radius);

#undef PYKEP3_EXPOSE_PLANET_GETTER

// We also support the various quantities as attributes for compatibility with pykep 2
// and because its nicer syntax to have them as attributes.
#define PYKEP3_EXPOSE_PLANET_ATTRIBUTE(name)                                                                           \
    planet_class.def_property_readonly(                                                                                \
        #name, [](const kep3::planet &pl) { return pl.get_##name(); },                                                 \
        pykep::planet_get_##name##_docstring().c_str());

    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(name);
    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(extra_info);
    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(mu_central_body);
    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(mu_self);
    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(radius);
    PYKEP3_EXPOSE_PLANET_ATTRIBUTE(safe_radius);

#undef PYKEP3_EXPOSE_PLANET_GETTER

    planet_class.def(
        "period",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when) {
            return std::visit([&](const auto &v) { return pl.period(v); }, when);
        },
        py::arg("when") = 0., pykep::planet_period_docstring().c_str());

    planet_class.def(
        "elements",
        [](const kep3::planet &pl, const std::variant<double, kep3::epoch> &when, kep3::elements_type el_ty) {
            return std::visit([&](const auto &v) { return pl.elements(v, el_ty); }, when);
        },
        py::arg("when") = 0., py::arg("el_type") = kep3::elements_type::KEP_F,
        pykep::planet_elements_docstring().c_str());

    // We now expose the cpp udplas. They will also add a constructor and the extract machinery to the planet_class
    // UDPLA module
    pykep::expose_all_udplas(m, planet_class);

    // Finalize (this constructor must be the last one of planet_class: else overload will fail with all the others)
    planet_class.def(py::init([](const py::object &o) { return kep3::planet{pk::python_udpla(o)}; }), py::arg("udpla"));

    // Exposing the Lambert problem class
    py::class_<kep3::lambert_problem> lambert_problem(m, "lambert_problem", pykep::lambert_problem_docstring().c_str());
    lambert_problem
        .def(py::init<const std::array<double, 3> &, const std::array<double, 3> &, double, double, bool, unsigned>(),
             py::arg("r0") = std::array<double, 3>{{1., 0., 0}}, py::arg("r1") = std::array<double, 3>{{0., 1., 0}},
             py::arg("tof") = kep3::pi / 2, py::arg("mu") = 1., py::arg("cw") = false, py::arg("multi_revs") = 1)
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::lambert_problem>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::lambert_problem>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::lambert_problem>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::lambert_problem>,
                        &pykep::pickle_setstate_wrapper<kep3::lambert_problem>))
        .def_property_readonly("v0", &kep3::lambert_problem::get_v0, "The velocity at the first point.")
        .def_property_readonly("v1", &kep3::lambert_problem::get_v1, "The velocity at the second point.")
        .def_property_readonly("r0", &kep3::lambert_problem::get_r0, "The first point.")
        .def_property_readonly("r1", &kep3::lambert_problem::get_r1, "The second point.")
        .def_property_readonly("tof", &kep3::lambert_problem::get_tof, "The time of flight between the two points.")
        .def_property_readonly("mu", &kep3::lambert_problem::get_mu,
                               "The gravitational parameter of the attracting body.")
        .def_property_readonly("x", &kep3::lambert_problem::get_x,
                               "The Battin variable x along the time of flight curves.")
        .def_property_readonly("iters", &kep3::lambert_problem::get_iters, "The number of iterations made.")
        .def_property_readonly("Nmax", &kep3::lambert_problem::get_Nmax, "The maximum number of iterations allowed.")
        .def_property_readonly("cw", &kep3::lambert_problem::get_cw, "The clockwise parameter.");

    // Exposing Taylor adaptive propagators
    // Create submodule "ta_cxx"
    py::module_ ta = m.def_submodule("ta_cxx", "Submodule for heyoka Taylor integrator related stuff");
    // Register the submodule so Python sees it as real (note we use the final python visible name for this)
    py::module_ sys = py::module_::import("sys");
    sys.attr("modules")["pykep.ta_cxx"] = ta;

    // KEP
    ta.def(
        "get_kep",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_kep(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_kep_docstring().c_str());
    ta.def(
        "get_kep_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_kep_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_kep_var_docstring().c_str());
    ta.def("kep_dyn", &kep3::ta::kep_dyn, pykep::kep_dyn_docstring().c_str());

    // BCP
    ta.def(
        "get_bcp",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_bcp(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_bcp_docstring().c_str());
    ta.def(
        "get_bcp_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_bcp_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_bcp_var_docstring().c_str());
    ta.def("bcp_dyn", &kep3::ta::bcp_dyn, pykep::bcp_dyn_docstring().c_str());

    // CR3BP
    // Add function to submodule
    ta.def("cr3bp_jacobi_C", &kep3::ta::cr3bp_jacobi_C, pykep::cr3bp_jacobi_C_docstring().c_str());
    ta.def("cr3bp_effective_potential_U", &kep3::ta::cr3bp_effective_potential_U,
           pykep::cr3bp_effective_potential_U_docstring().c_str());

    ta.def(
        "get_cr3bp",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_cr3bp(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_cr3bp_docstring().c_str());
    ta.def(
        "get_cr3bp_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_cr3bp_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_cr3bp_var_docstring().c_str());
    ta.def("cr3bp_dyn", &kep3::ta::cr3bp_dyn, pykep::cr3bp_dyn_docstring().c_str());
    // BCP
    ta.def(
        "get_bcp",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_bcp(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_bcp_docstring().c_str());
    ta.def(
        "get_bcp_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_bcp_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_bcp_var_docstring().c_str());
    ta.def("bcp_dyn", &kep3::ta::bcp_dyn, pykep::bcp_dyn_docstring().c_str());
    // Pontryagin Cartesian
    ta.def(
        "get_pc",
        [](double tol, kep3::optimality_type optimality) {
            auto ta_cache = kep3::ta::get_ta_pc(tol, optimality);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), py::arg("optimality"), pykep::get_pc_docstring().c_str());
    ta.def(
        "get_pc_var",
        [](double tol, kep3::optimality_type optimality) {
            auto ta_cache = kep3::ta::get_ta_pc_var(tol, optimality);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), py::arg("optimality"), pykep::get_pc_var_docstring().c_str());
    ta.def("pc_dyn", &kep3::ta::pc_dyn, pykep::pc_dyn_docstring().c_str());
    ta.def("get_pc_H_cfunc", &kep3::ta::get_pc_H_cfunc, pykep::get_pc_H_cfunc_docstring().c_str());
    ta.def("get_pc_SF_cfunc", &kep3::ta::get_pc_SF_cfunc, pykep::get_pc_SF_cfunc_docstring().c_str());
    ta.def("get_pc_u_cfunc", &kep3::ta::get_pc_u_cfunc, pykep::get_pc_u_cfunc_docstring().c_str());
    ta.def("get_pc_i_vers_cfunc", &kep3::ta::get_pc_i_vers_cfunc, pykep::get_pc_i_vers_cfunc_docstring().c_str());
    ta.def("get_pc_dyn_cfunc", &kep3::ta::get_pc_dyn_cfunc, pykep::get_pc_dyn_cfunc_docstring().c_str());

    // zero_hold_kep
    ta.def(
        "get_zero_hold_kep",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_zero_hold_kep(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_zero_hold_kep_docstring().c_str());
    ta.def(
        "get_zero_hold_kep_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_zero_hold_kep_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_zero_hold_kep_var_docstring().c_str());
    ta.def("zero_hold_kep_dyn", &kep3::ta::zero_hold_kep_dyn, pykep::zero_hold_kep_dyn_docstring().c_str());

    // zero_hold_cr3bp
    ta.def(
        "get_zero_hold_cr3bp",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_zero_hold_cr3bp(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_zero_hold_cr3bp_docstring().c_str());
    ta.def(
        "get_zero_hold_cr3bp_var",
        [](double tol) {
            auto ta_cache = kep3::ta::get_ta_zero_hold_cr3bp_var(tol);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), pykep::get_zero_hold_cr3bp_var_docstring().c_str());
    ta.def("zero_hold_cr3bp_dyn", &kep3::ta::zero_hold_cr3bp_dyn, pykep::zero_hold_cr3bp_dyn_docstring().c_str());

    // Pontryagin Equinoctial (TPBVP)
    ta.def(
        "get_peq",
        [](double tol, kep3::optimality_type optimality) {
            // retreive from cache
            auto ta_cache = kep3::ta::get_ta_peq(tol, optimality);
            // copy
            heyoka::taylor_adaptive<double> ta(ta_cache);
            // return a copy
            return ta;
        },
        py::arg("tol"), py::arg("optimality"), pykep::get_peq_docstring().c_str());
    ta.def(
        "get_peq_var",
        [](double tol, kep3::optimality_type optimality) {
            auto ta_cache = kep3::ta::get_ta_peq_var(tol, optimality);
            heyoka::taylor_adaptive<double> ta(ta_cache);
            return ta;
        },
        py::arg("tol"), py::arg("optimality"), pykep::get_peq_var_docstring().c_str());
    ta.def("peq_dyn", &kep3::ta::peq_dyn, pykep::peq_dyn_docstring().c_str());
    ta.def("get_peq_H_cfunc", &kep3::ta::get_peq_H_cfunc, pykep::get_peq_H_cfunc_docstring().c_str());
    ta.def("get_peq_SF_cfunc", &kep3::ta::get_peq_SF_cfunc, pykep::get_peq_SF_cfunc_docstring().c_str());
    ta.def("get_peq_u_cfunc", &kep3::ta::get_peq_u_cfunc, pykep::get_peq_u_cfunc_docstring().c_str());
    ta.def("get_peq_i_vers_cfunc", &kep3::ta::get_peq_i_vers_cfunc, pykep::get_peq_i_vers_cfunc_docstring().c_str());
    ta.def("get_peq_dyn_cfunc", &kep3::ta::get_peq_dyn_cfunc, pykep::get_peq_dyn_cfunc_docstring().c_str());

    // Exposing propagators
    m.def(
        "propagate_lagrangian",
        [](const std::array<std::array<double, 3>, 2> &pos_vel, double dt, double mu, bool request_stm) -> py::object {
        auto pl_retval = kep3::propagate_lagrangian(pos_vel, dt, mu, request_stm);
            if (pl_retval.second) {
                // The stm was requested lets transfer ownership to python
                const std::array<double, 36> &stm = pl_retval.second.value();

                // We create a capsule for the py::array_t to manage ownership change.
                auto vec_ptr = std::make_unique<std::array<double, 36>>(stm);

                py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                    const std::unique_ptr<std::array<double, 36>> vptr(static_cast<std::array<double, 36> *>(ptr));
                });

                // NOTE: at this point, the capsule has been created successfully (including
                // the registration of the destructor). We can thus release ownership from vec_ptr,
                // as now the capsule is responsible for destroying its contents. If the capsule constructor
                // throws, the destructor function is not registered/invoked, and the destructor
                // of vec_ptr will take care of cleaning up.
                auto *ptr = vec_ptr.release();

                auto computed_stm = py::array_t<double>(
                    py::array::ShapeContainer{static_cast<py::ssize_t>(6), static_cast<py::ssize_t>(6)}, // shape
                    ptr->data(), std::move(vec_caps));
                return py::make_tuple(py::make_tuple(pl_retval.first[0], pl_retval.first[1]), computed_stm);
            } else {
                return py::make_tuple(pl_retval.first[0], pl_retval.first[1]);
            }
        },
        py::arg("rv") = std::array<std::array<double, 3>, 2>{{{1, 0, 0}, {0, 1, 0}}}, py::arg("tof") = kep3::pi / 2,
        py::arg("mu") = 1, py::arg("stm") = false, pykep::propagate_lagrangian_docstring().c_str());

    m.def("propagate_lagrangian_grid", [](const std::array<std::array<double, 3>, 2> &pos_vel, const std::vector<double> &time_grid, double mu,
        bool request_stm){
            auto retval_c = kep3::propagate_lagrangian_grid(pos_vel, time_grid, mu, request_stm);      
            std::vector<py::tuple> retval_py{};
            retval_py.reserve(time_grid.size());
            for (decltype(time_grid.size()) i = 0u; i < time_grid.size(); ++i) {
                if (request_stm) {
                    // The stm was requested lets transfer ownership to python
                    const std::array<double, 36> &stm = retval_c[i].second.value();

                    // We create a capsule for the py::array_t to manage ownership change.
                    auto vec_ptr = std::make_unique<std::array<double, 36>>(stm);

                    py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                        const std::unique_ptr<std::array<double, 36>> vptr(static_cast<std::array<double, 36> *>(ptr));
                    });

                    // NOTE: at this point, the capsule has been created successfully (including
                    // the registration of the destructor). We can thus release ownership from vec_ptr,
                    // as now the capsule is responsible for destroying its contents. If the capsule constructor
                    // throws, the destructor function is not registered/invoked, and the destructor
                    // of vec_ptr will take care of cleaning up.
                    auto *ptr = vec_ptr.release();

                    auto computed_stm = py::array_t<double>(
                        py::array::ShapeContainer{static_cast<py::ssize_t>(6), static_cast<py::ssize_t>(6)}, // shape
                        ptr->data(), std::move(vec_caps));
                    retval_py.push_back(py::make_tuple(retval_c[i].first, computed_stm));
                } else {
                    retval_py.push_back(py::make_tuple(retval_c[i].first[0], retval_c[i].first[1]));
                }
            }
        return retval_py;
        }
        , py::arg("rv") = std::array<std::array<double, 3>, 2>{{{1, 0, 0}, {0, 1, 0}}}, py::arg("tofs") = std::vector<double>{kep3::pi / 2,},
        py::arg("mu") = 1, py::arg("stm") = false, pykep::propagate_lagrangian_grid_docstring().c_str());

    // Exposing the zero_hold_kep problem class
    py::class_<kep3::zero_hold_kep_problem> zero_hold_kep_problem(m, "zero_hold_kep_problem",
                                                                  pykep::zero_hold_kep_problem_docstring().c_str());
    zero_hold_kep_problem
        .def(py::init<double, double, double>(), py::arg("mu") = 1., py::arg("veff") = 1., py::arg("tol") = 1e-16)
        // repr().
        .def("__repr__", &pykep::ostream_repr<kep3::zero_hold_kep_problem>)
        // Copy and deepcopy.
        .def("__copy__", &pykep::generic_copy_wrapper<kep3::zero_hold_kep_problem>)
        .def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::zero_hold_kep_problem>)
        // Pickle support.
        .def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::zero_hold_kep_problem>,
                        &pykep::pickle_setstate_wrapper<kep3::zero_hold_kep_problem>))
        .def_property_readonly("mu", &kep3::zero_hold_kep_problem::get_mu, "The central body gravity parameter.")
        .def_property_readonly("veff", &kep3::zero_hold_kep_problem::get_veff, "The effective velocity (Isp g0)")
        .def_property_readonly("tol", &kep3::zero_hold_kep_problem::get_tol, "The Taylor integrator tolerance.")
        // The actual call to propagators. (we do not here care about copies and allocations as this is 20 times slower
        // than propagate lagrangian already on c++ side).
        .def("propagate", &kep3::zero_hold_kep_problem::propagate, py::arg("rvm_state"), py::arg("thrust"),
             py::arg("tof"), pykep::zero_hold_kep_problem_propagate_docstring().c_str())
        .def(
            "propagate_var",
            [](kep3::zero_hold_kep_problem &sp, const std::array<double, 7> &rvm_state, std::array<double, 3> thrust,
               double tof) {
                auto sp_retval = sp.propagate_var(rvm_state, thrust, tof);
                // Lets transfer ownership of dxdx to python (not sure this is actually needed to
                // get an efficient return value ... maybe its overkill here). It surely avoid one more copy /
                // allocation of 49+21 values, but in the overall algorithm maybe irrelevant.
                const std::array<double, 49> &dxdx = std::get<1>(sp_retval);

                // We create a capsule for the py::array_t to manage ownership change.
                auto vec_ptr = std::make_unique<std::array<double, 49>>(dxdx);

                py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                    const std::unique_ptr<std::array<double, 49>> vptr(static_cast<std::array<double, 49> *>(ptr));
                });

                // NOTE: at this point, the capsule has been created successfully (including
                // the registration of the destructor). We can thus release ownership from vec_ptr,
                // as now the capsule is responsible for destroying its contents. If the capsule constructor
                // throws, the destructor function is not registered/invoked, and the destructor
                // of vec_ptr will take care of cleaning up.
                auto *ptr = vec_ptr.release();

                auto computed_dxdx = py::array_t<double>(
                    py::array::ShapeContainer{static_cast<py::ssize_t>(7), static_cast<py::ssize_t>(7)}, // shape
                    ptr->data(), std::move(vec_caps));

                // Lets transfer ownership of dxdu to python
                const std::array<double, 21> &dxdu = std::get<2>(sp_retval);

                // We create a capsule for the py::array_t to manage ownership change.
                auto vec_ptr2 = std::make_unique<std::array<double, 21>>(dxdu);

                py::capsule vec_caps2(vec_ptr2.get(), [](void *ptr) {
                    const std::unique_ptr<std::array<double, 21>> vec_ptr2(static_cast<std::array<double, 21> *>(ptr));
                });

                // NOTE: at this point, the capsule has been created successfully (including
                // the registration of the destructor). We can thus release ownership from vec_ptr,
                // as now the capsule is responsible for destroying its contents. If the capsule constructor
                // throws, the destructor function is not registered/invoked, and the destructor
                // of vec_ptr will take care of cleaning up.
                auto *ptr2 = vec_ptr2.release();

                auto computed_dxdu = py::array_t<double>(
                    py::array::ShapeContainer{static_cast<py::ssize_t>(7), static_cast<py::ssize_t>(3)}, // shape
                    ptr2->data(), std::move(vec_caps2));
                return py::make_tuple(std::get<0>(sp_retval), computed_dxdx, computed_dxdu);
            },
            py::arg("rvm_state"), py::arg("thrust"), py::arg("tof"),
            pykep::zero_hold_kep_problem_propagate_var_docstring().c_str());

    // Exposing fly-by routines
    m.def("fb_con",
          py::overload_cast<const std::array<double, 3> &, const std::array<double, 3> &, const kep3::planet &>(
              &kep3::fb_con),
          py::arg("v_rel_in"), py::arg("v_rel_out"), py::arg("planet"), pykep::fb_con_docstring().c_str());
    m.def(
        "fb_con",
        py::overload_cast<const std::array<double, 3> &, const std::array<double, 3> &, double, double>(&kep3::fb_con),
        py::arg("v_rel_in"), py::arg("v_rel_out"), py::arg("mu"), py::arg("safe_radius"));

    m.def("fb_dv",
          py::overload_cast<const std::array<double, 3> &, const std::array<double, 3> &, const kep3::planet &>(
              &kep3::fb_dv),
          py::arg("v_rel_in"), py::arg("v_rel_out"), py::arg("planet"), pykep::fb_dv_docstring().c_str());
    m.def("fb_dv",
          py::overload_cast<const std::array<double, 3> &, const std::array<double, 3> &, double, double>(&kep3::fb_dv),
          py::arg("v_rel_in"), py::arg("v_rel_out"), py::arg("mu"), py::arg("safe_radius"));
    m.def("fb_vout", &kep3::fb_vout, py::arg("v_in"), py::arg("v_pla"), py::arg("rp"), py::arg("beta"), py::arg("mu"),
          pykep::fb_vout_docstring().c_str());

    // Exposing the sims_flanagan leg
    py::class_<kep3::leg::sims_flanagan> sims_flanagan(m, "_sims_flanagan", pykep::leg_sf_docstring().c_str());
    sims_flanagan.def(
        py::init<const std::array<std::array<double, 3>, 2> &, double, std::vector<double>,
                 const std::array<std::array<double, 3>, 2> &, double, double, double, double, double, double>(),
        py::arg("rvs") = std::array<std::array<double, 3>, 2>{{{1., 0, 0.}, {0., 1., 0.}}}, py::arg("ms") = 1.,
        py::arg("throttles") = std::vector<double>{0, 0, 0, 0, 0, 0},
        py::arg("rvf") = std::array<std::array<double, 3>, 2>{{{0., 1., 0.}, {-1., 0., 0.}}}, py::arg("mf") = 1.,
        py::arg("tof") = kep3::pi / 2, py::arg("max_thrust") = 1., py::arg("veff") = 1., py::arg("mu") = 1,
        py::arg("cut") = 0.5);
    // repr().
    sims_flanagan.def("__repr__", &pykep::ostream_repr<kep3::leg::sims_flanagan>);
    // Copy and deepcopy.
    sims_flanagan.def("__copy__", &pykep::generic_copy_wrapper<kep3::leg::sims_flanagan>);
    sims_flanagan.def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::leg::sims_flanagan>);
    // Pickle support.
    sims_flanagan.def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::leg::sims_flanagan>,
                                 &pykep::pickle_setstate_wrapper<kep3::leg::sims_flanagan>));
    // The rest
    sims_flanagan.def_property(
        "throttles", &kep3::leg::sims_flanagan::get_throttles,
        [](kep3::leg::sims_flanagan &sf, const std::vector<double> &throttles) { return sf.set_throttles(throttles); },
        pykep::leg_sf_throttles_docstring().c_str());

#define PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(name)                                                                          \
    sims_flanagan.def_property(#name, &kep3::leg::sims_flanagan::get_##name, &kep3::leg::sims_flanagan::set_##name,    \
                               pykep::leg_sf_##name##_docstring().c_str());
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(rvs);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(ms);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(rvf);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(mf);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(tof);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(max_thrust);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(veff);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(mu);
    PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES(cut);

#undef PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES

    sims_flanagan.def("compute_mismatch_constraints", &kep3::leg::sims_flanagan::compute_mismatch_constraints,
                      pykep::leg_sf_mc_docstring().c_str());
    sims_flanagan.def("compute_throttle_constraints", &kep3::leg::sims_flanagan::compute_throttle_constraints,
                      pykep::leg_sf_tc_docstring().c_str());
    sims_flanagan.def(
        "compute_mc_grad",
        [](const kep3::leg::sims_flanagan &leg) {
            auto tc_cpp = leg.compute_mc_grad();
            // Lets transfer ownership to python of the three
            const std::array<double, 49> &rs_addr = std::get<0>(tc_cpp);
            const std::array<double, 49> &rf_addr = std::get<1>(tc_cpp);
            const std::vector<double> &th_addr = std::get<2>(tc_cpp);

            // We create three separate capsules for the py::array_t to manage ownership change.
            auto vec_ptr_rs = std::make_unique<std::array<double, 49>>(rs_addr);
            py::capsule vec_caps_rs(vec_ptr_rs.get(), [](void *ptr) {
                const std::unique_ptr<std::array<double, 49>> vptr(static_cast<std::array<double, 49> *>(ptr));
            });
            auto vec_ptr_rf = std::make_unique<std::array<double, 49>>(rf_addr);
            py::capsule vec_caps_rf(vec_ptr_rf.get(), [](void *ptr) {
                const std::unique_ptr<std::array<double, 49>> vptr(static_cast<std::array<double, 49> *>(ptr));
            });
            auto vec_ptr_th = std::make_unique<std::vector<double>>(th_addr);
            py::capsule vec_caps_th(vec_ptr_th.get(), [](void *ptr) {
                const std::unique_ptr<std::vector<double>> vptr(static_cast<std::vector<double> *>(ptr));
            });
            // NOTE: at this point, the capsules have been created successfully (including
            // the registration of the destructor). We can thus release ownership from vec_ptr_xx,
            // as now the capsules are responsible for destroying its contents.
            auto *ptr_rs = vec_ptr_rs.release();
            auto *ptr_rf = vec_ptr_rf.release();
            auto *ptr_th = vec_ptr_th.release();
            auto rs_python = py::array_t<double>(
                py::array::ShapeContainer{static_cast<py::ssize_t>(7), static_cast<py::ssize_t>(7)}, // shape
                ptr_rs->data(), std::move(vec_caps_rs));
            auto rf_python = py::array_t<double>(
                py::array::ShapeContainer{static_cast<py::ssize_t>(7), static_cast<py::ssize_t>(7)}, // shape
                ptr_rf->data(), std::move(vec_caps_rf));
            auto th_python = py::array_t<double>(
                py::array::ShapeContainer{static_cast<py::ssize_t>(7),
                                          static_cast<py::ssize_t>(leg.get_nseg() * 3 + 1u)}, // shape
                ptr_th->data(), std::move(vec_caps_th));
            return py::make_tuple(rs_python, rf_python, th_python);
        },
        pykep::leg_sf_mc_grad_docstring().c_str());
    sims_flanagan.def(
        "compute_tc_grad",
        [](const kep3::leg::sims_flanagan &leg) {
            const std::vector<double> tc_cpp = leg.compute_tc_grad();
            // Lets transfer ownership to python
            const std::vector<double> &tc_cpp_addr = tc_cpp;
            // We create a capsule for the py::array_t to manage ownership change.
            auto vec_ptr = std::make_unique<std::vector<double>>(tc_cpp_addr);
            py::capsule vec_caps(vec_ptr.get(), [](void *ptr) {
                const std::unique_ptr<std::vector<double>> vptr(static_cast<std::vector<double> *>(ptr));
            });
            // NOTE: at this point, the capsule has been created successfully (including
            // the registration of the destructor). We can thus release ownership from vec_ptr,
            // as now the capsule is responsible for destroying its contents. If the capsule constructor
            // throws, the destructor function is not registered/invoked, and the destructor
            // of vec_ptr will take care of cleaning up.
            auto *ptr = vec_ptr.release();

            auto tc_python
                = py::array_t<double>(py::array::ShapeContainer{static_cast<py::ssize_t>(leg.get_nseg()),
                                                                static_cast<py::ssize_t>(leg.get_nseg() * 3)}, // shape
                                      ptr->data(), std::move(vec_caps));
            return tc_python;
        },
        pykep::leg_sf_tc_grad_docstring().c_str());
    sims_flanagan.def_property_readonly("nseg", &kep3::leg::sims_flanagan::get_nseg,
                                        pykep::leg_sf_nseg_docstring().c_str());
    sims_flanagan.def_property_readonly("nseg_fwd", &kep3::leg::sims_flanagan::get_nseg_fwd,
                                        pykep::leg_sf_nseg_fwd_docstring().c_str());
    sims_flanagan.def_property_readonly("nseg_bck", &kep3::leg::sims_flanagan::get_nseg_bck,
                                        pykep::leg_sf_nseg_bck_docstring().c_str());

    // Exposing the sims_flanagan_alpha leg
    py::class_<kep3::leg::sims_flanagan_alpha> sims_flanagan_alpha(m, "_sims_flanagan_alpha",
                                                                   pykep::leg_sf_alpha_docstring().c_str());
    sims_flanagan_alpha.def(
        py::init<const std::array<std::array<double, 3>, 2> &, double, std::vector<double>, std::vector<double>,
                 const std::array<std::array<double, 3>, 2> &, double, double, double, double, double, double>(),
        py::arg("rvs") = std::array<std::array<double, 3>, 2>{{{1., 0, 0.}, {0., 1., 0.}}}, py::arg("ms") = 1.,
        py::arg("throttles") = std::vector<double>{0, 0, 0, 0, 0, 0}, py::arg("talphas") = std::vector<double>{0, 0},
        py::arg("rvf") = std::array<std::array<double, 3>, 2>{{{0., 1., 0.}, {-1., 0., 0.}}}, py::arg("mf") = 1.,
        py::arg("tof") = kep3::pi / 2, py::arg("max_thrust") = 1., py::arg("veff") = 1., py::arg("mu") = 1,
        py::arg("cut") = 0.5);
    // repr().
    sims_flanagan_alpha.def("__repr__", &pykep::ostream_repr<kep3::leg::sims_flanagan_alpha>);
    // Copy and deepcopy.
    sims_flanagan_alpha.def("__copy__", &pykep::generic_copy_wrapper<kep3::leg::sims_flanagan_alpha>);
    sims_flanagan_alpha.def("__deepcopy__", &pykep::generic_deepcopy_wrapper<kep3::leg::sims_flanagan_alpha>);
    // Pickle support.
    sims_flanagan_alpha.def(py::pickle(&pykep::pickle_getstate_wrapper<kep3::leg::sims_flanagan_alpha>,
                                       &pykep::pickle_setstate_wrapper<kep3::leg::sims_flanagan_alpha>));
    // The rest
    sims_flanagan_alpha.def_property(
        "throttles", &kep3::leg::sims_flanagan_alpha::get_throttles,
        [](kep3::leg::sims_flanagan_alpha &sf, const std::vector<double> &throttles) {
            return sf.set_throttles(throttles);
        },
        pykep::leg_sf_throttles_docstring().c_str());

    // The rest
    sims_flanagan_alpha.def_property(
        "talphas", &kep3::leg::sims_flanagan_alpha::get_talphas,
        [](kep3::leg::sims_flanagan_alpha &sf, const std::vector<double> &talphas) { return sf.set_talphas(talphas); },
        pykep::leg_sf_talphas_docstring().c_str());

#define PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(name)                                                                    \
    sims_flanagan_alpha.def_property(#name, &kep3::leg::sims_flanagan_alpha::get_##name,                               \
                                     &kep3::leg::sims_flanagan_alpha::set_##name,                                      \
                                     pykep::leg_sf_##name##_docstring().c_str());
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(rvs);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(ms);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(rvf);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(mf);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(tof);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(max_thrust);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(veff);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(mu);
    PYKEP3_EXPOSE_LEG_SF_ALPHA_ATTRIBUTES(cut);

#undef PYKEP3_EXPOSE_LEG_SF_ATTRIBUTES

    using tas_type = std::pair<const heyoka::taylor_adaptive<double> &, const heyoka::taylor_adaptive<double> &>;

    sims_flanagan_alpha.def("compute_mismatch_constraints",
                            &kep3::leg::sims_flanagan_alpha::compute_mismatch_constraints,
                            pykep::leg_sf_mc_docstring().c_str());
    sims_flanagan_alpha.def("compute_throttle_constraints",
                            &kep3::leg::sims_flanagan_alpha::compute_throttle_constraints,
                            pykep::leg_sf_tc_docstring().c_str());
    sims_flanagan_alpha.def_property_readonly("nseg", &kep3::leg::sims_flanagan_alpha::get_nseg,
                                              pykep::leg_sf_nseg_docstring().c_str());
    sims_flanagan_alpha.def_property_readonly("nseg_fwd", &kep3::leg::sims_flanagan_alpha::get_nseg_fwd,
                                              pykep::leg_sf_nseg_fwd_docstring().c_str());
    sims_flanagan_alpha.def_property_readonly("nseg_bck", &kep3::leg::sims_flanagan_alpha::get_nseg_bck,
                                              pykep::leg_sf_nseg_bck_docstring().c_str());
}