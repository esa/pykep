// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <array>
#include <string>
#include <vector>

#include <fmt/core.h>
#include <fmt/std.h>

#include <kep3/core_astro/constants.hpp>
#include <kep3/planet.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "common_utils.hpp"
#include "python_udpla.hpp"

namespace pykep
{

python_udpla::python_udpla() = default;

python_udpla::python_udpla(py::object obj) : m_obj(std::move(obj))
{
    // Forbid the use of a pykep.planet as a UDPLA.
    if (pykep::type(m_obj).equal(py::module::import("pykep").attr("planet"))) {
        pykep::py_throw(PyExc_TypeError,
                        ("a pykep.planet cannot be used as a UDPLA to construct another pykep.planet (if you need to copy a "
                         "planet please use the standard Python copy()/deepcopy() functions)"));
    }
    // Check that o is an instance of a class, and not just a type.
    check_not_type(m_obj, "planet");
    // Check the presence of the mandatory methods
    check_mandatory_method(m_obj, "eph", "planet");
};

// Mandatory methods
[[nodiscard]] std::array<std::array<double, 3>, 2> python_udpla::eph(double mjd2000) const
{
    auto udpla_eph = pykep::callable_attribute(m_obj, "eph");
    if (udpla_eph.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError, ("the eph() method has been invoked, but it is not implemented "
                                                    "in the user-defined Python planet '"
                                                    + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                                                    + "': the method is either not present or not callable")
                                                       .c_str());
    }
    return py::cast<std::array<std::array<double, 3>, 2>>(udpla_eph(mjd2000));
}

// Optional methods
[[nodiscard]] std::vector<double> python_udpla::eph_v(const std::vector<double> &mjd2000s) const
{
    auto udpla_eph_v = pykep::callable_attribute(m_obj, "eph_v");
    if (!udpla_eph_v.is_none()) {
        auto ret = py::cast<py::array_t<double>>(udpla_eph_v(mjd2000s));
        auto retval = py::cast<std::vector<double>>(ret.attr("flatten")());
        return retval;
    } else {
        return kep3::detail::default_eph_vectorization(this, mjd2000s);
    }
}

[[nodiscard]] std::array<double, 3> python_udpla::acc(double mjd2000) const
{
    auto udpla_acc = pykep::callable_attribute(m_obj, "acc");
    if (udpla_acc.is_none()) {
        pykep::py_throw(PyExc_NotImplementedError, ("the acc() method has been invoked, but it is not implemented "
                                                    "in the user-defined Python planet '"
                                                    + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                                                    + "': the method is either not present or not callable")
                                                       .c_str());
    }
    return py::cast<std::array<double, 3>>(udpla_acc(mjd2000));
}

[[nodiscard]] std::vector<double> python_udpla::acc_v(const std::vector<double> &mjd2000s) const
{
    auto udpla_acc_v = pykep::callable_attribute(m_obj, "acc_v");
    if (!udpla_acc_v.is_none()) {
        auto ret = py::cast<py::array_t<double>>(udpla_acc_v(mjd2000s));
        auto retval = py::cast<std::vector<double>>(ret.attr("flatten")());
        return retval;
    } else {
        return kep3::detail::default_acc_vectorization(this, mjd2000s);
    }
}

[[nodiscard]] std::string python_udpla::get_name() const
{
    return getter_wrapper<std::string>(m_obj, "get_name", pykep::str(pykep::type(m_obj)));
}
[[nodiscard]] std::string python_udpla::get_extra_info() const
{
    return getter_wrapper<std::string>(m_obj, "get_extra_info", "");
}
[[nodiscard]] double python_udpla::get_mu_central_body() const
{
    return getter_wrapper<double>(m_obj, "get_mu_central_body", -1);
}
[[nodiscard]] double python_udpla::get_mu_self() const
{
    return getter_wrapper<double>(m_obj, "get_mu_self", -1);
}
[[nodiscard]] double python_udpla::get_radius() const
{
    return getter_wrapper<double>(m_obj, "get_radius", -1);
}
[[nodiscard]] double python_udpla::get_safe_radius() const
{
    return getter_wrapper<double>(m_obj, "get_safe_radius", -1);
}
[[nodiscard]] double python_udpla::period(double mjd2000) const
{
    auto udpla_period = pykep::callable_attribute(m_obj, "period");
    auto udpla_get_mu_central_body = pykep::callable_attribute(m_obj, "get_mu_central_body");

    if (udpla_period.is_none()) {
        if (udpla_get_mu_central_body.is_none()) {
            pykep::py_throw(
                PyExc_NotImplementedError,
                ("the period() method has been invoked, but "
                 "in the user-defined Python planet '"
                 + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                 + "': the methods period() or get_mu_central_body() are either not present or not callable")
                    .c_str());
        } else {
            // If the user provides the central body parameter, then compute the
            // period from the energy at epoch
            auto [r, v] = eph(mjd2000);
            double mu = get_mu_central_body();
            double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            double en = v2 / 2. - mu / R;

            if (en > 0) {
                // If the energy is positive we have an hyperbolae and we return nan
                return std::numeric_limits<double>::quiet_NaN();
            } else {
                double a = -mu / 2. / en;
                return kep3::pi * 2. * std::sqrt(a * a * a / mu);
            }
        }
    }
    return py::cast<double>(udpla_period(mjd2000));
}

[[nodiscard]] std::array<double, 6> python_udpla::elements(double mjd2000, kep3::elements_type el_type) const
{
    auto udpla_elements = pykep::callable_attribute(m_obj, "elements");
    auto udpla_get_mu_central_body = pykep::callable_attribute(m_obj, "get_mu_central_body");
    // If the user provides an efficient way to compute the orbital elements, then use it.
    if (!udpla_elements.is_none()) {
        return py::cast<std::array<double, 6>>(udpla_elements(mjd2000, el_type));
    } else if (!udpla_get_mu_central_body.is_none()) {
        // If the user provides the central body parameter, then compute the
        // elements using posvel computed at ep and converted.
        auto pos_vel = eph(mjd2000);
        double mu = get_mu_central_body();
        return kep3::detail::elements_from_posvel(pos_vel, mu, el_type);
    } else {
        pykep::py_throw(PyExc_NotImplementedError,
                        ("the elements() method has been invoked, but "
                         "in the user-defined Python planet '"
                         + pykep::str(m_obj) + "' of type '" + pykep::str(pykep::type(m_obj))
                         + "': the methods elements() or get_mu_central_body() are either not present or not callable")
                            .c_str());
    }
}

} // namespace pykep
