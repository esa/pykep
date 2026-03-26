// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#ifndef PYKEP_COMMON_UTILS_HPP
#define PYKEP_COMMON_UTILS_HPP

#include <sstream>
#include <string>

#include <kep3/detail/s11n.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "kep3/planet.hpp"
#include "python_udpla.hpp"

namespace pykep
{

namespace py = pybind11;

// Import and return the builtins module.
py::object builtins();

// Throw a Python exception.
[[noreturn]] void py_throw(PyObject *, const char *);

// Check if 'o' has a callable attribute (i.e., a method) named 's'. If so, it will
// return the attribute, otherwise it will return None.
py::object callable_attribute(const py::object &, const char *);

// Check if type is callable.
bool callable(const py::object &);

// Get the string representation of an object.
std::string str(const py::object &);

// Get the type of an object.
py::object type(const py::object &);

// Perform a deep copy of input object o.
py::object deepcopy(const py::object &);

// repr() via ostream.
template <typename T>
inline std::string ostream_repr(const T &x)
{
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

// Generic copy wrappers.
template <typename T>
inline T generic_copy_wrapper(const T &x)
{
    return x;
}

template <typename T>
inline T generic_deepcopy_wrapper(const T &x, const py::dict &)
{
    return x;
}

// Generic extract() wrappers.
template <typename C, typename T>
inline T *generic_cpp_extract(C &c, const T &)
{
    return c.template extract<T>();
}

template <typename C>
inline py::object generic_py_extract(C &c)
{
    auto ptr = c.template extract<pykep::python_udpla>();
    if (ptr) {
        return ptr->m_obj;
    }

    // The user-defined entity is not pythonic. Return None.
    return py::none();
}

// Check if a mandatory method is present in a user-defined entity.
void check_mandatory_method(const py::object &o, const char *s, const char *target);

// Check if the user is trying to construct a pykep object from a type, rather than from an object.
// This is an easy error to commit, and it is sneaky because the callable_attribute() machinery in
// check_mandatory_method() will detect the methods of the *class* (rather than instance methods), and it will thus not
// error out.
void check_not_type(const py::object &o, const char *target);

// A simple wrapper for getters. It will try to:
// - get the attribute "name" from the object o,
// - call it without arguments,
// - extract an instance from the ret value and return it.
// If the attribute is not there or it is not callable, the value "def_value" will be returned.
template <typename RetType>
static RetType getter_wrapper(const py::object &o, const char *name, const RetType &def_value)
{
    auto a = callable_attribute(o, name);
    if (a.is_none()) {
        return def_value;
    }
    return py::cast<RetType>(a());
}

// Helpers to implement pickling on top of Boost.Serialization.
template <typename T>
inline py::tuple pickle_getstate_wrapper(const T &x)
{
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oa(oss);
        oa << x;
    }

    return py::make_tuple(py::bytes(oss.str()));
}

template <typename T>
inline T pickle_setstate_wrapper(const py::tuple &state)
{
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("The state tuple passed to the deserialization wrapper "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto *ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "A bytes object is needed in the deserialization wrapper");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    T x;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> x;
    }

    return x;
}

// Specialize for kep3::planet to handle python_udpla stored inside
template <>
inline py::tuple pickle_getstate_wrapper<kep3::planet>(const kep3::planet &pl)
{
    // Use extract to get python_udpla pointer if it is stored
    auto *p0 = pl.extract<pykep::python_udpla>();

    if (p0) {
        return py::make_tuple(p0->m_obj);
    }

    // Else fallback to boost serialization
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oa(oss);
        oa << pl;
    }
    return py::make_tuple(py::bytes(oss.str()));
}

// Specialize for kep3::planet
template <>
inline kep3::planet pickle_setstate_wrapper<kep3::planet>(const py::tuple &state)
{
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, "Invalid state tuple size for kep3::planet");
    }

    py::handle st = state[0];

    if (PyBytes_Check(st.ptr())) {
        // This is a C++-backed planet, use Boost deserialization
        auto *ptr = PyBytes_AsString(st.ptr());
        if (!ptr) {
            py_throw(PyExc_TypeError, "Deserialization requires bytes");
        }
        std::istringstream iss(std::string(ptr, ptr + py::len(st)));
        kep3::planet pl;
        {
            boost::archive::binary_iarchive ia(iss);
            ia >> pl;
        }
        return pl;
    } else {
        // This is a Python-backed planet, reconstruct via python_udpla
        pykep::python_udpla pyudpla(st.cast<py::object>());
        return kep3::planet(pyudpla);
    }
    // Defensive: should not reach here
    throw std::runtime_error("Invalid pickle state for kep3::planet");
}

} // namespace pykep

#endif