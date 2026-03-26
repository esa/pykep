// Copyright © 2023–2025 Dario Izzo (dario.izzo@gmail.com), 
// Francesco Biscani (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// Licensed under the Mozilla Public License, version 2.0.
// You may obtain a copy of the MPL at https://www.mozilla.org/MPL/2.0/.

#include <cassert>
#include <initializer_list>
#include <iterator>
#include <string>
#include <tuple>
#include <utility> 
#include <vector>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "common_utils.hpp"

namespace pykep
{

namespace py = pybind11;

// Import and return the builtins module.
py::object builtins()
{
    return py::module::import("builtins");
}

// Throw a Python exception of type "type" with associated
// error message "msg".
void py_throw(PyObject *type, const char *msg)
{
    PyErr_SetString(type, msg);
    throw py::error_already_set();
}

// Check if 'o' has a callable attribute (i.e., a method) named 's'. If so, it will
// return the attribute, otherwise it will return None.
py::object callable_attribute(const py::object &o, const char *s)
{
    if (py::hasattr(o, s)) {
        auto retval = o.attr(s);
        if (callable(retval)) {
            return retval;
        }
    }
    return py::none();
}

// Check if type is callable.
bool callable(const py::object &o)
{
    return py::cast<bool>(builtins().attr("callable")(o));
}

// Get the string representation of an object.
std::string str(const py::object &o)
{
    return py::cast<std::string>(py::str(o));
}

// Get the type of an object.
py::object type(const py::object &o)
{
    return builtins().attr("type")(o);
}

// Perform a deep copy of input object o.
py::object deepcopy(const py::object &o)
{
    return py::module::import("copy").attr("deepcopy")(o);
}

// Check if a mandatory method is present in a user-defined entity.
void check_mandatory_method(const py::object &o, const char *s, const char *target)
{
    if (callable_attribute(o, s).is_none()) {
        py_throw(PyExc_NotImplementedError,
                 ("the mandatory '" + std::string(s) + "()' method has not been detected in the user-defined Python "
                  + std::string(target) + " '" + str(o) + "' of type '" + str(type(o))
                  + "': the method is either not present or not callable")
                     .c_str());
    }
}

// Check if the user is trying to construct a pagmo object from a type, rather than from an object.
// This is an easy error to commit, and it is sneaky because the callable_attribute() machinery will detect
// the methods of the *class* (rather than instance methods), and it will thus not error out.
void check_not_type(const py::object &o, const char *target)
{
    if (py::isinstance(o, builtins().attr("type"))) {
        py_throw(PyExc_TypeError, ("it seems like you are trying to instantiate a pykep " + std::string(target)
                                   + " using a type rather than an object instance: please construct an object "
                                     "and use that instead of the type in the "
                                   + std::string(target) + " constructor")
                                      .c_str());
    }
}

} // namespace pykep